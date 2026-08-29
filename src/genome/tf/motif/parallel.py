"""Scanning across processes: the shard plan, the pool, and the batches put back together.

**Processes and not threads.** The engine holds the GIL for the whole of a scan, so threads
would take turns rather than share the machine. Work is split by sequence across a process
pool, each worker building its own engine once — the setup a scan pays before it reads a
base — and keeping it for every shard it is handed.

**The start method is spawn, everywhere.** Fork is cheaper on Linux and would ship the
prepared matrices for nothing, but it is unsafe beside threads and is being retired; one
start method on every platform is worth more than the saving, when the worker then scans
for minutes.

**A long sequence is cut with an overlap, and each shard keeps only what it owns.** The
overlap is one less than the longest matrix, which is exactly enough for a hit *starting*
at the last position a shard owns to be scored whole; and each shard drops a hit starting
at or past the end of its own region, because the next shard owns that one. So a hit lying
across a boundary is reported once, by the shard whose region its start falls in — never
twice and never not at all. :func:`plan_shards` is that arithmetic on its own, so it can be
checked without starting anything.

**Serial and parallel produce the identical table**, and identical means row for row: the
shards of one sequence are put back into the order a serial scan would have emitted them —
motif by motif in the set's order, forward strand before reverse, ascending by position —
before the batch is yielded. That is the property the whole design is verified by.

Examples
--------
>>> plan_shards(10, overlap=2, shard_length=4)   # (offset, stop, owned)
[(0, 6, 4), (4, 10, 4), (8, 10, 2)]
>>> [stop - offset for offset, stop, _owned in plan_shards(10, 2, 4)]
[6, 6, 2]
"""

from __future__ import annotations

import multiprocessing
from collections import deque
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import Executor, Future, ProcessPoolExecutor
from dataclasses import dataclass
from itertools import islice
from typing import TYPE_CHECKING, Any

import numpy as np
import numpy.typing as npt
import pandas as pd

from genome.tf.motif.scan import _STRANDS, HIT_DTYPES, _batch, empty_hits, engine_for

if TYPE_CHECKING:  # pragma: no cover - typing only
    import MOODS.scan

    from genome.tf.motif.scan import _Prepared

#: The start method every pool here uses. See this module's own explanation for why it is
#: not the platform default.
START_METHOD = "spawn"

#: How long a piece of one sequence a worker is given at a time. Big enough that the
#: overlap and the round trip are noise beside the scan, small enough that a chromosome
#: becomes tens of pieces — which balances the workers and bounds what is in flight, since
#: every shard is copied to its worker.
_SHARD_LENGTH = 5_000_000

#: How many shards may be outstanding per worker. Two keeps every worker fed across the
#: hand-back of a result; more would only buy memory a claim on nothing.
_IN_FLIGHT = 2

#: The worker's own engine, built once by :func:`_install` and kept for every shard that
#: worker is given. A module global because that is the only thing a pool initializer can
#: hand to the tasks that follow it.
_ENGINE: tuple[_Prepared, MOODS.scan.Scanner] | None = None


@dataclass(frozen=True)
class _Shard:
    """One piece of one named sequence, and which part of it this piece is answerable for.

    ``sequence`` is the piece including its overlap; ``offset`` is where the piece starts in
    the whole sequence, so a position inside it reports as ``offset + position``; ``owned``
    is how many positions from ``offset`` this shard keeps hits for; and ``of`` is how many
    shards the sequence was cut into, which is how the parent knows when it has them all.
    """

    name: str
    sequence: str
    offset: int
    owned: int
    of: int


def plan_shards(length: int, overlap: int, shard_length: int) -> list[tuple[int, int, int]]:
    """Cut a sequence of ``length`` bases into ``(offset, stop, owned)`` shards.

    The piece is ``sequence[offset:stop]`` and the shard keeps hits starting in
    ``[offset, offset + owned)``. The owned regions **partition** ``[0, length)``, and each
    piece runs ``overlap`` bases past its own region — one less than the longest matrix, so
    a hit starting at the last owned position is still scored over every base it covers.

    Parameters
    ----------
    length : int
        How many bases the sequence holds.
    overlap : int
        One less than the longest matrix. Zero when there is nothing to scan with.
    shard_length : int
        The owned region to aim for. Raised to ``overlap + 1`` when it is smaller, since a
        shard owning less than that would be all overlap.

    Returns
    -------
    list of (int, int, int)
        ``(offset, stop, owned)`` per shard, in ascending order. Always at least one.

    Examples
    --------
    >>> plan_shards(100, overlap=14, shard_length=1000)      # one shard: nothing to cut
    [(0, 100, 100)]
    >>> plan_shards(100, overlap=14, shard_length=40)
    [(0, 54, 40), (40, 94, 40), (80, 100, 20)]
    >>> sum(owned for _offset, _stop, owned in plan_shards(100, 14, 40))
    100
    """
    step = max(shard_length, overlap + 1)
    if length <= step:
        return [(0, length, length)]
    plan = []
    for offset in range(0, length, step):
        owned = min(step, length - offset)
        plan.append((offset, min(offset + owned + overlap, length), owned))
    return plan


def parallel_batches(
    prepared: _Prepared, sequences: Iterable[tuple[str, str]], workers: int
) -> Iterator[pd.DataFrame]:
    """Scan across ``workers`` processes, yielding the batches a serial scan would yield.

    One **Hit table** batch per named sequence, in the order the sequences arrived, however
    many shards each was cut into — so the collector and the Parquet sink see exactly what
    they see from :func:`~genome.tf.motif.scan._batches`.

    Shards are submitted a bounded number ahead rather than all at once, so a FASTA is
    still streamed: what is in flight is ``workers`` times :data:`_IN_FLIGHT` pieces, never
    the file.

    Parameters
    ----------
    prepared : _Prepared
        The matrices, cutoffs and provenance of the scan. Pickled once per worker.
    sequences : iterable of (str, str)
        ``(name, sequence)`` pairs, drained once.
    workers : int
        How many processes. Two or more; one worker never reaches here.

    Yields
    ------
    pandas.DataFrame
        One batch per named sequence.

    Examples
    --------
    >>> import numpy as np
    >>> from genome.tf.motif import Motif, MotifSet
    >>> counts = np.zeros((4, 8))
    >>> for column, base in enumerate("GATTACAG"):
    ...     counts["ACGT".index(base), column] = 100.0
    >>> motifs = MotifSet([Motif("MA9999.1", "Gattacag", counts)])
    >>> serial = motifs.scan_sequences({"a": "TTGATTACAGTT", "b": "AAGATTACAGAA"})
    >>> shared = motifs.scan_sequences({"a": "TTGATTACAGTT", "b": "AAGATTACAGAA"}, workers=2)
    >>> list(serial["start"]) == list(shared["start"])
    True
    """
    overlap = max((len(motif) for motif in prepared.motifs), default=1) - 1
    context = multiprocessing.get_context(START_METHOD)
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=context,
        initializer=_install,
        initargs=(prepared,),
    ) as pool:
        gathered: list[pd.DataFrame] = []
        for shard, frame in _in_order(
            pool, _scan_shard, _shards(sequences, overlap), _IN_FLIGHT * workers
        ):
            gathered.append(frame)
            if len(gathered) == shard.of:
                yield _merge(prepared, gathered)
                gathered = []


def _shards(sequences: Iterable[tuple[str, str]], overlap: int) -> Iterator[_Shard]:
    """Cut each named sequence into the pieces :func:`plan_shards` says, lazily."""
    for name, sequence in sequences:
        plan = plan_shards(len(sequence), overlap, _SHARD_LENGTH)
        for offset, stop, owned in plan:
            yield _Shard(name, sequence[offset:stop], offset, owned, len(plan))


def _install(prepared: _Prepared) -> None:
    """Build this worker's engine, once, before it is given any shard to scan."""
    # A module global because a pool initializer has nowhere else to leave it: the tasks
    # that follow run in this process and reach it by name, and by nothing else.
    global _ENGINE
    _ENGINE = (prepared, engine_for(prepared))


def _scan_shard(shard: _Shard) -> pd.DataFrame:
    """Scan one shard with this worker's own engine, in the whole sequence's frame."""
    if _ENGINE is None:  # pragma: no cover - only reachable without the initializer
        raise RuntimeError(
            "the scan worker has no engine: _install() did not run. Reach this module "
            "through parallel_batches(), which is what installs one per worker."
        )
    prepared, engine = _ENGINE
    return _batch(
        prepared, engine, shard.name, shard.sequence, offset=shard.offset, owned=shard.owned
    )


def _in_order(
    pool: Executor,
    work: Callable[[_Shard], pd.DataFrame],
    shards: Iterable[_Shard],
    ahead: int,
) -> Iterator[tuple[_Shard, pd.DataFrame]]:
    """Run ``work`` over ``shards``, at most ``ahead`` outstanding, results in order.

    :meth:`concurrent.futures.Executor.map` submits the whole iterable at once, which for a
    genome would mean cutting every chromosome up before scanning any of it. This keeps a
    fixed window open instead, so the source stays a stream.
    """
    queue: deque[tuple[_Shard, Future[pd.DataFrame]]] = deque()
    remaining = iter(shards)
    for shard in islice(remaining, ahead):
        queue.append((shard, pool.submit(work, shard)))
    while queue:
        shard, future = queue.popleft()
        for following in islice(remaining, 1):
            queue.append((following, pool.submit(work, following)))
        yield shard, future.result()


def _merge(prepared: _Prepared, frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Put one sequence's shards back into the order a serial scan would have emitted them.

    Motif by motif in the set's own order, forward strand before reverse, ascending by
    position within each. The shards already arrive in ascending order, so the sort is
    stable and only regroups: it is what makes "identical to serial" a row-for-row claim
    rather than a set-equality one.
    """
    filled = [frame for frame in frames if not frame.empty]
    if not filled:
        return empty_hits()
    if len(filled) == 1:
        return filled[0]
    frame = pd.concat(filled, ignore_index=True).astype(dict(HIT_DTYPES))
    rank = {motif.motif_id: index for index, motif in enumerate(prepared.motifs)}
    order = np.lexsort(
        (
            frame["start"].to_numpy(),
            _ranked(frame["strand"], dict(zip(_STRANDS, range(len(_STRANDS)), strict=True))),
            _ranked(frame["motif_id"], rank),
        )
    )
    return frame.iloc[order].reset_index(drop=True)


def _ranked(column: Any, order: dict[str, int]) -> npt.NDArray[np.int64]:
    """Turn a categorical column into the sort key ``order`` gives its values.

    Through the categories rather than value by value, so the cost is the number of
    distinct values and not the number of rows — a merged chromosome is tens of millions of
    them. ``column`` is a ``pandas.Series``; indexing a frame is not typed as narrowly as
    that, and naming the narrower type here would only move the cast.
    """
    lookup = np.array([order[str(value)] for value in column.cat.categories], dtype=np.int64)
    return lookup[column.cat.codes.to_numpy()]
