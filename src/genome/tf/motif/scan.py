"""The scan engine adapter: MOODS driven for one **Motif set**, answering with a **Hit table**.

Every scan in this package goes through :func:`scan_stream`, whatever it was handed — a
single :class:`~genome.seq.DNA`, a mapping of named sequences, or a FASTA on disk. The
three entry points on :class:`~genome.tf.motif.motif.MotifSet` differ only in how they
turn their argument into the ``(name, sequence)`` pairs this module consumes, so the table
they hand back is one table with one schema.

**The engine has no strand.** MOODS scans a list of matrices along one sequence and knows
nothing about the second one. The adapter therefore doubles the list with the reverse
complement of every matrix and splits the results back by index: the first half is ``+``
and the second is ``-``. A position reported against a reverse-complement matrix is
already a **0-based half-open** start *in the forward frame* — the matrix moved, not the
sequence — so nothing is subtracted from it and both strands land in one coordinate
system. ``tests/test_scan.py`` asserts that against the engine rather than trusting it.

**The engine scores in nats and this package reports bits.** MOODS is handed
:meth:`~genome.tf.motif.motif.Motif.log_odds` scaled by ``ln 2`` and its scores are
divided by it again, so the log-odds arithmetic has exactly one implementation and the
``score`` column is in the same unit as **Information content**.

**Sequence is upper-cased before it is scanned**, so **Soft-masking** changes no answer and
there is no argument that would make it change one — ADR-0012. The engine's own alphabet
happens to fold case as well; the promise is this module's rather than the engine's, which
is why the call is made here.

**The Background is settled before the matrices are built, because it builds them.** That
is the whole reason :func:`scan_stream` hands its sequence source to
:func:`~genome.tf.motif.background.resolve_background` first and scans whatever comes back:
an automatic background has to see the input, and the input is an iterator that may be
drained only once. What comes back is the same records in the same order. The cutoffs those
matrices are called at are the engine's one slow step and a pure function of
``(matrices, background, p)``, so :func:`~genome.tf.motif.thresholds.cutoffs_for` answers
from disk on a repeat.

Batches, not one array: :func:`scan_stream` consumes an iterable of named sequences and
drains one frame per sequence, so the peak memory is set by the largest record rather than
by the file. That shape is deliberate, and both ends of it are now used —
:mod:`~genome.tf.motif.parquet`'s sink replaces the collector when an output path is given,
and :mod:`~genome.tf.motif.parallel`'s pool replaces the source above one worker — with
nothing between them changing, and the same batches in the same order either way.

Examples
--------
>>> import numpy as np
>>> from genome.tf.motif import Motif, MotifSet
>>> counts = np.zeros((4, 8))
>>> for column, base in enumerate("GATTACAG"):
...     counts["ACGT".index(base), column] = 100.0
>>> motifs = MotifSet([Motif("MA9999.1", "Gattacag", counts)])
>>> hits = motifs.scan("TTTTTGATTACAGTTTTT")
>>> hits[["motif_id", "sequence_name", "start", "end", "strand"]]
   motif_id sequence_name  start  end strand
0  MA9999.1      sequence      5   13      +
"""

from __future__ import annotations

import gzip
import math
from collections.abc import Hashable, Iterable, Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import IO, TYPE_CHECKING, Any, overload

import MOODS.scan
import MOODS.tools
import numpy as np
import pandas as pd

from genome.tf.motif.background import BackgroundArg, resolve_background
from genome.tf.motif.motif import DEFAULT_THRESHOLD, MIN_MOTIF_LENGTH, Motif
from genome.tf.motif.thresholds import cutoffs_for
from genome.tf.motif.workers import DEFAULT_WORKERS, resolve_workers

if TYPE_CHECKING:  # pragma: no cover - typing only, and the import runs the other way
    from genome.tf.motif.motif import MotifSet

#: The **Hit table**'s columns and their dtypes, in order. **This is the contract**, not an
#: optimisation: the compact forms cost 19 bytes a row where the dtypes pandas would infer
#: cost about 100, and a **Hit table** is the one thing every scan in this package returns,
#: so nothing downstream may have to ask which scan produced it.
HIT_DTYPES: Mapping[str, str] = MappingProxyType(
    {
        "motif_id": "category",
        "motif_name": "category",
        "sequence_name": "category",
        "start": "int32",
        "end": "int32",
        "strand": "category",
        "score": "float16",
    }
)

#: The **Hit table**'s columns, in order — the keys of :data:`HIT_DTYPES`.
HIT_COLUMNS: tuple[str, ...] = tuple(HIT_DTYPES)

#: What a **Hit table** records about the scan that made it, on ``frame.attrs``. Two tables
#: missing any of these six cannot be reconciled, which is why they travel with the rows.
HIT_PROVENANCE: tuple[str, ...] = (
    "background",
    "threshold",
    "release",
    "tax_group",
    "motifs_scanned",
    "motifs_skipped",
)

#: What :meth:`~genome.tf.motif.motif.Motif.log_odds` is given, so an unobserved base
#: scores very low rather than ``-inf``. The engine's own default, and the value the
#: measurement that chose the engine was taken at.
_PSEUDOCOUNT = 0.01

#: One nat in bits. MOODS scores in natural log and this package reports bits.
_BITS_PER_NAT = 1.0 / math.log(2.0)

#: The **Strand** each half of the doubled matrix list carries. Order matters: the forward
#: matrices come first and their reverse complements follow, which is what makes the split
#: an index comparison rather than anything stored.
_STRANDS: tuple[str, str] = ("+", "-")

#: How far ahead the engine's automaton looks. The engine's own default, and what
#: ``MOODS.scan.scan_dna`` passes when it builds a scanner per call — the same number, so
#: hoisting the construction out of the loop changes what is built and never what it says.
#: Unrelated to :data:`~genome.tf.motif.motif.MIN_MOTIF_LENGTH`, which happens to be 7 too.
_WINDOW_SIZE = 7


class FastaFormatError(ValueError):
    r"""A file handed to a scan is not FASTA, or holds a record with no name.

    A bad *file*, not a bad call. A **Hit table** is keyed by sequence name, so a record
    that has none could not be joined to anything and is refused rather than given one.

    Examples
    --------
    >>> from pathlib import Path
    >>> import tempfile
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     path = Path(directory) / "bad.fa"
    ...     _ = path.write_text("ACGTACGT\n")
    ...     try:
    ...         list(read_fasta(path))
    ...     except FastaFormatError as error:
    ...         print("'>'" in str(error))
    True
    """


@dataclass(frozen=True)
class _Prepared:
    """Everything a scan needs that does not depend on which sequence is being scanned.

    Built once per call and reused for every batch: the doubled matrix list, one cutoff per
    entry of it, and the provenance the **Hit table** will carry. Threshold computation is
    the engine's one slow step, which is why it is here and not inside the loop.
    """

    motifs: tuple[Motif, ...]
    skipped: tuple[str, ...]
    matrices: list[list[list[float]]]
    cutoffs: list[float]
    background: tuple[float, ...]
    threshold: float
    release: str | None
    tax_group: str | None


@overload
def scan_stream(
    motifs: MotifSet,
    sequences: Iterable[tuple[str, str]],
    *,
    threshold: float = ...,
    background: BackgroundArg = ...,
    output: None = ...,
    workers: int | None = ...,
) -> pd.DataFrame: ...


@overload
def scan_stream(
    motifs: MotifSet,
    sequences: Iterable[tuple[str, str]],
    *,
    threshold: float = ...,
    background: BackgroundArg = ...,
    output: str | Path,
    workers: int | None = ...,
) -> Path: ...


def scan_stream(
    motifs: MotifSet,
    sequences: Iterable[tuple[str, str]],
    *,
    threshold: float = DEFAULT_THRESHOLD,
    background: BackgroundArg = None,
    output: str | Path | None = None,
    workers: int | None = DEFAULT_WORKERS,
) -> pd.DataFrame | Path:
    """Scan named sequences with a **Motif set** and collect one **Hit table**.

    The one path every scan takes. ``sequences`` is drained once, one **Motif hit** batch
    per named sequence, and the batches are concatenated at the end — so the peak cost is
    the largest sequence rather than all of them.

    **Given an output path the batches go to Parquet instead** and the path comes back in
    place of the table, so the whole result is never held: 550 million rows is what hg38
    against a full vertebrate release comes to, and no guard here decides that is too many.

    Parameters
    ----------
    motifs : MotifSet
        The motifs to scan with. Those shorter than
        :data:`~genome.tf.motif.motif.MIN_MOTIF_LENGTH` are not scanned and are named in
        the result's ``motifs_skipped``.
    sequences : iterable of (str, str)
        ``(name, sequence)`` pairs. Names are what the ``sequence_name`` column carries and
        are not checked for uniqueness — a caller sharding one sequence wants them equal.
    threshold : float, default 1e-4
        The **Threshold**: one per-position p-value, converted per motif against
        ``background`` into the score that motif must clear. Must be in ``(0, 1)``.
    background : sequence of float or {"auto", "uniform", "derive"}, optional
        The **Background**: four frequencies over
        :data:`~genome.tf.motif.motif.BASES`, above zero and summing to 1, or one of the
        three modes. **Automatic when omitted** — derived from ``sequences`` when they
        hold at least :data:`~genome.tf.motif.background.BACKGROUND_FLOOR` unambiguous
        bases, uniform below that. Deciding reads a bounded prefix of ``sequences``, which
        is then scanned like the rest: the source is still drained exactly once.
    output : str or pathlib.Path, optional
        Where to stream the hits as Parquet. Omitted, the table comes back in memory.
    workers : int, default 1
        How many processes to shard the scan across. **One by default**, so importing this
        and calling it never starts a process unasked; ``None`` resolves the machine's
        allocation with :func:`~genome.tf.motif.parallel.resolve_workers`, which is what
        the command line passes. More than one produces the identical table.

    Returns
    -------
    pandas.DataFrame or pathlib.Path
        The **Hit table**: :data:`HIT_COLUMNS` with :data:`HIT_DTYPES`, and
        :data:`HIT_PROVENANCE` on ``frame.attrs``. Empty when nothing cleared its cutoff,
        which is a real answer and carries its provenance like any other. With ``output``,
        the path written — read it back with
        :func:`~genome.tf.motif.parquet.read_hits`, which restores both.

    Raises
    ------
    ValueError
        If ``threshold`` is not in ``(0, 1)``, or ``background`` is neither a mode nor four
        positive frequencies summing to 1.

    Examples
    --------
    >>> import numpy as np
    >>> from genome.tf.motif import Motif, MotifSet
    >>> counts = np.zeros((4, 8))
    >>> for column, base in enumerate("GATTACAG"):
    ...     counts["ACGT".index(base), column] = 100.0
    >>> motifs = MotifSet([Motif("MA9999.1", "Gattacag", counts)])
    >>> hits = scan_stream(motifs, [("chrTest", "TTTTTGATTACAGTTTTT")])
    >>> list(hits.itertuples(index=False, name=None))[0][:6]
    ('MA9999.1', 'Gattacag', 'chrTest', 5, 13, '+')
    >>> hits.attrs["background"]                  # under the floor, so uniform — and said so
    (0.25, 0.25, 0.25, 0.25)
    """
    # Imported here rather than at module scope: the sink is built on this module's schema,
    # so naming it up there would be a cycle — the same shape motif.py's import of this one
    # takes, and for the same reason.
    from genome.tf.motif.parquet import write_hits

    _check_threshold(threshold)
    count = resolve_workers(workers)
    frequencies, remaining = resolve_background(background, sequences)
    prepared = _prepare(motifs, threshold=threshold, background=frequencies)
    batches = _batches(prepared, remaining, count)
    if output is None:
        return _collect(batches, prepared)
    return write_hits(batches, output, _provenance(prepared))


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    r"""Yield ``(name, sequence)`` for each record of a FASTA, one record at a time.

    **The record name is the header up to its first whitespace**, which is what STAR and
    chromap write into an alignment produced from the same file — so a **Hit table** joins
    against that alignment without anyone renaming anything. Everything after the first
    whitespace is the record's description and is dropped.

    Plain or gzipped (``.gz``). One record is held at a time, never the file.

    Parameters
    ----------
    path : str or pathlib.Path
        The FASTA to read.

    Yields
    ------
    tuple of (str, str)
        The record name and its bases, in file order, exactly as written — case included,
        since it is the scan and not the reader that discards **Soft-masking**.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    FastaFormatError
        If bases appear before any header, or a header carries no name.

    Examples
    --------
    >>> from pathlib import Path
    >>> import tempfile
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     path = Path(directory) / "two.fa"
    ...     _ = path.write_text(">chrI  the first one\nACGT\nACGT\n>chrII\nTTTT\n")
    ...     list(read_fasta(path))
    [('chrI', 'ACGTACGT'), ('chrII', 'TTTT')]
    """
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"FASTA file not found: {source}")
    with _open_text(source) as handle:
        name: str | None = None
        chunks: list[str] = []
        for number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = _record_name(line, number, source), []
            elif line.strip():
                if name is None:
                    raise FastaFormatError(
                        f"{source} line {number}: bases before any header, so this is not "
                        f"FASTA. Every record opens with a '>' line naming it."
                    )
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks)


def _open_text(path: Path) -> IO[str]:
    """Open a FASTA for reading, transparently un-gzipping a ``.gz`` one."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _record_name(header: str, number: int, path: Path) -> str:
    """Read a header line's record name — everything before its first whitespace."""
    fields = header[1:].split(maxsplit=1)
    if not fields:
        raise FastaFormatError(
            f"{path} line {number}: a record header with no name. A hit table is keyed by "
            f"sequence name, so a nameless record could be joined to nothing — write the "
            f"name after the '>'."
        )
    return fields[0]


def _check_threshold(threshold: float) -> None:
    """Refuse a threshold that is not a p-value, before any sequence has been touched."""
    if not 0.0 < threshold < 1.0:
        raise ValueError(
            f"threshold must be a per-position p-value in (0, 1), got {threshold}. It is a "
            f"p-value and not a score — 1e-4 is the default, and a smaller number is "
            f"stricter."
        )


def _prepare(
    motifs: MotifSet,
    *,
    threshold: float,
    background: tuple[float, ...],
) -> _Prepared:
    """Convert a **Motif set** into matrices, cutoffs and the provenance of the scan.

    Motifs under :data:`~genome.tf.motif.motif.MIN_MOTIF_LENGTH` are set aside here rather
    than handed to the engine at a looser cutoff than was asked for; the engine would clamp
    a p-value it cannot reach to that matrix's best attainable score and over-call in
    silence.

    The **Background** arrives already settled and already quantised — this is downstream
    of the one decision, not a second copy of it.
    """
    frequencies = np.asarray(background, dtype=np.float64)
    scanned = tuple(motif for motif in motifs if len(motif) >= MIN_MOTIF_LENGTH)
    skipped = tuple(motif.motif_id for motif in motifs if len(motif) < MIN_MOTIF_LENGTH)
    # The engine works in nats; log_odds() is in bits, and this is the only place the two
    # meet on the way in, as the score column is the only place they meet coming back.
    forward = [
        (motif.log_odds(frequencies, _PSEUDOCOUNT) / _BITS_PER_NAT).tolist() for motif in scanned
    ]
    reverse = [[list(row) for row in MOODS.tools.reverse_complement(matrix)] for matrix in forward]
    matrices = forward + reverse
    return _Prepared(
        motifs=scanned,
        skipped=skipped,
        matrices=matrices,
        # Off the disk when this triple has been asked for before: the conversion is the
        # engine's one slow step and a pure function of it.
        cutoffs=cutoffs_for(matrices, background, threshold),
        background=background,
        threshold=float(threshold),
        # Read off the set rather than isinstance-tested: the release type imports this
        # module's neighbours, and filter() hands back a plain set that has neither.
        release=getattr(motifs, "release", None),
        tax_group=getattr(motifs, "tax_group", None),
    )


def engine_for(prepared: _Prepared) -> MOODS.scan.Scanner:
    """Build the scanner once for a whole scan — the setup ``scan_dna`` pays per call.

    Exactly what ``MOODS.scan.scan_dna`` constructs internally, hoisted out of the loop so
    that a thousand-record FASTA builds one automaton rather than a thousand, and so that a
    parallel worker can build its own once and keep it for every shard it is given.
    """
    scanner = MOODS.scan.Scanner(_WINDOW_SIZE)
    scanner.set_motifs(prepared.matrices, list(prepared.background), prepared.cutoffs)
    return scanner


def _batches(
    prepared: _Prepared, sequences: Iterable[tuple[str, str]], workers: int = 1
) -> Iterator[pd.DataFrame]:
    """Yield one **Hit table** batch per named sequence, in the order they arrive.

    The seam the Parquet sink and the parallel source both attach to: the source is
    whatever iterable was handed in, and the sink is whoever drains this. Above one worker
    the loop below is replaced wholesale by
    :func:`~genome.tf.motif.parallel.parallel_batches`, which yields the same batches in
    the same order — one per named sequence, however many shards it was cut into.
    """
    if workers > 1:
        # Imported here rather than at module scope: the parallel source is built on this
        # module's engine, so naming it up there would be a cycle.
        from genome.tf.motif.parallel import parallel_batches

        yield from parallel_batches(prepared, sequences, workers)
        return
    engine = engine_for(prepared)
    for name, sequence in sequences:
        yield _batch(prepared, engine, name, sequence)


def _batch(
    prepared: _Prepared,
    engine: MOODS.scan.Scanner,
    name: str,
    sequence: str,
    *,
    offset: int = 0,
    owned: int | None = None,
) -> pd.DataFrame:
    """Scan one named sequence and turn what the engine said into rows.

    Rows come out motif by motif in the set's own order, forward strand before reverse and
    ascending by position within each — the engine's own order, regrouped by the index
    split that recovers **Strand**.

    ``offset`` and ``owned`` are what make a *shard* of a sequence report in the whole
    sequence's frame: every position is reported as ``offset + position``, and a hit
    starting at or past ``owned`` belongs to the next shard and is dropped here — so a hit
    spanning a shard boundary is found once by the shard that owns its start and once only.
    """
    # Upper-cased here and nowhere else: soft-masking changes no answer (ADR-0012).
    results = engine.scan(sequence.upper())
    count = len(prepared.motifs)
    rows: list[tuple[str, str, str, int, int, str, float]] = []
    for index, motif in enumerate(prepared.motifs):
        length = len(motif)
        for half, strand in enumerate(_STRANDS):
            for hit in results[index + half * count]:
                start = int(hit.pos)
                if owned is not None and start >= owned:
                    continue
                rows.append(
                    (
                        motif.motif_id,
                        motif.motif_name,
                        name,
                        offset + start,
                        offset + start + length,
                        strand,
                        float(hit.score) * _BITS_PER_NAT,
                    )
                )
    if not rows:
        return empty_hits()
    return pd.DataFrame(rows, columns=pd.Index(list(HIT_COLUMNS))).astype(dict(HIT_DTYPES))


def _collect(batches: Iterable[pd.DataFrame], prepared: _Prepared) -> pd.DataFrame:
    """Concatenate the batches into one **Hit table** and attach the scan's provenance.

    Empty batches are dropped rather than concatenated: a sequence with no hit contributes
    no row, and pandas would otherwise have to decide the result's dtypes from a frame that
    holds none.
    """
    filled = [batch for batch in batches if not batch.empty]
    frame = pd.concat(filled, ignore_index=True) if filled else empty_hits()
    # Re-applied after the concatenation: two batches whose categories differ collapse to
    # object, and the dtypes are the contract rather than whatever survived the join.
    frame = frame.astype(dict(HIT_DTYPES))
    frame.attrs = _provenance(prepared)
    return frame


def _provenance(prepared: _Prepared) -> dict[Hashable, Any]:
    """Describe the scan that made a **Hit table**, for ``frame.attrs``."""
    return {
        "background": prepared.background,
        "threshold": prepared.threshold,
        "release": prepared.release,
        "tax_group": prepared.tax_group,
        "motifs_scanned": tuple(motif.motif_id for motif in prepared.motifs),
        "motifs_skipped": prepared.skipped,
    }


def empty_hits() -> pd.DataFrame:
    """Return a **Hit table** of no rows, with the columns and dtypes of a full one.

    What a scan that found nothing answers with, so a caller never branches on emptiness to
    learn the schema.

    Returns
    -------
    pandas.DataFrame
        Zero rows, :data:`HIT_COLUMNS` and :data:`HIT_DTYPES`. Carries no provenance —
        :func:`scan_stream` attaches that to whatever it returns.

    Examples
    --------
    >>> empty_hits().shape
    (0, 7)
    >>> empty_hits().dtypes["start"]
    dtype('int32')
    """
    return pd.DataFrame({name: pd.Series(dtype=dtype) for name, dtype in HIT_DTYPES.items()})
