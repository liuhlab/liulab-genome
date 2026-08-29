"""Where the **Background** a scan scores against comes from, and what shape it is left in.

The **Background** decides the answer more than any other scan parameter — moving from
uniform to one real chromosome's composition changed the hit count by 2.5% and turned over
26% of the hits — so it is neither assumed nor left to the caller to remember. It is
**automatic**: derived from the input when the input holds at least
:data:`BACKGROUND_FLOOR` unambiguous bases, uniform below that, since a composition
estimated from fewer would distort the very cutoffs it sets. At that floor the standard
error on each base frequency is about 0.004, which moves a per-position log-odds term by
under 0.02 nats.

**Derivation reads a bounded prefix of the input, not all of it.** A concrete background
has to exist before the matrices are built, and the matrices have to exist before the first
sequence is scanned — so the alternative would be a second pass over the source, which a
FASTA can afford only by being re-read and a generator cannot afford at all.
:func:`resolve_background` instead pulls whole records off the source until the floor is
reached, decides, and hands back the records it pulled chained in front of whatever is
left: **the source is consumed exactly once**, and what is held while the decision is made
is the same one record the scan loop holds anyway. The cost is that the estimate is a head
sample of the input rather than the whole of it — which is exactly the accuracy the floor
was chosen for, since 10,000 bases is by construction enough.

**Whatever the background turns out to be, it is quantised to a 0.001 grid** — three
decimals, well inside the 0.004 standard error above — and the quantised value is what
builds the matrices, what sets the cutoffs, what keys the threshold cache, and what the
**Hit table** records. One number, used everywhere and reported once, so a repeat scan
cannot answer differently depending on whose entry it found in the cache. Quantising is
the identity on any background written with three decimals or fewer, which is every
background anyone writes by hand.

Examples
--------
>>> resolve_background("uniform", [("chrTest", "AAAACCCC")])[0]
(0.25, 0.25, 0.25, 0.25)
>>> background, sequences = resolve_background("derive", [("chrTest", "AACC")])
>>> background
(0.375, 0.375, 0.125, 0.125)
>>> list(sequences)                       # the source is intact, not consumed
[('chrTest', 'AACC')]
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from itertools import chain
from typing import Literal

import numpy as np
import numpy.typing as npt

from genome.tf.motif.motif import BASES, _check_background

#: How many unambiguous bases the input must hold before a **Background** is derived from
#: it rather than assumed uniform. At this floor the standard error on each base frequency
#: is about 0.004, which moves a per-position log-odds term by under 0.02 nats; below it
#: the estimate would distort the cutoffs it exists to set.
BACKGROUND_FLOOR = 10_000

#: The three things a caller can ask for by name instead of handing over four frequencies.
#: ``"auto"`` is the default and is what ``None`` means: derive above
#: :data:`BACKGROUND_FLOOR`, uniform below it. ``"uniform"`` pins it. ``"derive"`` derives
#: whatever the input holds, floor or no floor.
BACKGROUND_MODES: tuple[str, str, str] = ("auto", "uniform", "derive")

#: The **Background** that assumes nothing: each base a quarter.
UNIFORM_BACKGROUND: tuple[float, float, float, float] = (0.25, 0.25, 0.25, 0.25)

#: What a caller may pass for a background: four frequencies, one of
#: :data:`BACKGROUND_MODES`, or ``None`` for ``"auto"``.
BackgroundMode = Literal["auto", "uniform", "derive"]
BackgroundArg = Sequence[float] | npt.NDArray[np.float64] | BackgroundMode | None

#: The grid a background is rounded onto — three decimals. Coarse enough that two peak sets
#: from one genome land on one entry of the threshold cache, and eight times finer than the
#: 0.004 standard error :data:`BACKGROUND_FLOOR` is chosen for.
_GRID = 1000

#: Added to each base's count before a background is derived, so a sequence holding no ``C``
#: at all still yields four positive frequencies rather than an infinite log-odds. At the
#: floor it moves a frequency by under 1e-4 — a tenth of the grid, so it is invisible.
_DERIVATION_PSEUDOCOUNT = 1.0

#: The most characters :func:`resolve_background` will pull off the source while deciding.
#: The floor counts unambiguous bases only, so without this a file of nothing but ``N``
#: would be buffered whole — and no scan may hold a genomic file in memory.
_SAMPLE_CAP = 100 * BACKGROUND_FLOOR


def resolve_background(
    background: BackgroundArg,
    sequences: Iterable[tuple[str, str]],
) -> tuple[tuple[float, ...], Iterator[tuple[str, str]]]:
    """Settle the **Background** for a scan, and hand back the sequences still to scan.

    Returns both because deciding may have had to look at the input: the records read
    while deciding come back chained in front of the rest, so the caller drains one
    iterator and the source is consumed exactly once. Nothing is read at all when the
    background was given outright or pinned to uniform.

    Parameters
    ----------
    background : sequence of float or {"auto", "uniform", "derive"} or None
        Four frequencies over :data:`~genome.tf.motif.motif.BASES`, above zero and summing
        to 1; or one of :data:`BACKGROUND_MODES`; or ``None``, which is ``"auto"``.
    sequences : iterable of (str, str)
        The ``(name, sequence)`` pairs the scan is about to consume.

    Returns
    -------
    tuple of float
        The background to scan with, quantised to the 0.001 grid — four frequencies in
        :data:`~genome.tf.motif.motif.BASES` order.
    iterator of (str, str)
        Every pair ``sequences`` would have yielded, in the same order, including any
        pulled off while deciding.

    Raises
    ------
    ValueError
        If ``background`` names a mode that does not exist, or is not four positive
        frequencies summing to 1.

    Examples
    --------
    >>> resolve_background(None, [("chrTest", "ACGTACGT")])[0]   # under the floor
    (0.25, 0.25, 0.25, 0.25)
    >>> resolve_background([0.3, 0.2, 0.2, 0.3], [])[0]
    (0.3, 0.2, 0.2, 0.3)
    >>> background, sequences = resolve_background("auto", [("chrX", "ACGTTTTTTT" * 4000)])
    >>> background                                              # over the floor, so derived
    (0.1, 0.1, 0.1, 0.7)
    >>> [(name, len(bases)) for name, bases in sequences]        # nothing was eaten
    [('chrX', 40000)]
    """
    if background is None:
        mode = "auto"
    elif isinstance(background, str):
        mode = _check_mode(background)
    else:
        return quantise_background(_check_background(background)), iter(sequences)
    if mode == "uniform":
        return UNIFORM_BACKGROUND, iter(sequences)
    counts, unambiguous, buffered, source = _sample(sequences)
    rest = chain(buffered, source)
    if mode == "auto" and unambiguous < BACKGROUND_FLOOR:
        return UNIFORM_BACKGROUND, rest
    return derive_background(counts), rest


def derive_background(counts: Sequence[float]) -> tuple[float, ...]:
    """Turn observed base counts into a quantised **Background**.

    Counts are in :data:`~genome.tf.motif.motif.BASES` order.
    :data:`_DERIVATION_PSEUDOCOUNT` is added to each, so a base that was never seen gets a
    small positive frequency instead of one that makes every log-odds against it infinite.

    Parameters
    ----------
    counts : sequence of float
        Four observed counts, in ``A``, ``C``, ``G``, ``T`` order.

    Returns
    -------
    tuple of float
        Four frequencies on the 0.001 grid, summing to 1.

    Examples
    --------
    >>> derive_background([250, 250, 250, 250])
    (0.25, 0.25, 0.25, 0.25)
    >>> derive_background([600, 200, 100, 100])
    (0.598, 0.2, 0.101, 0.101)
    >>> derive_background([90, 0, 0, 10])                # an unseen base still scores
    (0.874, 0.01, 0.01, 0.106)
    """
    weighted = np.asarray(counts, dtype=np.float64) + _DERIVATION_PSEUDOCOUNT
    return quantise_background(weighted / weighted.sum())


def quantise_background(
    frequencies: Sequence[float] | npt.NDArray[np.float64],
) -> tuple[float, ...]:
    """Round a background onto the 0.001 grid, keeping it four positive frequencies.

    Rounding alone would leave the four summing to something near 1 rather than 1, so the
    residual goes to the largest of them — which is what keeps the returned value a
    background and not merely four numbers near one. No frequency is ever rounded to zero.

    Parameters
    ----------
    frequencies : sequence of float
        Four frequencies summing to 1.

    Returns
    -------
    tuple of float
        Four frequencies, each a multiple of 0.001, summing to 1.

    Examples
    --------
    >>> quantise_background([0.3, 0.2, 0.2, 0.3])       # already on the grid
    (0.3, 0.2, 0.2, 0.3)
    >>> quantise_background([1 / 3, 1 / 3, 1 / 6, 1 / 6])
    (0.333, 0.333, 0.167, 0.167)
    """
    scaled = np.maximum(np.rint(np.asarray(frequencies, dtype=np.float64) * _GRID), 1.0)
    grid_counts = scaled.astype(np.int64)
    grid_counts[int(np.argmax(grid_counts))] += _GRID - int(grid_counts.sum())
    return tuple(float(count) / _GRID for count in grid_counts)


def _check_mode(mode: str) -> str:
    """Return the mode asked for, refusing a name that is not one of the three."""
    if mode not in BACKGROUND_MODES:
        raise ValueError(
            f"background must be four frequencies or one of {', '.join(BACKGROUND_MODES)}, "
            f"got {mode!r}. Omit it, or pass 'auto', to derive above {BACKGROUND_FLOOR} "
            f"unambiguous bases and stay uniform below."
        )
    return mode


def _sample(
    sequences: Iterable[tuple[str, str]],
) -> tuple[list[int], int, list[tuple[str, str]], Iterator[tuple[str, str]]]:
    """Pull records off the source until the floor is reached, counting bases as they go.

    Stops at :data:`_SAMPLE_CAP` characters whatever the counts have reached, so nothing
    unbounded is ever buffered. Returns the counts, how many of them were unambiguous, the
    records pulled, and the source they were pulled from.
    """
    counts = [0, 0, 0, 0]
    unambiguous = 0
    examined = 0
    buffered: list[tuple[str, str]] = []
    source = iter(sequences)
    for record in source:
        buffered.append(record)
        unambiguous += _count_bases(record[1], counts)
        examined += len(record[1])
        if unambiguous >= BACKGROUND_FLOOR or examined >= _SAMPLE_CAP:
            break
    return counts, unambiguous, buffered, source


def _count_bases(sequence: str, counts: list[int]) -> int:
    """Add one sequence's unambiguous bases to ``counts``, returning how many there were."""
    found = 0
    for index, base in enumerate(BASES):
        seen = sequence.count(base) + sequence.count(base.lower())
        counts[index] += seen
        found += seen
    return found
