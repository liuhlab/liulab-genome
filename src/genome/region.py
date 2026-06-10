"""Genomic coordinate primitives — the shared region/interval abstraction.

This module is the home for region, interval, and BED-style coordinate logic
that the rest of the package builds on. The core type is :class:`Region`: a
single genomic interval stored in the package's canonical internal convention —
**0-based, half-open** ``[start, end)`` — with an explicit strand.

All coordinates everywhere in this package are 0-based half-open (the BED
convention).

Examples
--------
>>> from genome.region import Region, parse_region
>>> Region("chr1", 0, 10)
Region(chrom='chr1', start=0, end=10, strand='.')
>>> len(Region("chr1", 0, 10))
10
>>> str(Region("chr1", 0, 10))                # 0-based, as stored
'chr1:0-10'
>>> parse_region("chr1:0-10")
('chr1', 0, 10)
>>> parse_region("chrM")                       # bare chromosome, no coordinates
('chrM', None, None)
"""

from __future__ import annotations

import re
from dataclasses import dataclass

#: A ``chrom:start-end`` locus, with optional thousands separators in the
#: coordinates (e.g. ``chr1:1,000-2,000``). Coordinates are **0-based**.
_RANGE_RE = re.compile(r"^(?P<chrom>.+):(?P<start>[\d,]+)-(?P<end>[\d,]+)$")

#: Allowed strand markers: forward, reverse, or unknown. Never defaulted to ``+``.
_STRANDS = frozenset({"+", "-", "."})


def parse_region(text: str) -> tuple[str, int | None, int | None]:
    """Parse a ``chrom:start-end`` locus string into its parts.

    Coordinates are interpreted as **0-based, half-open** ``[start, end)`` (the
    BED convention). Thousands separators in the numbers are tolerated. A bare
    chromosome name (no ``:start-end``) returns ``(chrom, None, None)`` so the
    caller can resolve it to the whole sequence against known chromosome sizes.

    Parameters
    ----------
    text : str
        ``"chrom:start-end"`` (0-based) or a bare ``"chrom"``.

    Returns
    -------
    tuple[str, int | None, int | None]
        ``(chrom, start, end)`` — ``start``/``end`` are ``None`` for a bare
        chromosome.

    Raises
    ------
    ValueError
        If ``text`` looks like a coordinate locus (contains ``:``) but is not a
        well-formed ``chrom:start-end``.

    Examples
    --------
    >>> parse_region("chr1:0-10")
    ('chr1', 0, 10)
    >>> parse_region("chr1:1,000-2,000")
    ('chr1', 1000, 2000)
    >>> parse_region("chrM")
    ('chrM', None, None)
    """
    text = text.strip()
    match = _RANGE_RE.match(text)
    if match:
        start = int(match["start"].replace(",", ""))
        end = int(match["end"].replace(",", ""))
        return match["chrom"], start, end
    if ":" in text:
        raise ValueError(
            f"malformed region {text!r}: expected 'chrom:start-end' (0-based, "
            f"half-open) or a bare chromosome name."
        )
    return text, None, None


@dataclass(frozen=True)
class Region:
    """A single genomic interval in 0-based, half-open coordinates.

    Parameters
    ----------
    chrom : str
        Chromosome / sequence name.
    start : int
        0-based start, inclusive. Must be ``>= 0``.
    end : int
        0-based end, exclusive (half-open). Must be ``>= start``.
    strand : str, default ``"."``
        ``"+"``, ``"-"``, or ``"."`` (unknown). Never silently defaulted to a
        real strand.

    Raises
    ------
    ValueError
        If ``start < 0``, ``end < start``, or ``strand`` is not one of
        ``"+"``, ``"-"``, ``"."``.

    Examples
    --------
    >>> r = Region("chr1", 0, 10)
    >>> len(r), r.length
    (10, 10)
    >>> Region.from_string("chr2:100-200", strand="-")
    Region(chrom='chr2', start=100, end=200, strand='-')
    """

    chrom: str
    start: int
    end: int
    strand: str = "."

    def __post_init__(self) -> None:
        """Validate the coordinate and strand invariants."""
        if self.start < 0:
            raise ValueError(f"start must be >= 0 (0-based), got {self.start}.")
        if self.end < self.start:
            raise ValueError(f"end ({self.end}) must be >= start ({self.start}).")
        if self.strand not in _STRANDS:
            raise ValueError(f"strand must be one of '+', '-', '.', got {self.strand!r}.")

    def __len__(self) -> int:
        """Return the number of bases spanned, ``end - start``."""
        return self.end - self.start

    @property
    def length(self) -> int:
        """Number of bases spanned, ``end - start`` (alias of ``len(self)``)."""
        return self.end - self.start

    def __str__(self) -> str:
        """Return the 0-based locus string ``chrom:start-end`` (as stored)."""
        return f"{self.chrom}:{self.start}-{self.end}"

    @classmethod
    def from_string(cls, text: str, *, strand: str = ".") -> Region:
        """Build a :class:`Region` from a 0-based ``chrom:start-end`` string.

        Parameters
        ----------
        text : str
            A ``chrom:start-end`` locus in 0-based half-open coordinates.
        strand : str, default ``"."``
            Strand to attach (the string itself carries no strand).

        Raises
        ------
        ValueError
            If ``text`` is malformed or carries no coordinates (a bare
            chromosome cannot be sized without a chromosome-length table).

        Examples
        --------
        >>> Region.from_string("chr1:0-10")
        Region(chrom='chr1', start=0, end=10, strand='.')
        """
        chrom, start, end = parse_region(text)
        if start is None or end is None:
            raise ValueError(
                f"region {text!r} has no coordinates; Region needs an explicit "
                f"start-end (resolve a whole chromosome via Genome instead)."
            )
        return cls(chrom, start, end, strand)
