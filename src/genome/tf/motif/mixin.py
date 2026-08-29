"""Mixin giving a :class:`~genome.genome.Genome` one scan over **Region**s.

**The one place region-local positions become Chromosome coordinates.** Scanning a
fetched **Region** answers in the frame of the bases that were handed to the engine,
which starts at zero for every region; turning that into the assembly's frame means
adding the region's start, and for a ``-`` **Strand** region flipping the interval and
the hit strand as well. That arithmetic is written once, here, rather than in every
notebook that ever scanned a peak set — it is where the off-by-ones live.

The dependency runs **Genome to motif and never back**: this module names
:class:`~genome.genome.Genome` under type checking alone, so a **Motif set** stays usable
with no genome open and a motif goes on belonging to no **Assembly**. The scan itself is
the one every other entry point takes, called once with every region's bases
(:meth:`~genome.tf.motif.motif.MotifSet.scan_sequences`) — so the **Hit table** that comes
back is that table, with its coordinates lifted and its **Assembly** named on top of the
provenance already there.

Examples
--------
>>> from genome import Genome, Region
>>> from genome.tf.motif import JasparDatabase
>>> sacCer3 = Genome("sacCer3")                                   # doctest: +SKIP
>>> ctcf = JasparDatabase().filter(lambda motif: motif.motif_name == "CTCF")  # doctest: +SKIP
>>> hits = sacCer3.scan_regions(ctcf, Region("chrI", 5000, 6000, "-"))       # doctest: +SKIP
>>> hits[["sequence_name", "start", "end", "strand"]].head(1)     # doctest: +SKIP
  sequence_name  start   end strand
0          chrI   5140  5155      +
>>> hits.attrs["assembly"]                                        # doctest: +SKIP
'sacCer3'
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import TYPE_CHECKING, cast

import numpy as np
import numpy.typing as npt
import pandas as pd

from genome.region import Region
from genome.tf.motif.motif import DEFAULT_THRESHOLD
from genome.tf.motif.scan import HIT_DTYPES, HIT_PROVENANCE

if TYPE_CHECKING:  # pragma: no cover - typing only, and the dependency runs one way
    from genome.genome import Genome
    from genome.tf.motif.motif import MotifSet

#: What a **Hit table** scanned from **Region**s carries on ``frame.attrs``: everything
#: :data:`~genome.tf.motif.scan.HIT_PROVENANCE` records about the scan, and the
#: **Assembly** its coordinates are in. Extended rather than replaced — one table, one
#: mechanism — and the assembly is recorded because it is never inferable from a row: a
#: chromosome name and an interval mean nothing until something says which reference they
#: are in (ADR-0003).
REGION_HIT_PROVENANCE: tuple[str, ...] = (*HIT_PROVENANCE, "assembly")


class MotifScanMixin:
    """Scan **Region**s of this **Genome** and get hits in the assembly's own frame.

    One method, because there is one thing a genome adds to a scan: the coordinates. The
    raw forms stay exactly as they were — a :class:`~genome.tf.motif.motif.MotifSet`
    handed a string, a mapping of sequences or a FASTA still answers in the frame of what
    it was given and names no assembly.
    """

    def scan_regions(
        self,
        motifs: MotifSet,
        regions: Region | Iterable[Region],
        *,
        threshold: float = DEFAULT_THRESHOLD,
        background: Sequence[float] | npt.NDArray[np.float64] | None = None,
    ) -> pd.DataFrame:
        """Scan **Region**s of this genome and return the **Hit table**, in its coordinates.

        Each region's bases are fetched exactly as
        :meth:`~genome.genome.Genome.fetch_sequence` returns them — reverse-complemented
        for a ``-`` region — scanned in one call, and the hits lifted into the assembly's
        frame. ``sequence_name`` carries the **Chromosome** name as *this assembly* spells
        it, since a region whose chromosome the assembly does not carry raises rather than
        being reconciled.

        **The arithmetic, so it can be checked without being run.** Write the region as
        ``[S, E)`` on a chromosome, ``L = E - S``, and a hit as the **0-based half-open**
        ``[s, e)`` the scan found in the fetched bases:

        - ``+``, and ``.`` with it: the bases run along the chromosome, so local ``i`` is
          chromosome ``S + i``. The hit is ``[S + s, S + e)`` and its **Strand** is
          unchanged. An unknown strand is *not* promoted to ``+`` — it is that the fetch
          returns forward bases for it, so there is nothing to flip.
        - ``-``: the bases are the reverse complement of ``[S, E)``, so local ``i`` is
          chromosome ``E - 1 - i``. A hit covering local ``s .. e - 1`` therefore covers
          chromosome ``E - e .. E - s - 1``, which as a half-open interval is
          ``[E - e, E - s)``: **the two ends swap**, and the ``- 1``s cancel exactly
          because the interval is half-open — in 1-based-inclusive coordinates they would
          not. Its strand flips too: what matched the forward strand of the fetched bases
          matched the reverse strand of the chromosome.

        Two checks on the ``-`` case: a hit spanning the whole region, ``s = 0`` and
        ``e = L``, comes back as ``[E - L, E) = [S, E)``; and ``(E - s) - (E - e) = e - s``,
        so the length is preserved.

        Parameters
        ----------
        motifs : MotifSet
            The motifs to scan with — a **Release**, a filtered set, or de novo matrices.
            Those shorter than :data:`~genome.tf.motif.motif.MIN_MOTIF_LENGTH` are not
            scanned and are named in the result's ``motifs_skipped``.
        regions : Region or iterable of Region
            One region or many, on any chromosomes, in any order and overlapping freely:
            two regions on one chromosome are ordinary here, and each is scanned in its
            own right. A locus *string* is not accepted, because it carries no **Strand**
            and the strand is the whole question — write
            ``Region.from_string("chrI:0-600", strand="-")``.
        threshold : float, default 1e-4
            The **Threshold**: one per-position p-value, converted per motif against
            ``background`` into the score that motif must clear. In ``(0, 1)``.
        background : sequence of float, optional
            The **Background**: four frequencies over
            :data:`~genome.tf.motif.motif.BASES`, above zero and summing to 1. Uniform
            when omitted.

        Returns
        -------
        pandas.DataFrame
            The **Hit table** — :data:`~genome.tf.motif.scan.HIT_COLUMNS` with
            :data:`~genome.tf.motif.scan.HIT_DTYPES`, one row per **Motif hit** — with
            ``sequence_name``, ``start``, ``end`` and ``strand`` in this assembly's frame.
            :data:`REGION_HIT_PROVENANCE` is on ``frame.attrs``: the scan's own provenance
            and the **Assembly** these coordinates belong to. Empty when nothing cleared
            its cutoff, or when no region was given, and carrying its provenance either way.

        Raises
        ------
        TypeError
            If ``regions`` holds anything that is not a :class:`~genome.region.Region`.
        ValueError
            If a region names a chromosome this assembly does not carry or falls outside
            it; or if ``threshold`` is not in ``(0, 1)`` or ``background`` is not four
            positive frequencies summing to 1.

        See Also
        --------
        genome.tf.motif.motif.MotifSet.scan_sequences : the same scan, region-local and
            naming no assembly.

        Examples
        --------
        >>> from genome import Genome, Region
        >>> from genome.tf.motif import JasparDatabase
        >>> sacCer3 = Genome("sacCer3")                              # doctest: +SKIP
        >>> peaks = [Region("chrI", 0, 5000, "+"), Region("chrII", 100, 900, "-")]
        >>> hits = sacCer3.scan_regions(JasparDatabase(), peaks)     # doctest: +SKIP
        >>> sorted(set(hits["sequence_name"]))                       # doctest: +SKIP
        ['chrI', 'chrII']
        """
        genome = cast("Genome", self)
        wanted = _as_regions(regions)
        # Keyed by position and not by chromosome: two regions on one chromosome are the
        # everyday case, and a mapping keyed by name could hold only one of them. The key
        # is replaced by the chromosome on the way out, so nothing outside this call ever
        # sees it — and one call means the engine's slow threshold step is paid once.
        sequences = {
            str(index): genome.fetch_sequence(region) for index, region in enumerate(wanted)
        }
        hits = motifs.scan_sequences(sequences, threshold=threshold, background=background)
        return _lifted(hits, wanted, genome.assembly)


def _as_regions(regions: Region | Iterable[Region]) -> tuple[Region, ...]:
    """Read the ``regions`` argument as the **Region**s it names — one, or an iterable."""
    wanted = (regions,) if isinstance(regions, Region) else tuple(regions)
    offender = next((item for item in wanted if not isinstance(item, Region)), None)
    if offender is not None:
        raise TypeError(
            f"scan_regions takes a Region or an iterable of them, and got {offender!r}. A "
            f"locus string is not one: it carries no strand, and the strand is what decides "
            f"whether the hits are flipped — build one with "
            f"Region.from_string('chrI:0-600', strand='-')."
        )
    return wanted


def _lifted(hits: pd.DataFrame, regions: tuple[Region, ...], assembly: str) -> pd.DataFrame:
    """Lift a region-local **Hit table** into ``assembly``'s frame and name the assembly.

    Vectorised over the whole table: every hit knows which region it came from by the key
    the scan was given, and the four columns that depend on the region are recomputed at
    once. The dtypes are re-applied because they are the contract rather than whatever
    survived the arithmetic, and the provenance is carried over rather than propagated,
    since a frame method is free to drop ``attrs``.
    """
    at = hits["sequence_name"].astype("int64").to_numpy()
    count = len(regions)
    starts = np.fromiter((region.start for region in regions), dtype="int64", count=count)
    ends = np.fromiter((region.end for region in regions), dtype="int64", count=count)
    reverse = np.fromiter((region.strand == "-" for region in regions), dtype=bool, count=count)
    chroms = np.array([region.chrom for region in regions], dtype=object)

    region_start, region_end, flip = starts[at], ends[at], reverse[at]
    local_start = hits["start"].to_numpy("int64")
    local_end = hits["end"].to_numpy("int64")
    strand = hits["strand"].astype("object").to_numpy()

    frame = hits.assign(
        sequence_name=chroms[at],
        # The ends swap for a reverse region, and only there; see scan_regions' docstring
        # for why the half-open convention is what makes this symmetric.
        start=np.where(flip, region_end - local_end, region_start + local_start),
        end=np.where(flip, region_end - local_start, region_start + local_end),
        strand=np.where(flip, np.where(strand == "+", "-", "+"), strand),
    ).astype(dict(HIT_DTYPES))
    frame.attrs = {**hits.attrs, "assembly": assembly}
    return frame
