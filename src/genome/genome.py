"""The :class:`Genome` class — the package's main entry point.

A :class:`Genome` ties an assembly name to its on-disk reference files and the
operations over them. Constructing one downloads and prepares everything needed
(FASTA, ``.fai`` index, ``.2bit``, ``chrom.sizes``) behind the scenes; from
there you query it directly — e.g. fetch the sequence of a locus as a
:class:`~genome.seq.DNA`.

Coordinates follow the package's canonical internal convention everywhere:
**0-based, half-open** ``[start, end)`` (the BED convention). See
:mod:`genome.region`.

This is a pilot surface: sequence retrieval today, more genome operations later.

Examples
--------
>>> from genome import Genome
>>> sacCer3 = Genome("sacCer3")                  # download + prepare on first use  # doctest: +SKIP
>>> sacCer3.fetch_sequence("chrIV:0-10")         # 0-based, half-open  # doctest: +SKIP
DNA('ACACCACACC')
>>> sacCer3["chrIV:0-10"]                         # indexing is sugar for fetch_sequence  # doctest: +SKIP
DNA('ACACCACACC')
"""

from __future__ import annotations

from pathlib import Path
from types import TracebackType
from typing import Self

import pandas as pd

from genome.io.download import UCSCGenomeDownloader
from genome.io.fasta import GenomeFiles, read_chrom_sizes
from genome.io.twobit import TwoBit
from genome.region import Region, parse_region
from genome.seq import DNA


class Genome:
    """A reference genome and the operations over it.

    Constructing a ``Genome`` ensures the assembly's reference files exist
    locally: the FASTA is downloaded from the UCSC golden path (if not already
    cached) and its ``.fai`` index, ``.2bit`` encoding, and ``chrom.sizes`` are
    prepared. All of this is cached under
    ``<LIULAB_DATA>/genome/<assembly>/`` (see
    :func:`~genome.io.download.assembly_data_dir`), so repeat constructions are
    cheap and offline.

    Sequence is read from the ``.2bit`` file via ``py2bit``; coordinates are
    **0-based, half-open** throughout.

    Parameters
    ----------
    assembly : str
        UCSC assembly name, e.g. ``"sacCer3"``, ``"hg38"``, ``"mm39"``. Validated
        against UCSC before any download, so a typo fails fast.
    cache_dir : str or pathlib.Path, optional
        Override the storage directory for this assembly's files. Defaults to the
        shared per-assembly reference directory.
    progressbar : bool, default True
        Show a download progress bar on first fetch (requires ``tqdm``).

    Attributes
    ----------
    assembly : str
        The assembly name.
    files : genome.io.fasta.GenomeFiles
        Paths to the prepared FASTA and its derived index/companion files.

    Raises
    ------
    ValueError
        If ``assembly`` is unknown to UCSC.
    genome.external.ToolNotFoundError
        If a required native tool (``samtools``, ``faToTwoBit``, ``twoBitInfo``)
        is not on ``PATH``.

    Examples
    --------
    >>> sacCer3 = Genome("sacCer3")               # doctest: +SKIP
    >>> sacCer3.fetch_sequence("chrIV:0-10")      # doctest: +SKIP
    DNA('ACACCACACC')
    """

    def __init__(
        self,
        assembly: str,
        *,
        cache_dir: str | Path | None = None,
        progressbar: bool = True,
    ) -> None:
        self.assembly = assembly
        self._downloader = UCSCGenomeDownloader(assembly, cache_dir)
        self.files: GenomeFiles = self._downloader.fetch_genome(progressbar=progressbar)
        self._chrom_sizes: pd.Series = read_chrom_sizes(self.files.chrom_sizes)
        self._twobit = TwoBit(self.files.twobit)

    def __repr__(self) -> str:
        """Return e.g. ``Genome('sacCer3', 17 sequences)``."""
        return f"{type(self).__name__}({self.assembly!r}, {len(self._chrom_sizes)} sequences)"

    def close(self) -> None:
        """Release the open 2bit file handle (idempotent)."""
        twobit = getattr(self, "_twobit", None)
        if twobit is not None:
            twobit.close()

    def __enter__(self) -> Self:
        """Return ``self`` for use as a context manager."""
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        """Close the 2bit handle on context-manager exit."""
        self.close()

    @property
    def chrom_sizes(self) -> pd.Series:
        """Chromosome lengths as a pandas Series (a defensive copy).

        Integer lengths indexed by chromosome name, in reference order.
        """
        return self._chrom_sizes.copy()

    @property
    def chromosomes(self) -> list[str]:
        """Chromosome names, in the order the reference declares them."""
        return list(self._chrom_sizes.index)

    def fetch_sequence(self, region: str | Region) -> DNA:
        """Return the reference sequence for ``region`` as a :class:`~genome.seq.DNA`.

        Parameters
        ----------
        region : str or genome.region.Region
            Either a locus string ``chrom:start-end`` with **0-based, half-open**
            coordinates (``chr1:0-10`` is the first ten bases; thousands
            separators tolerated), a bare ``chrom`` for the whole sequence, or a
            :class:`~genome.region.Region`. When a ``Region`` carries strand
            ``"-"``, the reverse complement is returned.

        Returns
        -------
        genome.seq.DNA
            The sequence, with soft-masking case preserved. May contain ``N``
            runs where the reference is unknown.

        Raises
        ------
        ValueError
            If ``region`` is malformed, names an unknown chromosome, or its
            coordinates fall outside ``[0, chromosome length]``.

        Examples
        --------
        >>> genome = Genome("sacCer3")            # doctest: +SKIP
        >>> genome.fetch_sequence("chrIV:0-10")   # doctest: +SKIP
        DNA('ACACCACACC')
        """
        resolved = self._resolve_region(region)
        seq = DNA(self._twobit.nocheck_sequence(resolved.chrom, resolved.start, resolved.end))
        return seq.reverse_complement() if resolved.strand == "-" else seq

    def __getitem__(self, region: str | Region) -> DNA:
        """Index by locus string or :class:`~genome.region.Region` — sugar for :meth:`fetch_sequence`."""
        return self.fetch_sequence(region)

    def _resolve_region(self, region: str | Region) -> Region:
        """Validate ``region`` against the chrom sizes and return a concrete :class:`Region`.

        Accepts a 0-based locus string (a bare chromosome expands to the whole
        sequence) or an existing ``Region``. Raises :class:`ValueError` with an
        actionable message on an unknown chromosome or out-of-range coordinates.
        """
        if isinstance(region, Region):
            chrom, start, end, strand = region.chrom, region.start, region.end, region.strand
        else:
            chrom, start, end = parse_region(region)
            strand = "."

        if chrom not in self._chrom_sizes.index:
            known = ", ".join(str(name) for name in list(self._chrom_sizes.index)[:5])
            raise ValueError(
                f"unknown chromosome {chrom!r}; known sequences include: {known}, ... "
                f"(see Genome.chromosomes for the full list)."
            )
        size = int(self._chrom_sizes[chrom])

        if start is None or end is None:
            start, end = 0, size
        if start < 0:
            raise ValueError(f"region {region!s}: start must be >= 0 (0-based), got {start}.")
        if start > end:
            raise ValueError(f"region {region!s}: start ({start}) is past end ({end}).")
        if end > size:
            raise ValueError(
                f"region {region!s}: end ({end}) exceeds {chrom} length ({size}). "
                f"Coordinates are 0-based half-open, so the maximum valid end is {size}."
            )
        return Region(chrom, start, end, strand)
