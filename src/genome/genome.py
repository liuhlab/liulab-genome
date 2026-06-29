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

from genome.aligner.mixin import AlignerMixin
from genome.io.download import UCSCGenomeDownloader
from genome.io.fasta import GenomeFiles, read_chrom_sizes
from genome.io.gtf import GtfAnnotation, list_annotations, register_gtf
from genome.io.twobit import TwoBit
from genome.metadata import lookup_assembly
from genome.region import Region, parse_region
from genome.seq import DNA


class Genome(AlignerMixin):
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
        UCSC assembly name, e.g. ``"sacCer3"``, ``"hg38"``, ``"mm39"``. When
        ``path_or_url`` is omitted, the FASTA is downloaded from UCSC and the name
        is validated against UCSC first, so a typo fails fast. When
        ``path_or_url`` is given, the name only labels the cache directory and
        files; UCSC is not contacted.
    path_or_url : str or pathlib.Path, optional
        Seed the assembly from your own FASTA instead of downloading from UCSC —
        either a local file path (copied into the cache) or an http(s)/ftp URL
        (downloaded with ``curl``). Gzipped (``.gz``) sources are decompressed.
        Useful when UCSC is unreachable (firewall/proxy) or for a custom
        reference. See :meth:`~genome.io.download.UCSCGenomeDownloader.fetch_genome_from`.
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
        path_or_url: str | Path | None = None,
        cache_dir: str | Path | None = None,
        progressbar: bool = True,
        assembly_name: str | None = None,
        species: str | None = None,
        ucsc_name: str | None = None,
        ncbi_name: str | None = None,
        ncbi_assembly_id: str | None = None,
        ncbi_taxid: int | None = None,
        default_gtf: str | None = None,
    ) -> None:
        self.assembly = assembly
        self._set_metadata(
            assembly_name=assembly_name,
            species=species,
            ucsc_name=ucsc_name,
            ncbi_name=ncbi_name,
            ncbi_assembly_id=ncbi_assembly_id,
            ncbi_taxid=ncbi_taxid,
        )
        self._downloader = UCSCGenomeDownloader(assembly, cache_dir)
        self._assembly_dir: Path = self._downloader.cache_dir
        self.files: GenomeFiles = (
            self._downloader.fetch_genome_from(path_or_url, progressbar=progressbar)
            if path_or_url is not None
            else self._downloader.fetch_genome(progressbar=progressbar)
        )
        self._chrom_sizes: pd.Series = read_chrom_sizes(self.files.chrom_sizes)
        self._twobit = TwoBit(self.files.twobit)
        self._set_default_gtf(default_gtf)

    def _set_metadata(
        self,
        *,
        assembly_name: str | None,
        species: str | None,
        ucsc_name: str | None,
        ncbi_name: str | None,
        ncbi_assembly_id: str | None,
        ncbi_taxid: int | None,
    ) -> None:
        """Set metadata attributes, filling unset ones from the curated table when known."""
        table = lookup_assembly(self.assembly)
        self.assembly_name: str | None = assembly_name or (table.assembly_name if table else None)
        self.species: str | None = species or (table.species if table else None)
        self.ucsc_name: str | None = ucsc_name or (table.ucsc_name if table else None)
        self.ncbi_name: str | None = ncbi_name or (table.ncbi_name if table else None)
        self.ncbi_assembly_id: str | None = ncbi_assembly_id or (
            table.ncbi_assembly_id if table else None
        )
        self.ncbi_taxid: int | None = (
            ncbi_taxid if ncbi_taxid is not None else (table.ncbi_taxid if table else None)
        )

    def _set_default_gtf(self, default_gtf: str | None) -> None:
        """Discover registered annotations and pick the default GTF."""
        self._annotations: dict[str, GtfAnnotation] = list_annotations(self._assembly_dir)
        if default_gtf is not None:
            if default_gtf not in self._annotations:
                known = ", ".join(self._annotations) or "(none registered)"
                raise ValueError(
                    f"default_gtf {default_gtf!r} is not registered for {self.assembly!r}; "
                    f"registered annotations: {known}."
                )
            self.default_gtf: str | None = default_gtf
        elif len(self._annotations) == 1:
            self.default_gtf = next(iter(self._annotations))
        else:
            self.default_gtf = None

    @property
    def annotations(self) -> list[str]:
        """Names of the GTF annotations registered for this assembly."""
        return list(self._annotations)

    def register_gtf(
        self,
        gtf: str | Path,
        name: str,
        *,
        force: bool = False,
        disable_infer_genes: bool = True,
        disable_infer_transcripts: bool = True,
    ) -> GtfAnnotation:
        """Register a GTF under ``name`` and build its gffutils database.

        The GTF is placed under ``<assembly dir>/gtf/<name>/`` (a gzipped
        ``.gz`` source is decompressed automatically) and a gffutils database is
        built beside it. If no default GTF is set and this becomes the only
        annotation, it is adopted as :attr:`default_gtf`.
        """
        annotation = register_gtf(
            self._assembly_dir,
            gtf,
            name,
            force=force,
            disable_infer_genes=disable_infer_genes,
            disable_infer_transcripts=disable_infer_transcripts,
        )
        self._annotations[name] = annotation
        if self.default_gtf is None and len(self._annotations) == 1:
            self.default_gtf = name
        return annotation

    def get_gtf_path(self, name: str) -> Path:
        """Return the GTF file path of the annotation registered as ``name``."""
        if name not in self._annotations:
            known = ", ".join(self._annotations) or "(none registered)"
            raise KeyError(f"no annotation {name!r} for {self.assembly!r}; registered: {known}.")
        return self._annotations[name].gtf

    @property
    def default_gtf_path(self) -> Path | None:
        """GTF file path of the default annotation, or ``None`` when no default is set."""
        if self.default_gtf is None:
            return None
        return self.get_gtf_path(self.default_gtf)

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
