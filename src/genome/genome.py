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
from genome.io.gtf import (
    AnnotationNotRegisteredError,
    BrokenAnnotation,
    GtfAnnotation,
    default_annotation,
    fetch_annotation,
    list_annotations,
    list_broken_annotations,
    register_gtf,
)
from genome.io.twobit import TwoBit
from genome.metadata import (
    AnnotationMetadata,
    AssemblyMetadata,
    list_annotation_metadata,
    lookup_assembly,
)
from genome.region import Region, parse_region
from genome.seq import DNA


class Genome(AlignerMixin):
    """A reference genome and the operations over it.

    Constructing a ``Genome`` ensures the assembly's reference files exist
    locally: the FASTA is fetched from the **Source** its metadata row pins —
    UCSC for most assemblies, WormBase or NCBI for others — or, for an assembly
    no row lists, from a URL derived from the UCSC golden path. Its ``.fai``
    index, ``.2bit`` encoding and ``chrom.sizes`` are then prepared, and a
    completion record is written recording what was done. Everything lands under
    ``<LIULAB_DATA>/genome/<assembly>/`` (see
    :func:`~genome.io.download.assembly_data_dir`). A later construction reads
    that record, confirms every file it claims is present at the size it claims,
    and opens the ``.2bit`` — so it is instant and works offline. Nothing is
    downloaded twice.

    Sequence is read from the ``.2bit`` file via ``py2bit``; coordinates are
    **0-based, half-open** throughout.

    Parameters
    ----------
    assembly : str
        Assembly name, e.g. ``"sacCer3"``, ``"hg38"``, ``"mm39"``. A free-form
        local key (ADR-0003), not necessarily a UCSC one — ``"ecHT115"`` is a
        reference UCSC has never carried. When ``path_or_url`` is omitted and no
        row pins a source, the FASTA is downloaded from UCSC and the name is
        validated against UCSC first, so a typo fails fast; a pinned source *is*
        the source, so there is nothing to guess and that check is skipped. When
        ``path_or_url`` is given, the name only labels the cache directory and
        files; UCSC is not contacted.
    path_or_url : str or pathlib.Path, optional
        Seed the assembly from your own FASTA instead of downloading from UCSC —
        either a local file path (copied into the cache) or an http(s)/ftp/sftp
        URL (fetched with pooch). Gzipped (``.gz``) sources are decompressed.
        Useful when UCSC is unreachable (firewall/proxy) or for a custom
        reference. See :meth:`~genome.io.download.UCSCGenomeDownloader.fetch_genome_from`.
    cache_dir : str or pathlib.Path, optional
        Override the storage directory for this assembly's files. Defaults to the
        shared per-assembly reference directory.
    progressbar : bool, default True
        Show a download progress bar on first fetch (requires ``tqdm``).
    metadata : genome.metadata.AssemblyMetadata, optional
        A complete metadata record, used *instead of* the curated table's row for
        ``assembly``. All-or-nothing: pass a record and every identifier comes
        from it; omit it and every identifier comes from the table — or is
        ``None`` when the table does not list ``assembly``, which is legal, since
        the table is a cross-reference rather than an allow-list.
    default_gtf : str, optional
        Name of the annotation to serve as :attr:`default_gtf`, overruling the one the
        annotation table flags. It need not be registered yet — see
        :attr:`default_gtf_path`.

    Attributes
    ----------
    assembly : str
        The assembly name.
    default_gtf : str or None
        Name of this genome's **Default annotation**: the ``default_gtf`` argument if
        one was given, else the annotation the table flags for this assembly, else the
        sole registered annotation, else ``None``. It names an annotation that may not
        be registered here, which on a fresh machine is the normal state and not an
        error.
    files : genome.io.fasta.GenomeFiles
        Paths to the prepared FASTA and its derived index/companion files.
    metadata : genome.metadata.AssemblyMetadata or None
        The assembly's metadata record — the one passed in, else the curated
        table's row, else ``None`` for an assembly the table does not list. It is
        also what says where this assembly's FASTA is fetched from and which
        checksum it must match. Its fields are read directly off the genome too, as
        :attr:`assembly_name`, :attr:`species`, :attr:`ucsc_name`,
        :attr:`ncbi_name`, :attr:`ncbi_assembly_id`, :attr:`ncbi_taxid`,
        :attr:`source_url` and :attr:`sha256`.

    Raises
    ------
    ValueError
        If ``assembly`` is unknown to UCSC.
    genome.io.completion.RegistrationError
        If the assembly's directory holds a registration that cannot be trusted —
        files with no record (an interrupted run), or a record that disagrees with
        what is on disk. The message names the file and
        ``genome register <assembly> --force``, which repairs it (ADR-0007). An
        absent or empty directory is not this: that is a fresh registration.
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
        metadata: AssemblyMetadata | None = None,
        default_gtf: str | None = None,
    ) -> None:
        self.assembly = assembly
        self.metadata: AssemblyMetadata | None = (
            metadata if metadata is not None else lookup_assembly(assembly)
        )
        self._downloader = UCSCGenomeDownloader(assembly, cache_dir, metadata=self.metadata)
        self._assembly_dir: Path = self._downloader.cache_dir
        self.files: GenomeFiles = (
            self._downloader.fetch_genome_from(path_or_url, progressbar=progressbar)
            if path_or_url is not None
            else self._downloader.fetch_genome(progressbar=progressbar)
        )
        self._chrom_sizes: pd.Series = read_chrom_sizes(self.files.chrom_sizes)
        self._twobit = TwoBit(self.files.twobit)
        self._set_default_gtf(default_gtf)

    @property
    def assembly_name(self) -> str | None:
        """Canonical name of the assembly, or ``None`` when its metadata is unknown."""
        return self.metadata.assembly_name if self.metadata else None

    @property
    def species(self) -> str | None:
        """Species this assembly is a reference for, or ``None`` when unknown."""
        return self.metadata.species if self.metadata else None

    @property
    def ucsc_name(self) -> str | None:
        """UCSC's name for the assembly, or ``None`` when it has none.

        ``None`` covers both an assembly the metadata table does not list and one it
        lists that UCSC has never carried — the assembly id is a local key and UCSC is
        only the default source, so a reference can be fully supported here and have no
        name in that namespace at all (ADR-0003).
        """
        return self.metadata.ucsc_name if self.metadata else None

    @property
    def ncbi_name(self) -> str | None:
        """NCBI's name for the assembly (e.g. ``"GRCh38"``), or ``None`` when unknown."""
        return self.metadata.ncbi_name if self.metadata else None

    @property
    def ncbi_assembly_id(self) -> str | None:
        """NCBI assembly accession (e.g. ``"GCF_000001405.40"``), or ``None`` when unknown."""
        return self.metadata.ncbi_assembly_id if self.metadata else None

    @property
    def ncbi_taxid(self) -> int | None:
        """NCBI taxonomy id of the species, or ``None`` when unknown."""
        return self.metadata.ncbi_taxid if self.metadata else None

    @property
    def source_url(self) -> str | None:
        """URL this assembly's FASTA is pinned to, or ``None`` when nothing is pinned."""
        return self.metadata.source_url if self.metadata else None

    @property
    def sha256(self) -> str | None:
        """Pinned sha256 of the *unpacked* FASTA, or ``None`` when nothing is pinned.

        The value the metadata records and the download is checked against — not a
        digest of the files on disk, which :meth:`verify_fasta
        <genome.io.download.UCSCGenomeDownloader.verify_fasta>` computes.
        """
        return self.metadata.sha256 if self.metadata else None

    def _set_default_gtf(self, default_gtf: str | None) -> None:
        """Read the annotation lists and settle which one is the default.

        All three are read, none is acted on: the table is looked up and the ``gtf/``
        subtree is listed both ways, and nothing is fetched, built or created. Opening a
        genome must never start a registration — for a human annotation that is a
        gigabyte download and a database build running many minutes — and an annotation
        it cannot vouch for is recorded to report rather than raised over, so one broken
        annotation never costs the genome.
        """
        self._annotations: dict[str, GtfAnnotation] = list_annotations(self._assembly_dir)
        self._broken: dict[str, BrokenAnnotation] = list_broken_annotations(
            self._assembly_dir, self.assembly
        )
        self._offered: list[AnnotationMetadata] = list_annotation_metadata(self.assembly)
        self.default_gtf: str | None = default_annotation(
            self._offered, self._annotations, explicit=default_gtf
        )

    @property
    def annotations(self) -> list[str]:
        """Names of the GTF annotations registered for this assembly **on this machine**.

        What is here, as against :attr:`offered_annotations`, which is what the lab
        supports, and :attr:`broken_annotations`, which is what is here and cannot be
        trusted.
        """
        return list(self._annotations)

    @property
    def broken_annotations(self) -> list[BrokenAnnotation]:
        """The annotation directories here that cannot be trusted as finished.

        What :attr:`annotations` leaves out. A build killed part-way leaves a database
        with most of the genes missing, so it is deliberately not registered — but it is
        also not nothing, and a caller who never re-registers it would otherwise never
        hear that it is there. Each entry says what is wrong and names the one command
        that repairs it.

        Reading this raises nothing, whatever state the directory is in: it was settled
        when the genome opened, which is why one broken annotation never stopped it.

        Examples
        --------
        >>> sacCer3 = Genome("sacCer3")                              # doctest: +SKIP
        >>> [broken.name for broken in sacCer3.broken_annotations]   # doctest: +SKIP
        ['ensgene_v101']
        """
        return list(self._broken.values())

    @property
    def offered_annotations(self) -> list[AnnotationMetadata]:
        """The annotations the curated table offers for this assembly, in table order.

        What the lab supports for this assembly — each row saying who publishes it,
        which release it is, where it is fetched from and what it must hash to — as
        against :attr:`annotations`, which is what is registered here. A row appears
        whether or not anyone has registered it, and registering one is
        :meth:`register_annotation`; nothing here reads the disk. Empty for an assembly
        the table offers nothing for, which is legal: the table is a cross-reference
        rather than an allow-list (ADR-0003).

        A fresh list each call, so a caller may sort or filter it.

        Examples
        --------
        >>> sacCer3 = Genome("sacCer3")                            # doctest: +SKIP
        >>> [record.name for record in sacCer3.offered_annotations]  # doctest: +SKIP
        ['ensgene_v101']
        """
        return list(self._offered)

    def register_annotation(
        self,
        name: str,
        *,
        force: bool = False,
        progressbar: bool = True,
        metadata: AnnotationMetadata | None = None,
        check_chromosomes: bool = True,
        disable_infer_genes: bool = True,
        disable_infer_transcripts: bool = True,
    ) -> GtfAnnotation:
        """Register the annotation this assembly's table row lists as ``name``.

        Naming it is enough: the curated annotation table says where the GTF comes
        from and what it must hash to. It is fetched, verified against that digest,
        checked against this assembly's chromosome names, placed under
        ``<assembly dir>/gtf/<name>/``, built into a gffutils database and recorded. If
        no default GTF is set and this becomes the only annotation, it is adopted as
        :attr:`default_gtf`.

        One that is already registered is returned silently — nothing is fetched and
        nothing is rebuilt — while a directory that cannot be trusted raises and names
        its repair. See :func:`~genome.io.gtf.fetch_annotation`.

        Parameters
        ----------
        name : str
            The **Registered name** the table lists for this assembly.
        force : bool, default False
            Register again from scratch — the repair for a directory that raises.
        progressbar : bool, default True
            Show a download progress bar (requires ``tqdm``).
        metadata : genome.metadata.AnnotationMetadata, optional
            A complete annotation record, used *instead of* the curated table's row for
            ``name`` — the same all-or-nothing override the constructor takes for the
            assembly's own metadata.
        check_chromosomes : bool, default True
            Refuse a GTF naming sequences this assembly does not carry, before paying
            for the database build. Pass ``False`` to register one whose mismatch you
            have inspected and accept.
        disable_infer_genes : bool, default True
            Do not reconstruct ``gene`` features from exon lines.
        disable_infer_transcripts : bool, default True
            Do not reconstruct ``transcript`` features from exon lines.

        Returns
        -------
        genome.io.gtf.GtfAnnotation
            The registered annotation's name and its two file paths.

        Raises
        ------
        ValueError
            If the table lists no annotation ``name`` for this assembly.
        genome.io.gtf.ChromosomeMismatchError
            If the GTF names sequences this assembly does not carry.
        genome.io.utils.ChecksumMismatchError
            If the fetched GTF is not the digest the row pins.
        genome.io.completion.RegistrationError
            If the annotation's directory cannot be trusted as finished.

        Examples
        --------
        >>> sacCer3 = Genome("sacCer3")                        # doctest: +SKIP
        >>> sacCer3.register_annotation("ensgene_v101")        # doctest: +SKIP
        GtfAnnotation(name='ensgene_v101', ...)
        """
        return self._adopt(
            fetch_annotation(
                self._assembly_dir,
                self.assembly,
                name,
                force=force,
                progressbar=progressbar,
                metadata=metadata,
                check_chromosomes=check_chromosomes,
                disable_infer_genes=disable_infer_genes,
                disable_infer_transcripts=disable_infer_transcripts,
            )
        )

    def register_gtf(
        self,
        gtf: str | Path,
        name: str,
        *,
        force: bool = False,
        check_chromosomes: bool = True,
        disable_infer_genes: bool = True,
        disable_infer_transcripts: bool = True,
    ) -> GtfAnnotation:
        """Register the GTF at ``gtf`` under ``name`` and build its gffutils database.

        The escape hatch for an annotation the curated table does not list —
        :meth:`register_annotation` is the way in for one it does. Its chromosome names
        are checked against this genome's own before anything is created, the GTF is
        placed under ``<assembly dir>/gtf/<name>/`` (a gzipped ``.gz`` source is
        decompressed automatically), a gffutils database is built beside it, and the
        record that says so is written last. If no default GTF is set and this becomes
        the only annotation, it is adopted as :attr:`default_gtf`. ``check_chromosomes``
        is as it is on :meth:`register_annotation`; see
        :func:`~genome.io.gtf.register_gtf` for the rest.
        """
        return self._adopt(
            register_gtf(
                self._assembly_dir,
                gtf,
                name,
                force=force,
                chrom_sizes=self.chrom_sizes_path,
                check_chromosomes=check_chromosomes,
                disable_infer_genes=disable_infer_genes,
                disable_infer_transcripts=disable_infer_transcripts,
            )
        )

    def _adopt(self, annotation: GtfAnnotation) -> GtfAnnotation:
        """Add a freshly registered annotation to the registry, adopting it if it is alone.

        The sole-registered clause of the default rule, applied the moment it becomes
        true. A default already decided — the caller's choice, or the table's flag —
        is never displaced by one being registered.

        Registering over a broken directory is what repairs it, so the name stops being
        reported as broken here rather than only on the next open.
        """
        self._annotations[annotation.name] = annotation
        self._broken.pop(annotation.name, None)
        if self.default_gtf is None and len(self._annotations) == 1:
            self.default_gtf = annotation.name
        return annotation

    def get_gtf_path(self, name: str) -> Path:
        """Return the GTF file path of the annotation registered as ``name``.

        Parameters
        ----------
        name : str
            The **Registered name** to resolve.

        Returns
        -------
        pathlib.Path
            Path to the placed ``<name>.gtf``.

        Raises
        ------
        genome.io.gtf.AnnotationNotRegisteredError
            If nothing of that name is registered here. The message says what is
            registered, and names the command that closes the gap: the one that
            registers ``name`` when the table offers it, the path-based way in when it
            does not, and — when a directory of that name is there but broken — the one
            that registers it again from scratch, so the command named is one that runs
            rather than one that raises in turn. It is a :class:`KeyError`.

        Examples
        --------
        >>> sacCer3 = Genome("sacCer3")                     # doctest: +SKIP
        >>> sacCer3.get_gtf_path("ensgene_v101")            # doctest: +SKIP
        PosixPath('/data/genome/sacCer3/gtf/ensgene_v101/ensgene_v101.gtf')
        """
        if name not in self._annotations:
            raise AnnotationNotRegisteredError(
                self.assembly,
                name,
                self._annotations,
                [record.name for record in self._offered],
                broken=self._broken.get(name),
            )
        return self._annotations[name].gtf

    @property
    def default_gtf_path(self) -> Path | None:
        """GTF file path of the **Default annotation**, or ``None`` when there is no default.

        Where the default stops being an intention and has to exist. :attr:`default_gtf`
        may name an annotation nobody has registered on this machine — the table's
        choice, on a machine that has not fetched it yet, or one a caller named at
        construction ahead of registering it — and asking for its path is what says so,
        naming the command that registers it, or the one that repairs it when a broken
        directory of that name is there. ``None`` means no default was decided at all,
        which is a different answer from one that is not registered.

        Raises
        ------
        genome.io.gtf.AnnotationNotRegisteredError
            If the default annotation is not registered here.

        Examples
        --------
        >>> sacCer3 = Genome("sacCer3")                     # doctest: +SKIP
        >>> sacCer3.default_gtf_path                        # doctest: +SKIP
        PosixPath('/data/genome/sacCer3/gtf/ensgene_v101/ensgene_v101.gtf')
        """
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
    def fasta_path(self) -> Path:
        """Path to the reference FASTA file."""
        return self.files.fasta

    @property
    def twobit_path(self) -> Path:
        """Path to the ``.2bit`` encoding of the reference."""
        return self.files.twobit

    @property
    def chrom_sizes_path(self) -> Path:
        r"""Path to the ``chrom.sizes`` file (``<name>\t<length>`` per sequence)."""
        return self.files.chrom_sizes

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
