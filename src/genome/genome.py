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
from genome.chimera import ChimeraNamingError, split_suffixed
from genome.io.chimera import ChimeraBuilder, ChimeraDetails, read_chimera_details
from genome.io.download import UCSCGenomeDownloader
from genome.io.fasta import GenomeFiles, read_chrom_sizes
from genome.io.gtf import AnnotationRegistry, BrokenAnnotation, GtfAnnotation
from genome.io.registration import AssemblyDir
from genome.io.twobit import TwoBit
from genome.metadata import AnnotationMetadata, AssemblyMetadata, assembly_metadata
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
        unknown when the table does not list ``assembly``, which is legal, since
        the table is a cross-reference rather than an allow-list.
    default_gtf : str, optional
        Name of the annotation to serve as :attr:`default_gtf`, overruling the one the
        annotation table flags. It need not be registered yet — see
        :attr:`default_gtf_path`.

    Attributes
    ----------
    assembly : str
        The assembly name.
    files : genome.io.fasta.GenomeFiles
        Paths to the prepared FASTA and its derived index/companion files.
    metadata : genome.metadata.AssemblyMetadata
        The assembly's metadata record, and always a record: the one passed in,
        else the curated table's row, else — for an assembly the table does not
        list — one carrying the name with every identifier unknown. Read a field
        off it (``genome.metadata.species``) without asking whether it is there.
        It is also what says where this assembly's FASTA is fetched from and which
        checksum it must match. Whether the table lists this assembly at all is a
        different question, and :func:`~genome.metadata.lookup_assembly`'s.

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
        # Total: an assembly the table does not list has a record whose fields are
        # unknown, never no record. Whether the table *lists* a name is a separate
        # question with a separate function, and it is what tells a chimera's derived
        # name from a free-form local key (ADR-0003).
        self.metadata: AssemblyMetadata = (
            metadata if metadata is not None else assembly_metadata(assembly)
        )
        # The override alone, and never the record above: given none the downloader
        # derives the same record for itself, so passing this one back would only hide
        # which of the two a caller supplied.
        self._downloader = UCSCGenomeDownloader(assembly, cache_dir, metadata=metadata)
        self._dir: AssemblyDir = self._downloader.dir
        self.files: GenomeFiles = (
            self._downloader.fetch_genome_from(path_or_url, progressbar=progressbar)
            if path_or_url is not None
            else self._downloader.fetch_genome(progressbar=progressbar)
        )
        self._chrom_sizes: pd.Series = read_chrom_sizes(self.files.chrom_sizes)
        self._twobit = TwoBit(self.files.twobit)
        # The record, never the metadata row, is what says an assembly is a chimera — and
        # it is already on disk by now, since the registration above wrote or confirmed
        # it. One small JSON read at open, and both accessors are then answered from
        # memory.
        self._chimera: ChimeraDetails | None = read_chimera_details(self._dir.path)
        # Reads the annotation table and lists the `gtf/` subtree both ways; acts on
        # none of it. Opening a genome must never start a registration — for a human
        # annotation that is a gigabyte download and a database build running many
        # minutes — and an annotation it cannot vouch for is recorded to report rather
        # than raised over, so one broken annotation never costs the genome.
        self._registry: AnnotationRegistry = AnnotationRegistry(
            self._dir, chrom_sizes=self.files.chrom_sizes, default=default_gtf
        )

    @classmethod
    def chimera(
        cls,
        *components: Genome,
        cache_dir: str | Path | None = None,
        force: bool = False,
    ) -> Self:
        """Build a **Chimera** of ``components`` and return it open, like any other genome.

        A second constructor rather than a second type (ADR-0008): the reference it
        produces is an assembly, so everything an assembly can do — fetch sequence,
        register an annotation, build an index — it can do, by one code path and not two.
        Its FASTA is its components' bytes with every chromosome name suffixed by the
        component it came from (ADR-0009); see :mod:`genome.io.chimera` for how it is
        written.

        The name is **derived**, never given — the component names sorted and joined by
        ``_``, so ``ce11`` and ``ecHT115`` in either order build and reopen the one
        ``ce11_ecHT115``. A chimera whose record says it finished is opened without
        rewriting anything; a directory that cannot be trusted raises and names
        ``genome register <name> --force``, which is what ``force`` is.

        Nothing is fetched, and no annotation argument is needed or accepted: each
        component carries its own :attr:`default_gtf`, so a caller's ``(assembly, gtf)``
        pairs split at the door and the annotation half travels with the components. Those
        defaults are **merged in the same act**, registered under the ``+``-join of their
        names in sorted-component order — so a built chimera arrives annotated and
        ``force`` repairs the annotation and the FASTA together. A component with no
        annotation contributes nothing, and components that contribute nothing between
        them leave the chimera with no annotation rather than an empty one; a component
        whose default is named but not registered here, or which has several registered
        and no default, raises before anything is written.

        Parameters
        ----------
        *components : Genome
            Two or more prepared component assemblies, in any order, each given once.
            None may itself be a chimera — a **Component** is always a canonical
            assembly, so nesting is forbidden by the model rather than deferred.
        cache_dir : str or pathlib.Path, optional
            Override the directory the chimera is built and opened in, exactly as the
            constructor's ``cache_dir`` does for an ordinary assembly. Defaults to the
            shared per-assembly reference directory, under the derived name.
        force : bool, default False
            Build again from scratch — the repair for a directory that raises.

        Returns
        -------
        Genome
            The chimera, opened under its derived name.

        Raises
        ------
        genome.chimera.ChimeraNamingError
            If fewer than two components are given, a component repeats, a component's
            name is not alphanumeric, a component is itself a chimera, or a component's
            FASTA carries a header that names no sequence for the suffix to ride on.
        genome.io.gtf.AnnotationNotRegisteredError
            If a component's default annotation is named but not registered here; the
            message names the command that registers it.
        genome.io.chimera.AmbiguousDefaultAnnotationError
            If a component carries several annotations and none is its default; the
            message names ``default_gtf=``.
        genome.io.completion.RegistrationError
            If the chimera's directory holds a build that cannot be trusted as finished,
            or the FASTA just built does not carry the sequences its components predict.
        genome.io.gtf.ChromosomeMismatchError
            If the merged annotation names a sequence the built FASTA does not carry.
        genome.external.ToolNotFoundError
            If ``samtools``, ``faToTwoBit`` or ``twoBitInfo`` are not on ``PATH``.

        Examples
        --------
        >>> worm, food = Genome("ce11"), Genome("ecHT115")     # doctest: +SKIP
        >>> chimera = Genome.chimera(worm, food)               # doctest: +SKIP
        >>> chimera.assembly                                   # doctest: +SKIP
        'ce11_ecHT115'
        >>> chimera["I__ce11:0-10"]                            # doctest: +SKIP
        DNA('GCCTAAGCCT')
        >>> chimera.default_gtf                                # doctest: +SKIP
        'wormbase_ws298+refseq_rs_2025_06_26'
        """
        builder = ChimeraBuilder(components, cache_dir)
        builder.build_genome(overwrite=force)
        return cls(builder.assembly, cache_dir=builder.cache_dir, progressbar=False)

    @property
    def default_gtf(self) -> str | None:
        """Name of this genome's **Default annotation**, or ``None`` when there is none.

        The ``default_gtf`` argument if one was given, else the annotation the table flags
        for this assembly, else the sole registered annotation, else ``None``. It names an
        annotation that may not be registered here, which on a fresh machine is the normal
        state and not an error — :attr:`default_gtf_path` is where it has to exist.
        """
        return self._registry.default

    @property
    def annotations(self) -> list[str]:
        """Names of the GTF annotations registered for this assembly **on this machine**.

        What is here, as against :attr:`offered_annotations`, which is what the lab
        supports, and :attr:`broken_annotations`, which is what is here and cannot be
        trusted.
        """
        return self._registry.registered

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
        return self._registry.broken

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
        return self._registry.offered

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
        return self._registry.register(
            name,
            force=force,
            progressbar=progressbar,
            metadata=metadata,
            check_chromosomes=check_chromosomes,
            disable_infer_genes=disable_infer_genes,
            disable_infer_transcripts=disable_infer_transcripts,
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
        :meth:`~genome.io.gtf.AnnotationRegistry.register_path` for the rest.
        """
        return self._registry.register_path(
            gtf,
            name,
            force=force,
            check_chromosomes=check_chromosomes,
            disable_infer_genes=disable_infer_genes,
            disable_infer_transcripts=disable_infer_transcripts,
        )

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
        return self._registry.path(name)

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
        default = self._registry.default
        return None if default is None else self._registry.path(default)

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
    def assembly_dir(self) -> AssemblyDir:
        """The **Assembly dir** this genome was opened in, and the layout inside it.

        Where everything tied to this assembly lives — its own files, the ``gtf/``
        subtree its annotations are filed under, the ``index/`` subtree its indexes are
        built into. Public because it is what an **Index** is derived from: an index
        belongs inside the assembly it indexes, and asking the genome is the only way to
        get the directory this genome was actually opened in rather than the one the
        **Data dir** layout would name for its assembly.

        Examples
        --------
        >>> sacCer3 = Genome("sacCer3")                   # doctest: +SKIP
        >>> sacCer3.assembly_dir.index_dir("chromap")     # doctest: +SKIP
        PosixPath('/data/genome/sacCer3/index/chromap')
        """
        return self._dir

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

    @property
    def components(self) -> list[str] | None:
        """The **Component** assembly names this is a **Chimera** of, or ``None``.

        The single test of whether an assembly is a chimera, and the completion record
        is what answers it — never the metadata table, which lists a chimera as a
        cross-reference and would answer the same question differently on a machine
        where the row is stale or absent (ADR-0008).

        Sorted, which for a chimera is the order its own derived name spells them in.
        ``None`` — not an empty list — for an ordinary assembly: it is not a chimera of
        nothing, it is not a chimera.

        Examples
        --------
        >>> chimera = Genome.chimera(Genome("ce11"), Genome("ecHT115"))  # doctest: +SKIP
        >>> chimera.components                                           # doctest: +SKIP
        ['ce11', 'ecHT115']
        >>> Genome("ce11").components is None                            # doctest: +SKIP
        True
        """
        return None if self._chimera is None else self._chimera.components

    @property
    def chrom_components(self) -> pd.Series:
        """Which assembly each chromosome came from, as a Series mirroring :attr:`chrom_sizes`.

        Attribution, and **total**: every chromosome this reference carries gets an
        answer. For a chimera each name is split at its recorded separator, so the answer
        is read out of the name itself and no mapping was ever stored (ADR-0009); for an
        assembly that is not a chimera every chromosome maps to that assembly's own name,
        which is true rather than merely convenient — and it is what leaves
        :attr:`components` as the single is-chimera test, since no caller has to read
        this one to find out.

        Returns
        -------
        pandas.Series
            Component assembly name per chromosome, indexed and ordered exactly as
            :attr:`chrom_sizes` is.

        Examples
        --------
        >>> chimera = Genome.chimera(Genome("ce11"), Genome("ecHT115"))  # doctest: +SKIP
        >>> chimera.chrom_components["I__ce11"]                          # doctest: +SKIP
        'ce11'
        >>> Genome("ce11").chrom_components["I"]                         # doctest: +SKIP
        'ce11'
        """
        if self._chimera is None:
            components = [self.assembly] * len(self._chrom_sizes)
        else:
            components = [
                split_suffixed(str(name), self._chimera.separator)[1]
                for name in self._chrom_sizes.index
            ]
        return pd.Series(components, index=self._chrom_sizes.index, name="component")

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
            coordinates fall outside ``[0, chromosome length]``. Against a chimera a
            bare chromosome name is one of the unknown ones (ADR-0009), and the message
            names the suffixed spellings that do resolve.

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
            raise self._unknown_chromosome(chrom)
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

    def _unknown_chromosome(self, chrom: str) -> ValueError:
        """Return the error a name this reference does not carry earns, unraised.

        Two messages, because a chimera has one more thing to say. A bare name that this
        chimera carries under one or more suffixed spellings is not merely unknown — it is
        the name of a real sequence, spelled the way a component spells it — and the
        refusal that ADR-0009 accepted the cost of is only bearable if it hands back the
        spellings that do resolve. Every other unknown name, on a chimera or not, gets the
        general message, which names a few sequences and where the rest are.

        A :class:`ValueError` in both cases: it is the type this raises today and the one
        callers catch.
        """
        spellings = self._suffixed_spellings(chrom)
        if spellings:
            listed = ", ".join(spellings)
            return ValueError(
                f"unknown chromosome {chrom!r}; {self.assembly} is a chimera, and every "
                f"chromosome name in one carries the component it came from, so a bare "
                f"name never resolves (ADR-0009). It carries {chrom!r} as: {listed}. Ask "
                f"for the one you meant."
            )
        known = ", ".join(str(name) for name in list(self._chrom_sizes.index)[:5])
        return ValueError(
            f"unknown chromosome {chrom!r}; known sequences include: {known}, ... "
            f"(see Genome.chromosomes for the full list)."
        )

    def _suffixed_spellings(self, chrom: str) -> list[str]:
        """Return this chimera's names for the bare chromosome ``chrom``, in reference order.

        Empty for an assembly that is not a chimera, and empty for a bare name none of its
        components contributed — so an ordinary unknown name reads the same either way.
        Computed here rather than kept, since nothing but a failed lookup ever asks.
        """
        if self._chimera is None:
            return []
        spellings: list[str] = []
        for name in self._chrom_sizes.index:
            spelled = str(name)
            try:
                bare, _component = split_suffixed(spelled, self._chimera.separator)
            except ChimeraNamingError:
                # A name no chimera build wrote — skipped rather than raised on, since
                # this is already the error path and a second failure helps nobody.
                continue
            if bare == chrom:
                spellings.append(spelled)
        return spellings
