"""Annotation registration — place a GTF, build its database, and record that it finished.

I/O boundary module. A reference assembly may carry several gene annotations (GENCODE,
RefSeq, WormBase, …). Each is registered under a **Registered name** and lives in its own
directory beside the assembly's sequence files::

    <LIULAB_DATA>/genome/<assembly>/gtf/<name>/
        <name>.gtf          # the annotation, kept decompressed
        <name>.db           # the gffutils SQLite database built from it
        .completion.json    # the record saying all of that finished
        .work/              # the disposable working area a fetch downloads into

There are two ways in. :func:`fetch_annotation` takes the **name** the curated annotation
table lists for this assembly, fetches that row's URL, checks the unpacked GTF against
the sha256 the row pins (ADR-0006), builds the database and writes the record;
:func:`register_annotation` is the same thing addressed by assembly name, returning the
record itself, and is what the CLI drives. :func:`register_gtf` takes a **path** instead
and is the escape hatch for a GTF no row lists.

**A GTF belongs to its assembly or to nothing.** Either way in checks that every
**Chromosome** the GTF names is one the assembly's ``chrom.sizes`` carries, and raises
:class:`ChromosomeMismatchError` when it is not — the Ensembl-versus-UCSC spelling
(``1`` against ``chr1``) is what this catches, and an annotation whose every feature sits
on a sequence the assembly never heard of is worse than no annotation at all. The check
is strict one way only, since an assembly may carry scaffolds the annotation never
mentions; it streams the GTF; and it runs *before* the database build and before the GTF
is placed, so a mismatch costs seconds and leaves nothing behind.

**A record is the only thing that says an annotation is registered.** A database file's
mere existence never is — a `gffutils` build killed half-way leaves a partial database
that answers queries with most of the genes missing, which is exactly what the
**Completion marker** exists to distrust. So :func:`list_annotations` reports what has a
record that agrees with disk, re-registering something that already has one returns it
silently, and a directory holding files but no record raises and names its repair
(ADR-0007).

Examples
--------
>>> from pathlib import Path
>>> from genome.io.gtf import annotation_dir
>>> annotation_dir(Path("/data/genome/sacCer3"), "ensgene_v101").name
'ensgene_v101'
"""

from __future__ import annotations

import gzip
import shutil
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import IO, Any

import gffutils
import pooch

from genome.io import download
from genome.io.completion import (
    RECORD_NAME,
    RegistrationError,
    build_record,
    check_registration,
    clear_work_dir,
    disagreements,
    read_record,
    work_dir,
    write_record,
)
from genome.io.fasta import read_chrom_sizes
from genome.io.utils import ChecksumMismatchError, _gunzip, sha256_file
from genome.metadata import AnnotationMetadata, list_annotation_metadata, lookup_annotation

#: Subdirectory under an assembly's data dir holding all its GTF annotations.
_GTF_SUBDIR = "gtf"

#: How many names an error lists before saying how many it left out. A whole-genome
#: mismatch offends in the thousands, and a message that long is one nobody reads.
_MAX_LISTED_NAMES = 10


class ChromosomeMismatchError(ValueError):
    """A GTF names **Chromosome**s its assembly does not carry, so the two do not line up.

    Registering it would build an annotation where nothing matches: every feature would
    sit on a sequence the assembly has never heard of, and every query over it would
    answer nothing while looking perfectly healthy. The usual cause is a spelling
    difference rather than a wrong file, so the message says which one out loud and
    names the argument that registers it anyway.

    The check behind it is strict in one direction only. An assembly carrying scaffolds
    the annotation never mentions is normal and is not this.

    Parameters
    ----------
    name : str
        The **Registered name** the annotation was being registered under.
    missing : iterable of str
        Every name the GTF uses that the assembly's ``chrom.sizes`` does not list.
    known : iterable of str
        The names the assembly does carry, for the message to contrast against.

    Attributes
    ----------
    name : str
        The registered name.
    missing : tuple of str
        **Every** offending name, sorted — the message lists at most ten of them and
        counts the rest, this is the whole set.
    known : tuple of str
        The names the assembly carries, as they were passed in.

    Examples
    --------
    >>> raise ChromosomeMismatchError("gencode_v44", ["1", "2"], ["chr1", "chr2"])
    Traceback (most recent call last):
    genome.io.gtf.ChromosomeMismatchError: the GTF for 'gencode_v44' names 2 ...
    """

    def __init__(self, name: str, missing: Iterable[str], known: Iterable[str]) -> None:
        self.name = name
        self.missing: tuple[str, ...] = tuple(sorted(missing))
        self.known: tuple[str, ...] = tuple(known)
        count = len(self.missing)
        super().__init__(
            f"the GTF for {name!r} names {count} chromosome{'' if count == 1 else 's'} the "
            f"assembly does not carry: {_elide(self.missing)}. An annotation and its assembly "
            f"must spell chromosomes the same way, and the usual cause is a UCSC-versus-Ensembl "
            f"mismatch ('chr1' against '1', 'chrM' against 'MtDNA'). The assembly carries: "
            f"{_elide(self.known)}. Register the annotation built for this assembly, or pass "
            f"check_chromosomes=False to register this one anyway."
        )


@dataclass(frozen=True)
class GtfAnnotation:
    """A registered GTF annotation: its name and the on-disk GTF + database paths."""

    name: str
    gtf: Path
    db: Path


def _annotations_root(assembly_dir: Path) -> Path:
    """Return ``<assembly_dir>/gtf``, the parent of every annotation directory."""
    return assembly_dir / _GTF_SUBDIR


def annotation_dir(assembly_dir: Path, name: str) -> Path:
    """Return the directory holding the annotation registered as ``name``."""
    return _annotations_root(assembly_dir) / name


def _annotation_files(assembly_dir: Path, name: str) -> GtfAnnotation:
    """Resolve the GTF + database paths for ``name`` (without checking existence)."""
    directory = annotation_dir(assembly_dir, name)
    return GtfAnnotation(name=name, gtf=directory / f"{name}.gtf", db=directory / f"{name}.db")


def _repair_command(assembly: str, name: str) -> str:
    """Return the command that registers ``name`` again from scratch.

    Quoted verbatim into every error a broken annotation directory raises, so it has to
    be a command that exists and does the job.
    """
    return f"genome register-annotation {assembly} {name} --force"


def _path_repair_command(source: Path, name: str) -> str:
    """Return the call that registers an unlisted GTF again from scratch.

    A GTF no row lists cannot be repaired by name — nothing knows where it came from
    but the caller who passed the path — so the repair names the call, with the path.
    """
    return f"register_gtf(<assembly dir>, {str(source)!r}, {name!r}, force=True)"


def list_annotations(assembly_dir: Path) -> dict[str, GtfAnnotation]:
    """Return the annotations registered under ``assembly_dir``, keyed by name.

    Registered means *a record is there and agrees with what is on disk* — never that a
    database file exists, which is true of a build killed half-way through as well as of
    a finished one. Anything else in the ``gtf/`` subtree is left out rather than raised
    over: listing is a question about this machine, and one unfinished annotation must
    not stop a genome from opening. Registering that name again is what reports it as
    broken.

    Parameters
    ----------
    assembly_dir : pathlib.Path
        The assembly directory whose ``gtf/`` subtree is listed.

    Returns
    -------
    dict of str to GtfAnnotation
        Registered name to its files, in directory-name order. Empty when nothing is
        registered.

    Examples
    --------
    >>> from pathlib import Path
    >>> list_annotations(Path("/tmp/definitely-not-an-assembly"))
    {}
    """
    root = _annotations_root(assembly_dir)
    if not root.is_dir():
        return {}
    found: dict[str, GtfAnnotation] = {}
    for directory in sorted(root.iterdir()):
        if not directory.is_dir():
            continue
        record = read_record(directory)
        if record is None or disagreements(directory, record):
            continue
        annotation = _annotation_files(assembly_dir, directory.name)
        found[annotation.name] = annotation
    return found


def register_gtf(
    assembly_dir: Path,
    gtf: str | Path,
    name: str,
    *,
    force: bool = False,
    chrom_sizes: str | Path | None = None,
    check_chromosomes: bool = True,
    disable_infer_genes: bool = True,
    disable_infer_transcripts: bool = True,
) -> GtfAnnotation:
    """Register the GTF at ``gtf`` under ``name`` and build its gffutils database.

    The escape hatch for an annotation the curated table does not list: the caller says
    where the file is, and it is placed, built and recorded exactly as a listed one is.
    A gzipped (``.gz``) source is decompressed into the registered ``<name>.gtf``; a
    plain GTF is copied as-is. The digest recorded is of the placed GTF — nothing is
    compared against, because an unlisted annotation has no pinned digest to compare to.

    The source's chromosome names are checked against ``chrom_sizes`` before anything is
    created, so a GTF that does not line up with the assembly leaves the annotation
    directory exactly as it was found. Unlike :func:`fetch_annotation`, this form is
    given no assembly name and so cannot derive where that file is — pass it, or accept
    that nothing is checked. What was decided is written into the record either way, as
    ``details["chromosomes_checked"]``.

    Gene/transcript inference is disabled by default — standard annotation GTFs
    (GENCODE, Ensembl, RefSeq) already declare ``gene``/``transcript`` features, and
    inferring them is the classic gffutils slow path. Enable it only for a bare
    exon-level GTF.

    Parameters
    ----------
    assembly_dir : pathlib.Path
        The assembly directory this annotation is filed under.
    gtf : str or pathlib.Path
        Path to the source GTF, plain or ``.gz``.
    name : str
        The **Registered name** to address it by, unique within the assembly.
    force : bool, default False
        Register again from scratch — the repair for a directory that raises.
    chrom_sizes : str or pathlib.Path, optional
        The assembly's ``chrom.sizes``, whose names the GTF's must be among. Omitted —
        or pointing at a file that is not there — there is nothing to check against and
        the names are not checked.
    check_chromosomes : bool, default True
        Check the GTF's chromosome names against the assembly's. Pass ``False`` to
        register a GTF whose mismatch you have inspected and accept.
    disable_infer_genes : bool, default True
        Do not reconstruct ``gene`` features from exon lines.
    disable_infer_transcripts : bool, default True
        Do not reconstruct ``transcript`` features from exon lines.

    Returns
    -------
    GtfAnnotation
        The registered annotation's name and its two file paths.

    Raises
    ------
    FileNotFoundError
        If ``gtf`` is not a file.
    ChromosomeMismatchError
        If the GTF names sequences the assembly does not carry; the message lists them
        and names the usual cause.
    genome.io.completion.UnfinishedRegistrationError
        If the annotation's directory holds files but no record.
    genome.io.completion.RegistrationMismatchError
        If its record disagrees with what is on disk.

    Examples
    --------
    >>> from pathlib import Path
    >>> register_gtf(                                    # doctest: +SKIP
    ...     Path("/data/genome/sacCer3"), "custom.gtf.gz", "custom"
    ... )
    GtfAnnotation(name='custom', ...)
    """
    source = Path(gtf)
    if not source.is_file():
        raise FileNotFoundError(
            f"GTF file not found: {source}. Pass the path of an existing .gtf or "
            f".gtf.gz, or register a listed annotation by name instead."
        )

    annotation = _annotation_files(assembly_dir, name)
    directory = annotation_dir(assembly_dir, name)
    if _already_registered(directory, force=force, repair=_path_repair_command(source, name)):
        return annotation

    known = (
        _assembly_chromosomes(Path(chrom_sizes))
        if check_chromosomes and chrom_sizes is not None
        else None
    )
    # Against the source, before the directory is even created: a mismatch then leaves
    # nothing behind, so the next call reports the same problem rather than the files
    # of an interrupted registration. A .gz source is streamed twice — once here and
    # once to place it — which buys that, and costs a fraction of the database build.
    _reject_unknown_chromosomes(source, known, name=name)

    directory.mkdir(parents=True, exist_ok=True)
    # A gzipped source is stream-decompressed into the registered <name>.gtf;
    # a plain GTF is copied as-is (skipping a copy onto itself).
    if source.suffix == ".gz":
        _gunzip(source, annotation.gtf)
    elif source.resolve() != annotation.gtf.resolve():
        shutil.copy2(source, annotation.gtf)

    return _build_and_record(
        annotation,
        source_url=str(source),
        sha256=sha256_file(annotation.gtf),
        details={"chromosomes_checked": known is not None},
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )


def fetch_annotation(
    assembly_dir: Path,
    assembly: str,
    name: str,
    *,
    force: bool = False,
    progressbar: bool = True,
    metadata: AnnotationMetadata | None = None,
    check_chromosomes: bool = True,
    disable_infer_genes: bool = True,
    disable_infer_transcripts: bool = True,
) -> GtfAnnotation:
    """Register the annotation the table lists for ``assembly`` as ``name``.

    Naming an annotation is enough: where its GTF comes from and which digest it must
    match are the curated table's to know. The row's URL is fetched into the working
    area, the **unpacked** GTF is checked against the sha256 the row pins (ADR-0006) —
    so a GTF that is not the pinned one never reaches the annotation directory — the
    gffutils database is built, and the record is written last.

    Its chromosome names are checked too, against ``<assembly_dir>/<assembly>.chrom.sizes``
    and while the GTF is still in the working area: every name the GTF uses must be one
    the assembly carries, so an Ensembl-spelled GTF registered against a UCSC-spelled
    assembly fails in seconds rather than after the minutes the database build takes.
    The reverse is not required — an assembly may carry scaffolds the annotation never
    mentions. An assembly with no ``chrom.sizes`` yet has nothing to check against, and
    the record says so in ``details["chromosomes_checked"]``.

    An annotation that already has a valid record is returned silently: nothing is
    fetched, nothing is rebuilt and nothing is warned about. A directory that cannot be
    trusted — files with no record, or a record that disagrees with disk — **raises**,
    naming ``genome register-annotation <assembly> <name> --force`` (ADR-0007). That is
    what ``force=True`` is: it skips the question, keeps a GTF whose digest can be shown
    to be the pinned one, and fetches the source again when it cannot.

    Parameters
    ----------
    assembly_dir : pathlib.Path
        The assembly directory this annotation is filed under.
    assembly : str
        The assembly whose table row is looked up, e.g. ``"hg38"``.
    name : str
        The **Registered name** the table lists, e.g. ``"gencode_v50"``.
    force : bool, default False
        Register again from scratch, repairing a directory that raises.
    progressbar : bool, default True
        Show a download progress bar (requires ``tqdm``).
    metadata : genome.metadata.AnnotationMetadata, optional
        A complete annotation record to use *instead of* the curated table's row. Omit
        it and the row is looked up here.
    check_chromosomes : bool, default True
        Check the GTF's chromosome names against the assembly's. Pass ``False`` to
        register an annotation whose mismatch you have inspected and accept.
    disable_infer_genes : bool, default True
        Do not reconstruct ``gene`` features from exon lines.
    disable_infer_transcripts : bool, default True
        Do not reconstruct ``transcript`` features from exon lines.

    Returns
    -------
    GtfAnnotation
        The registered annotation's name and its two file paths.

    Raises
    ------
    ValueError
        If the table lists no annotation ``name`` for ``assembly``; the message lists
        what it does offer and points at the path-based form for an unlisted GTF.
    ChromosomeMismatchError
        If the GTF names sequences the assembly does not carry; the message lists them
        and names the usual cause.
    genome.io.utils.ChecksumMismatchError
        If the row pins a sha256 and the unpacked GTF is not it; the message names both
        digests.
    genome.io.completion.UnfinishedRegistrationError
        If the annotation's directory holds files but no record.
    genome.io.completion.RegistrationMismatchError
        If its record disagrees with what is on disk.

    Examples
    --------
    >>> from pathlib import Path
    >>> fetch_annotation(                                # doctest: +SKIP
    ...     Path("/data/genome/sacCer3"), "sacCer3", "ensgene_v101"
    ... )
    GtfAnnotation(name='ensgene_v101', ...)
    """
    row = metadata if metadata is not None else lookup_annotation(assembly, name)
    if row is None:
        offered = ", ".join(r.name for r in list_annotation_metadata(assembly)) or "(none)"
        raise ValueError(
            f"no annotation named {name!r} is listed for {assembly!r}. Listed for it: "
            f"{offered}. An annotation the table does not list is registered by path "
            f"instead — Genome.register_gtf(<path>, {name!r})."
        )

    annotation = _annotation_files(assembly_dir, name)
    directory = annotation_dir(assembly_dir, name)
    if _already_registered(directory, force=force, repair=_repair_command(assembly, name)):
        return annotation

    known = (
        _assembly_chromosomes(assembly_dir / f"{assembly}.chrom.sizes")
        if check_chromosomes
        else None
    )
    digest = _proven_gtf(annotation.gtf, row.sha256)
    if digest is None:
        digest = _fetch_gtf(annotation, row, progressbar=progressbar, known=known)
    else:
        # Kept from a previous run, so it is already placed; there is nothing to undo.
        _reject_unknown_chromosomes(annotation.gtf, known, name=name)

    return _build_and_record(
        annotation,
        source_url=row.url,
        sha256=digest,
        details={
            "provider": row.provider,
            "version": row.version,
            "chromosomes_checked": known is not None,
        },
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )


def register_annotation(
    assembly: str,
    name: str,
    *,
    force: bool = False,
    cache_dir: str | Path | None = None,
    progressbar: bool = True,
    metadata: AnnotationMetadata | None = None,
    check_chromosomes: bool = True,
) -> dict[str, object]:
    """Register ``name`` for ``assembly`` and return the record of what that did.

    :func:`fetch_annotation` addressed by assembly name rather than by directory, and
    answering with the record rather than the paths — the one call the CLI makes, and
    the one a script makes when it wants to serialize what happened. An annotation that
    is already registered is returned from its record without fetching anything.

    Parameters
    ----------
    assembly : str
        The assembly the annotation belongs to, e.g. ``"hg38"``.
    name : str
        The **Registered name** the table lists, e.g. ``"gencode_v50"``.
    force : bool, default False
        Register again from scratch, repairing a directory that raises.
    cache_dir : str or pathlib.Path, optional
        Override which **assembly** directory the annotation is filed under. Defaults to
        :func:`assembly_data_dir(assembly) <genome.io.download.assembly_data_dir>`.
    progressbar : bool, default True
        Show a download progress bar (requires ``tqdm``).
    metadata : genome.metadata.AnnotationMetadata, optional
        A complete annotation record to use instead of the curated table's row.
    check_chromosomes : bool, default True
        Check the GTF's chromosome names against the assembly's. Pass ``False`` to
        register an annotation whose mismatch you have inspected and accept.

    Returns
    -------
    dict
        The completion record's own fields — ``files``, ``source_url``, ``sha256``,
        ``details``, ``completed_at`` and the rest — plus the ``assembly`` and the
        ``directory`` they live in. Ready to serialize as it is.

    Raises
    ------
    ValueError
        If the table lists no annotation ``name`` for ``assembly``.
    ChromosomeMismatchError
        If the GTF names sequences the assembly does not carry.
    genome.io.completion.RegistrationError
        If the directory holds a build that cannot be trusted as finished, or (with
        ``force``) if the run somehow left no record behind.
    genome.io.utils.ChecksumMismatchError
        If the row pins a sha256 and the unpacked GTF is not it.

    Examples
    --------
    >>> register_annotation("sacCer3", "ensgene_v101")   # doctest: +SKIP
    {'kind': 'annotation', 'name': 'ensgene_v101', 'files': {...}, ...}
    """
    assembly_dir = (
        Path(cache_dir).expanduser()
        if cache_dir is not None
        else download.assembly_data_dir(assembly)
    )
    fetch_annotation(
        assembly_dir,
        assembly,
        name,
        force=force,
        progressbar=progressbar,
        metadata=metadata,
        check_chromosomes=check_chromosomes,
    )
    directory = annotation_dir(assembly_dir, name)
    record = read_record(directory)
    if record is None:
        raise RegistrationError(
            f"{name} was registered for {assembly} in {directory} but no {RECORD_NAME} is "
            f"there, so nothing can vouch for it. Register it again with "
            f"`{_repair_command(assembly, name)}`."
        )
    payload: dict[str, object] = dict(asdict(record))
    payload["assembly"] = assembly
    payload["directory"] = str(directory)
    return payload


def _already_registered(directory: Path, *, force: bool, repair: str) -> bool:
    """Return whether ``directory`` holds a finished annotation that needs no work.

    The record is the only thing consulted, and the two ways a directory can contradict
    it raise from :func:`~genome.io.completion.check_registration`. ``force`` skips the
    question entirely, which is what makes it the repair.
    """
    if force:
        return False
    return check_registration(directory, repair=repair) is not None


def _proven_gtf(gtf: Path, expected: str | None) -> str | None:
    """Return the placed GTF's digest when it is provably the pinned one, else ``None``.

    What makes repairing cheap: a re-registration that can prove the GTF on disk is the
    pinned one keeps it and rebuilds only the database. ``None`` — fetch the source
    again — in all three of the cases where it cannot be proven: the GTF is missing, its
    digest is a different one, or **the row pins no digest at all**, since with nothing
    to compare against there is no way to show that what is there is right.
    """
    if expected is None or not gtf.is_file():
        return None
    actual = sha256_file(gtf)
    return actual if actual == expected else None


def _elide(names: Sequence[str], limit: int = _MAX_LISTED_NAMES) -> str:
    """Return ``names`` comma-joined, cut to ``limit`` and counting what was cut."""
    listed = ", ".join(names[:limit])
    hidden = len(names) - limit
    return listed if hidden <= 0 else f"{listed} (and {hidden} more)"


def _open_text(path: Path) -> IO[str]:
    """Open ``path`` for line-by-line reading, decompressing a ``.gz`` as it goes."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def _gtf_chromosomes(gtf: Path) -> set[str]:
    """Return the distinct sequence names ``gtf`` uses, reading it one line at a time.

    A GENCODE GTF is well over a gigabyte unpacked, so this is a single streaming pass
    that keeps only the distinct values of the first column — the file is never held in
    memory, and a ``.gz`` is decompressed as it streams rather than unpacked first.
    Comment lines and anything without a column separator are skipped, so a header never
    becomes a chromosome.
    """
    names: set[str] = set()
    with _open_text(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, separator, _ = line.partition("\t")
            if separator:
                names.add(chrom)
    return names


def _assembly_chromosomes(chrom_sizes: Path | None) -> frozenset[str] | None:
    """Return the names in ``chrom_sizes``, or ``None`` when there is none to read.

    ``None`` means *nothing to check against* rather than *no chromosomes*: an
    annotation can be registered before its assembly is prepared, and then no
    ``chrom.sizes`` exists yet. The registration then records that the names were not
    checked rather than claiming they were.
    """
    if chrom_sizes is None or not chrom_sizes.is_file():
        return None
    return frozenset(str(name) for name in read_chrom_sizes(chrom_sizes).index)


def _reject_unknown_chromosomes(gtf: Path, known: frozenset[str] | None, *, name: str) -> None:
    """Raise :class:`ChromosomeMismatchError` if ``gtf`` names anything outside ``known``.

    One-directional on purpose: only the names the GTF uses are looked up, since an
    assembly carrying scaffolds the annotation never mentions is normal. ``known`` of
    ``None`` is the check turned off, or no ``chrom.sizes`` to run it against, and
    nothing happens. Call it *before* the database build and before the GTF is placed —
    a mismatch then costs a streaming pass rather than the minutes a build takes, and
    leaves the annotation directory as it was found.
    """
    if known is None:
        return
    missing = _gtf_chromosomes(gtf) - known
    if missing:
        raise ChromosomeMismatchError(name, missing, sorted(known))


def _fetch_gtf(
    annotation: GtfAnnotation,
    row: AnnotationMetadata,
    *,
    progressbar: bool = True,
    known: frozenset[str] | None = None,
) -> str:
    """Download ``row``'s GTF, verify the unpacked file, place it, and return its digest.

    Both the archive and the file it unpacks to land in the annotation's working area,
    which is on the same filesystem, so placing the GTF is a rename rather than a copy
    and the archive survives an interrupted run. Both verifications — the pinned digest
    and the chromosome names — happen while the GTF is still in the working area, so a
    GTF this assembly cannot use never reaches the annotation directory.
    """
    directory = annotation.gtf.parent
    gzipped = row.url.endswith(".gz")
    # Named after the annotation, not after the URL: a provider's file name says
    # nothing about the name the lab registered it under.
    fetched = download.fetch_url(
        row.url,
        work_dir(directory),
        fname=f"{annotation.name}.gtf.gz" if gzipped else annotation.gtf.name,
        processor=pooch.Decompress(method="gzip", name=annotation.gtf.name) if gzipped else None,
        progressbar=progressbar,
    )
    digest = sha256_file(fetched)
    if row.sha256 is not None and digest != row.sha256:
        raise ChecksumMismatchError(fetched, row.sha256, digest)
    _reject_unknown_chromosomes(fetched, known, name=annotation.name)
    directory.mkdir(parents=True, exist_ok=True)
    fetched.replace(annotation.gtf)
    return digest


def _build_and_record(
    annotation: GtfAnnotation,
    *,
    source_url: str,
    sha256: str,
    details: dict[str, Any],
    disable_infer_genes: bool,
    disable_infer_transcripts: bool,
) -> GtfAnnotation:
    """Build the database beside the placed GTF, then write the record that ends the job.

    The record is written last, once both files exist, and the working area goes only
    after it — so an interrupted run leaves its download in place and repairs without
    fetching the GTF again.

    gffutils' version is recorded in ``details`` rather than in ``tool_versions``,
    which is for **External tool**s: a tool resolved on ``PATH``, version-detected by
    running it, and installable by a command an error can name. gffutils is an installed
    Python library and none of that applies to it, so recording it there would blur the
    one word the package keeps sharp for binaries it shells out to.
    """
    database = gffutils.create_db(
        str(annotation.gtf),
        str(annotation.db),
        # Reached only when the annotation is being built, so an older database left by
        # an interrupted or forced re-registration is replaced rather than refused.
        force=True,
        keep_order=True,
        merge_strategy="create_unique",
        sort_attribute_values=True,
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )
    # The on-disk database is now fully written; release the SQLite connection
    # so we don't leak an open file handle (the build is the only thing we need).
    database.conn.close()

    directory = annotation.gtf.parent
    record = build_record(
        directory,
        kind="annotation",
        name=annotation.name,
        files=[annotation.gtf, annotation.db],
        source_url=source_url,
        sha256=sha256,
        details={**details, "gffutils_version": gffutils.__version__},
    )
    write_record(directory, record)
    clear_work_dir(directory)
    return annotation
