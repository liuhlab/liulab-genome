"""Putting an **Annotation** on disk — the fetch, the check, the build, and the record.

I/O boundary module, and the only one of the four that writes. Everything that places an
annotation is here: the fetch, the placement under its **Registered name**, the check of
its **Chromosome** names against the assembly's ``chrom.sizes``, the repair-command
strings a broken directory names, and the write of the **Completion marker** that is the
only thing saying an annotation is registered. Where those files go is here too, since
placing one is the first thing that needs to know.

There are two ways to add an annotation, and
:class:`~genome.io.annotation.registry.AnnotationRegistry` is the way in to both. By
**name**: :meth:`~genome.io.annotation.registry.AnnotationRegistry.register` takes the
name the curated annotation table lists for this assembly, fetches that row's URL, checks
the unpacked GTF against the sha256 the row pins (ADR-0006), builds the database and
writes the record. By **path**:
:meth:`~genome.io.annotation.registry.AnnotationRegistry.register_path` is the escape
hatch for a GTF no row lists — the caller says where the file is, and it is placed, built
and recorded the same way. :func:`register_annotation` and :func:`register_gtf` are those
same two addressed by assembly name and answering with the record rather than the paths,
one apiece, matching ``genome register-annotation`` and ``genome register-gtf`` exactly;
both build a registry for the length of the call, so they add no second code path. They
register, which is why they are here rather than beside the other assembly-addressed
functions. What they answer *with* — :class:`RegisteredAnnotation` — is here too, beside
whatever returns it (ADR-0022). A registration carries the **Completion marker** whole
rather than copying it out field by field, so every question a surface then asks is
answered off the record in hand rather than by reading the directory again.

A third way in has exactly one caller. :func:`register_merged_gtf` writes the **Merged
annotation** a **Chimera** build derives from its components' own annotations, inside the
act that writes the chimera's FASTA (ADR-0008). Nothing is fetched, so its record pins no
source and no table row describes it; and nothing on disk is ever adopted, so the build
that owns it writes it every time it runs. Owning it cuts both ways: the merged name is
derived from what contributed, so a rebuild whose contributors changed writes a *different*
name, and :func:`discard_merged_annotation` takes the one that build no longer owns. Both
take a bare **Assembly dir** and adopt nothing, so neither has ever needed the registry's
four-way state, and :mod:`genome.io.chimera` is the only caller of either.

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
**Completion marker** exists to distrust. So re-registering something that already has a
record returns it silently, and a directory holding files but no record raises and names
its repair (ADR-0007). Which command that is depends on which route registered it, and the
three spellings of it live together here so a renamed command is renamed once.

Examples
--------
>>> from pathlib import Path
>>> from genome.io.annotation.registration import annotation_dir
>>> annotation_dir(Path("/data/genome/sacCer3"), "ensgene_v101").name
'ensgene_v101'
"""

from __future__ import annotations

import gzip
import hashlib
import shlex
import shutil
from collections.abc import Container, Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any

import pooch

from genome.chimera import suffixed
from genome.io import fetch
from genome.io.annotation.database import LIBRARY_VERSION_KEY, build_database
from genome.io.completion import (
    RECORD_NAME,
    CompletionRecord,
    RegistrationError,
    build_record,
    check_registration,
    clear_work_dir,
    read_record,
    work_dir,
    write_record,
)
from genome.io.fasta import read_chrom_sizes
from genome.io.registration import ANNOTATIONS_SUBDIR, assembly_repair_command
from genome.io.utils import ChecksumMismatchError, _gunzip, sha256_file
from genome.metadata import AnnotationMetadata

if TYPE_CHECKING:
    from genome.io.annotation.registry import AnnotationRegistry


#: Subdirectory under an assembly's data dir holding all its GTF annotations.
#: The Assembly context owns the layout, so the name is read from there.
_GTF_SUBDIR = ANNOTATIONS_SUBDIR


#: How many names an error lists before saying how many it left out. A whole-genome
#: mismatch offends in the thousands, and a message that long is one nobody reads.
_MAX_LISTED_NAMES = 10


#: What a repair command puts where a GTF's path belongs when nothing on disk remembers
#: it. Deliberately not a path: a command naming a file that is not there is one that
#: fails when it is pasted, which is worse than one visibly asking to be filled in.
_UNKNOWN_PATH = "<path>"


#: ``details`` key marking a **Merged annotation**: one entry per contributing component,
#: naming the component and the annotation of its own that went in. It is what tells a
#: reader — and :func:`_annotation_repair` — that this annotation was derived here rather
#: than fetched or handed in by path, so neither of those commands would repair it.
_MERGED_FROM_KEY = "merged_from"


#: Keys of one entry under :data:`_MERGED_FROM_KEY`.
_MERGED_COMPONENT_KEY = "component"


_MERGED_ANNOTATION_KEY = "annotation"


#: What ``details["chromosomes_unchecked_because"]`` says when the caller stood the check
#: down — ``check_chromosomes=False``, or ``--no-check-chromosomes`` from a shell. There is
#: no advice to give about it: the assembly may be registered and the names deliberately
#: accepted, so all a surface can say is what the record therefore does not vouch for.
UNCHECKED_CALLER_OVERRIDE = "caller-override"


#: …and when the check was asked for but had nothing to run against, the assembly having
#: no ``chrom.sizes`` yet. Registering the assembly is what makes the check possible, so
#: this is the one of the two states where saying so is useful advice rather than noise.
UNCHECKED_NO_CHROM_SIZES = "no-chrom-sizes"


#: What each state of the chromosome check reads as, one sentence apiece — including the
#: one where it ran and passed, since a surface that says nothing about it reads as a pass.
#: Keyed by ``details["chromosomes_unchecked_because"]``; ``None`` is the check that ran.
_CHECK_SUMMARIES = {
    None: "chromosomes checked — every name the GTF uses is one the assembly carries",
    UNCHECKED_NO_CHROM_SIZES: (
        "chromosomes not checked — nothing to check against; register the assembly first "
        "to verify them"
    ),
    UNCHECKED_CALLER_OVERRIDE: (
        "chromosomes not checked — the check was stood down, so the record does not vouch "
        "for the names"
    ),
}


#: What a record written before the reason was recorded reads as. Its bare ``False`` was
#: written for either reason and nothing on disk says which, so it is reported as neither.
_UNKNOWN_REASON_SUMMARY = (
    "chromosomes not checked — this record does not say why, so whether the names match "
    "the assembly is unknown"
)


def annotation_register_command(assembly: str, name: str) -> str:
    """Return the command that registers the annotation ``name`` for ``assembly``.

    The one spelling of it. Errors quote it, the repair adds ``--force`` to it, and
    :attr:`~genome.io.annotation.registry.AnnotationStatus.default_summary` names it for a **Default annotation** nobody
    has fetched yet — so a renamed command is renamed once.

    Parameters
    ----------
    assembly : str
        The **Assembly** the annotation belongs to, e.g. ``"hg38"``.
    name : str
        The **Registered name** to address it by.

    Returns
    -------
    str
        A shell command, unquoted and unfenced — the caller decides how to set it.

    Examples
    --------
    >>> annotation_register_command("hg38", "gencode_v50")
    'genome register-annotation hg38 gencode_v50'
    """
    return f"genome register-annotation {assembly} {name}"


def chromosome_check_summary(details: Mapping[str, Any]) -> str:
    """Return the one line a surface prints about an annotation's chromosome-name check.

    Four states, four sentences, and one of them is always returned: the check ran and
    passed; it had nothing to run against, and registering the assembly is what fixes
    that; the caller stood it down, which is not something to advise about; or the record
    does not say which, and none of the three may be claimed. Silence is not a fifth
    state — a surface that prints nothing about the check reads as one that passed.

    ``details`` is a registration record's ``details``; a caller holding what a
    registration answered with asks :attr:`RegisteredAnnotation.chromosome_check` instead
    and never spells the two fields. Those are ``chromosomes_checked`` — the check ran and
    the GTF's names were all the assembly's — and ``chromosomes_unchecked_because``, which
    says which of the two reasons it did not, and is ``None`` when it did.

    A record written before the second field existed carries a bare
    ``chromosomes_checked: false`` that was written for either reason, and nothing on disk
    says which. It reads as *unknown* rather than as either one, and rather than raising:
    the reason is a fact that was never gathered, which is what an absent entry in
    ``tool_versions`` means too.

    Parameters
    ----------
    details : mapping of str to object
        A registration record's ``details``. Anything else it holds is ignored, and a
        mapping holding neither field reads as unknown.

    Returns
    -------
    str
        One sentence, with no trailing punctuation and no leading indent — the caller
        decides how to set it.

    Examples
    --------
    >>> chromosome_check_summary({"chromosomes_checked": True})
    'chromosomes checked — every name the GTF uses is one the assembly carries'
    >>> print(chromosome_check_summary({"chromosomes_unchecked_because": "caller-override"}))
    chromosomes not checked — the check was stood down, so the record does not vouch for the names
    """
    if details.get("chromosomes_checked") is True:
        return _CHECK_SUMMARIES[None]
    # Anything else — the field absent, or a reason a later version writes and this one
    # has never heard of — is a reason that cannot be reported, which is the unknown.
    because = details.get("chromosomes_unchecked_because")
    if isinstance(because, str) and because in _CHECK_SUMMARIES:
        return _CHECK_SUMMARIES[because]
    return _UNKNOWN_REASON_SUMMARY


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
    genome.io.annotation.registration.ChromosomeMismatchError: the GTF for 'gencode_v44' ...
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
            f"check_chromosomes=False — --no-check-chromosomes from a shell — to register "
            f"this one anyway."
        )


@dataclass(frozen=True)
class GtfAnnotation:
    """A registered GTF annotation: its name and the on-disk GTF + database paths."""

    name: str
    gtf: Path
    db: Path


@dataclass(frozen=True)
class MergeSource:
    """One component's contribution to a **Merged annotation**.

    What :func:`register_merged_gtf` needs about a single component: whose sequences the
    features sit on, which of that component's annotations was taken, and where its GTF
    is. The component name is not decoration — it is the suffix every seqname the merge
    writes carries (ADR-0009), so the merged features land on the chimera's own
    chromosome names.

    Attributes
    ----------
    component : str
        The **Component** assembly name, alphanumeric.
    annotation : str
        The **Registered name** of that component's contributing annotation.
    gtf : pathlib.Path
        That annotation's placed GTF, read one line at a time and never in full.

    Examples
    --------
    >>> from pathlib import Path
    >>> MergeSource("ce11", "wormbase_ws298", Path("/data/ce11.gtf")).component
    'ce11'
    """

    component: str
    annotation: str
    gtf: Path


@dataclass(frozen=True)
class RegisteredAnnotation:
    """What registering one annotation produced: its record, and where that landed.

    :func:`register_annotation`'s answer and :func:`register_gtf`'s — what ``genome
    register-annotation`` and ``genome register-gtf`` print, and what their ``--json``
    serializes. A :class:`GtfAnnotation` says where an annotation's two files are; this
    says what the run that wrote them did, which is the **Completion marker** itself,
    carried whole. Every question a surface then asks — the digest, the source, the files
    claimed, whether the chromosome names were actually checked — is answered from that one
    record rather than by reading the directory again.

    Attributes
    ----------
    assembly : str
        The **Assembly** the annotation belongs to. It is not in the record, which names
        the annotation rather than what it annotates.
    directory : pathlib.Path
        The annotation's own directory, ``<assembly dir>/gtf/<name>/``.
    record : genome.io.completion.CompletionRecord
        The record the registration wrote, read back.

    Examples
    --------
    >>> from pathlib import Path
    >>> from genome.io.completion import CompletionRecord
    >>> registered = RegisteredAnnotation(
    ...     assembly="hg38",
    ...     directory=Path("/data/genome/hg38/gtf/gencode_v50"),
    ...     record=CompletionRecord(
    ...         kind="annotation",
    ...         name="gencode_v50",
    ...         files={"gencode_v50.gtf": 12, "gencode_v50.db": 34},
    ...         source_url="https://example.org/gencode_v50.gtf.gz",
    ...         sha256="1a2b3c",
    ...         tool_versions={},
    ...         package_version="2026.8.0",
    ...         completed_at="2026-08-12T09:00:00+00:00",
    ...         details={"chromosomes_checked": True},
    ...     ),
    ... )
    >>> registered.name, registered.file_names
    ('gencode_v50', ['gencode_v50.db', 'gencode_v50.gtf'])
    >>> print(registered.chromosome_check)
    chromosomes checked — every name the GTF uses is one the assembly carries
    """

    assembly: str
    directory: Path
    record: CompletionRecord

    @property
    def name(self) -> str:
        """The **Registered name** it is addressed by — the record's own name."""
        return self.record.name

    @property
    def source_url(self) -> str | None:
        """The URL fetched, or the path a GTF was handed over at; ``None`` for a merge."""
        return self.record.source_url

    @property
    def sha256(self) -> str | None:
        """Digest of the placed GTF, or ``None`` when none was computed."""
        return self.record.sha256

    @property
    def file_names(self) -> list[str]:
        """Every file the record claims, sorted — a fresh list each call."""
        return sorted(self.record.files)

    @property
    def chromosome_check(self) -> str:
        """The one line saying what the chromosome-name check settled for this annotation.

        :func:`chromosome_check_summary` over the record this registration wrote, so the
        surface that prints it never reads the record's own keys. Always a sentence:
        silence would read as a pass.
        """
        return chromosome_check_summary(self.record.details)

    def as_json(self) -> dict[str, Any]:
        """Return this registration as ``--json`` serializes it.

        The record's own fields under the record's own names, then the ``assembly`` it
        belongs to and the ``directory`` it landed in — the two facts a record does not
        hold about itself. The names are the ones written on disk and are never respelled
        here.

        Returns
        -------
        dict
            The record's fields, followed by ``assembly`` and ``directory``.
        """
        return {**asdict(self.record), "assembly": self.assembly, "directory": str(self.directory)}


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


def _register_gtf_command(assembly: str, source: str, name: str) -> str:
    """Return the command that registers the GTF at ``source`` as ``name``.

    ``source`` is rendered by the caller, so a message about a GTF nobody has yet named
    can say ``<path>`` where one about a real file says the file, shell-quoted.
    """
    return f"genome register-gtf {assembly} {source} {name}"


def _repair_command(assembly: str, name: str) -> str:
    """Return the command that registers ``name`` again from scratch.

    Quoted verbatim into every error a broken annotation directory raises, so it has to
    be a command that exists and does the job.
    """
    return f"{annotation_register_command(assembly, name)} --force"


def _path_repair_command(assembly: str, source: str, name: str) -> str:
    """Return the command that registers the GTF at ``source`` again from scratch.

    ``source`` is rendered by the caller, as :func:`_register_gtf_command` takes it: a
    file that is there is shell-quoted, and one nothing remembers the path of is
    :data:`_UNKNOWN_PATH`.
    """
    return f"{_register_gtf_command(assembly, source, name)} --force"


def _annotation_repair(directory: Path, *, assembly: str, offered: Container[str]) -> str:
    """Return the command that registers ``directory``'s annotation again from scratch.

    Which route repairs a broken annotation is decided by which route registered it, and
    the three differ. A **Merged annotation** was written by a chimera build and by
    nothing else, so neither registering it by name nor handing it a GTF would rebuild it:
    what repairs it is rebuilding the chimera, and its record is asked first because that
    is the only place the fact is written down. A listed one is fetched again from the row
    that lists it. An unlisted one has to be handed the GTF it was built from — and the
    record is what remembers that path, so an annotation whose record is gone, or whose
    source has since moved, can only name the command with the path left to fill in. That
    is the honest answer, and the alternative is printing a path that is not there.
    """
    name = directory.name
    record = read_record(directory)
    if record is not None and record.details.get(_MERGED_FROM_KEY):
        return assembly_repair_command(assembly)
    if name in offered:
        return _repair_command(assembly, name)
    source = Path(record.source_url) if record is not None and record.source_url else None
    if source is not None and source.is_file():
        return _path_repair_command(assembly, shlex.quote(str(source)), name)
    return _path_repair_command(assembly, _UNKNOWN_PATH, name)


def _register_gtf(
    assembly_dir: Path,
    source: Path,
    name: str,
    *,
    repair: str,
    force: bool,
    chrom_sizes: str | Path | None,
    check_chromosomes: bool,
    disable_infer_genes: bool,
    disable_infer_transcripts: bool,
) -> GtfAnnotation:
    """Place, build and record the GTF at ``source``, as ``register_path`` describes.

    Addressed by directory, because placing a file needs one and nothing here needs the
    assembly's name — which is why ``repair``, the command a broken annotation directory
    tells the caller to run, is handed in by the caller that does know it.
    """
    if not source.is_file():
        raise FileNotFoundError(
            f"GTF file not found: {source}. Pass the path of an existing .gtf or "
            f".gtf.gz, or register a listed annotation by name instead."
        )

    annotation = _annotation_files(assembly_dir, name)
    directory = annotation_dir(assembly_dir, name)
    if _already_registered(directory, force=force, repair=repair):
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
        details=_chromosome_check_details(known, requested=check_chromosomes),
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )


def register_merged_gtf(
    assembly_dir: Path,
    name: str,
    sources: Sequence[MergeSource],
    *,
    separator: str,
    chrom_sizes: str | Path,
    disable_infer_genes: bool = True,
    disable_infer_transcripts: bool = True,
) -> GtfAnnotation:
    """Write the **Merged annotation** of ``sources`` under ``name`` and build its database.

    The annotation half of a **Chimera** build, called from
    :meth:`~genome.io.chimera.ChimeraBuilder.build_genome` and from nowhere else. Each
    source's GTF is streamed a line at a time into one file whose seqnames carry the
    component suffix the chimera's FASTA already carries (ADR-0009), which is then placed,
    checked, built and recorded exactly as any other annotation is.

    **No coordinate is converted.** Only the first column of each data line is rewritten;
    every byte after the first tab — including both position fields, which are 1-based and
    inclusive as GTF has them — is copied through untouched. The features are the
    components' own features on the components' own sequences, under a new spelling of the
    sequence name and nothing else.

    Comment lines are **dropped**, all of them. A ``#!genome-build`` pragma names the
    single assembly its file was built for, and several of those concatenated would each
    be false about the chimera; the ordinary ``#`` comment beside them describes a file
    that no longer exists as such. Nothing else is dropped: a line carrying a tab is a data
    line and survives, unsorted and in component order.

    The chromosome-name check is **not optional here** and has no argument that stands it
    down. Everything else in the build derives the chimera's names twice — once for the
    FASTA and once for this — and the check is the one place those two answers are set
    against each other, so a merge that misspelled a name raises
    :class:`ChromosomeMismatchError` rather than registering an annotation that queries
    empty.

    Nothing on disk is adopted: unlike the other ways in, this one never asks whether the
    annotation is already registered. It is written by the build that owns it, every time
    that build runs, which is what makes a stale database impossible to hand back — the
    name is derived from the contributing annotations, so it changes when they do, but it
    cannot say *which* components contributed and would otherwise be reusable under a
    meaning it no longer has.

    Parameters
    ----------
    assembly_dir : pathlib.Path
        The chimera's **Assembly dir**, which this annotation is filed under.
    name : str
        The **Registered name** to write it as — derived by the caller from the
        contributing annotations' names.
    sources : sequence of MergeSource
        One entry per contributing component, in the order their sequences are written.
        Must not be empty: no contributors means no annotation, which the caller decides
        rather than registering an empty one.
    separator : str
        The run of underscores this chimera's chromosome names carry, as
        :func:`~genome.chimera.derive_separator` gave it.
    chrom_sizes : str or pathlib.Path
        The chimera's ``chrom.sizes``, whose names every merged seqname must be among.
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
        If ``sources`` is empty.
    ChromosomeMismatchError
        If a merged seqname is not one the chimera carries — the merge and the FASTA
        build disagreeing, which nothing else would catch.
    genome.chimera.ChimeraNamingError
        If a component name or the separator does not obey the naming contract.

    Examples
    --------
    >>> from pathlib import Path
    >>> register_merged_gtf(                             # doctest: +SKIP
    ...     Path("/data/genome/ce11_ecHT115"),
    ...     "wormbase_ws298+refseq_rs_2025_06_26",
    ...     [MergeSource("ce11", "wormbase_ws298", Path("/data/ce11.gtf"))],
    ...     separator="__",
    ...     chrom_sizes=Path("/data/genome/ce11_ecHT115/ce11_ecHT115.chrom.sizes"),
    ... )
    GtfAnnotation(name='wormbase_ws298+refseq_rs_2025_06_26', ...)
    """
    if not sources:
        raise ValueError(
            f"a merged annotation for {assembly_dir.name!r} needs at least one contributing "
            f"component, got none. A chimera whose components carry no annotation registers "
            f"none at all rather than an empty one — do not call this with an empty list."
        )
    annotation = _annotation_files(assembly_dir, name)
    known = _assembly_chromosomes(Path(chrom_sizes))
    # Written into the working area and checked there, so a merge the chimera cannot use
    # never reaches the annotation directory. Same filesystem, so placing it is a rename.
    staged = work_dir(annotation_dir(assembly_dir, name)) / annotation.gtf.name
    staged.parent.mkdir(parents=True, exist_ok=True)
    digest = _write_merged_gtf(sources, staged, separator=separator)
    _reject_unknown_chromosomes(staged, known, name=name)
    staged.replace(annotation.gtf)

    return _build_and_record(
        annotation,
        source_url=None,
        sha256=digest,
        details={
            _MERGED_FROM_KEY: [
                {
                    _MERGED_COMPONENT_KEY: source.component,
                    _MERGED_ANNOTATION_KEY: source.annotation,
                }
                for source in sources
            ],
            **_chromosome_check_details(known, requested=True),
        },
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )


def _write_merged_gtf(sources: Sequence[MergeSource], destination: Path, *, separator: str) -> str:
    """Write ``sources`` into one GTF at ``destination`` and return the sha256 of it.

    One streaming pass per source, a line at a time, hashing what is written as it is
    written — so a GENCODE-sized annotation is neither held in memory nor read again to
    produce the digest the record carries. A data line's seqname is **extended** rather
    than rebuilt: the bytes before the first tab get the suffix appended and every byte
    after it is copied verbatim, which is why no coordinate is touched. Comment lines and
    anything carrying no tab — a blank line, a stray fragment — are dropped, since neither
    names a sequence and a line the chromosome check never saw must not reach the file.

    The suffix is spelled once per source by :func:`~genome.chimera.suffixed`, the same
    function the FASTA build spells its headers with, rather than assembled here: that is
    what stops the two halves of one chimera drifting apart about a name.
    """
    digest = hashlib.sha256()
    with destination.open("wb") as output:
        for source in sources:
            # suffixed("", …) is exactly the tail every name of this component gains, and
            # it validates the component and the separator once instead of per line.
            suffix = suffixed("", source.component, separator).encode()
            with source.gtf.open("rb") as handle:
                for line in handle:
                    if line.startswith(b"#"):
                        continue
                    chromosome, tab, rest = line.partition(b"\t")
                    if not tab:
                        continue
                    written = chromosome + suffix + tab + rest
                    if not written.endswith(b"\n"):
                        written += b"\n"
                    output.write(written)
                    digest.update(written)
    return digest.hexdigest()


def discard_merged_annotation(assembly_dir: Path, name: str) -> bool:
    """Remove the **Merged annotation** registered as ``name``, when that is what it is.

    The other half of a chimera build owning its annotation. The merged name is the
    ``+``-join of the contributing annotations' names, so a rebuild whose contributing set
    changed registers the merge under a *new* name — and the previous one, which nothing
    else will ever write again, would otherwise stay registered beside it. Two derived
    annotations with nothing to choose between them is a chimera whose **Default
    annotation** is suddenly none, which is how an annotated chimera comes back from a
    legitimate repair with none at all. So the build removes what it no longer owns, and
    :meth:`~genome.io.chimera.ChimeraBuilder.build_genome` is the only caller.

    Owning it is **proved, not assumed**: only a directory whose record carries the
    ``merged_from`` marker a merge writes is removed, so an annotation a caller registered
    by hand — and a directory nothing vouches for — is left exactly where it is, whatever
    it is called.

    Parameters
    ----------
    assembly_dir : pathlib.Path
        The chimera's **Assembly dir**, which the annotation is filed under.
    name : str
        The **Registered name** to remove, as the previous build's completion record
        names it.

    Returns
    -------
    bool
        Whether an annotation was removed. ``False`` for a name nothing is registered
        under, and for one whose record does not show a merge wrote it.

    Examples
    --------
    >>> from pathlib import Path
    >>> discard_merged_annotation(Path("/tmp/definitely-not-an-assembly"), "a+b")
    False
    """
    directory = annotation_dir(assembly_dir, name)
    record = read_record(directory)
    if record is None or not record.details.get(_MERGED_FROM_KEY):
        return False
    shutil.rmtree(directory)
    # A chimera that merged nothing carries no `gtf/` tree at all, and one whose last
    # derived annotation has just gone is in exactly that state.
    root = _annotations_root(assembly_dir)
    if not any(root.iterdir()):
        root.rmdir()
    return True


def _registry_for(assembly: str, cache_dir: str | Path | None) -> AnnotationRegistry:
    """Return the registry the two assembly-addressed registrars build for one call.

    The one place this module reaches back into the registry, and it reaches at call time
    rather than at import time. The registry is built out of what placement writes — the
    scans, the repair commands, the four states — so it imports this module at the top of
    its file; holding the class here would make that a cycle. Both registrars go through
    here, so the crossing is spelled once and stays visible.
    """
    from genome.io.annotation.registry import AnnotationRegistry

    return AnnotationRegistry.locate(assembly, cache_dir)


def register_annotation(
    assembly: str,
    name: str,
    *,
    force: bool = False,
    cache_dir: str | Path | None = None,
    progressbar: bool = True,
    metadata: AnnotationMetadata | None = None,
    check_chromosomes: bool = True,
    disable_infer_genes: bool = True,
    disable_infer_transcripts: bool = True,
) -> RegisteredAnnotation:
    """Register ``name`` for ``assembly`` and return the record of what that did.

    :meth:`~genome.io.annotation.registry.AnnotationRegistry.register` addressed by assembly
    name, and answering with the
    record rather than the paths — the call ``genome register-annotation`` makes, and the
    one a script makes when it wants to serialize what happened. An annotation that is
    already registered is returned from its record without fetching anything.
    :func:`register_gtf` is the same shape for a GTF the table does not
    list.

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
        register an annotation whose mismatch you have inspected and accept; the record
        then says the check was stood down, rather than merely that it did not run.
    disable_infer_genes : bool, default True
        Do not reconstruct ``gene`` features from exon lines.
    disable_infer_transcripts : bool, default True
        Do not reconstruct ``transcript`` features from exon lines.

    Returns
    -------
    RegisteredAnnotation
        The completion record the run wrote — ``files``, ``source_url``, ``sha256``,
        ``details``, ``completed_at`` and the rest — with the ``assembly`` it belongs to
        and the ``directory`` it lives in. :meth:`RegisteredAnnotation.as_json`
        serializes it.

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
    RegisteredAnnotation(assembly='sacCer3', directory=PosixPath('...'), record=...)
    """
    annotation = _registry_for(assembly, cache_dir).register(
        name,
        force=force,
        progressbar=progressbar,
        metadata=metadata,
        check_chromosomes=check_chromosomes,
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )
    return _registered_annotation(
        annotation, assembly=assembly, repair=_repair_command(assembly, name)
    )


def register_gtf(
    assembly: str,
    gtf: str | Path,
    name: str,
    *,
    force: bool = False,
    cache_dir: str | Path | None = None,
    check_chromosomes: bool = True,
    disable_infer_genes: bool = True,
    disable_infer_transcripts: bool = True,
) -> RegisteredAnnotation:
    """Register the GTF at ``gtf`` for ``assembly`` and return the record of what that did.

    :meth:`~genome.io.annotation.registry.AnnotationRegistry.register_path` addressed by
    assembly name, and answering
    with the record rather than the paths — the call ``genome register-gtf`` makes, and
    the way a script registers an annotation the curated table does not list and then
    serializes what happened. :func:`register_annotation` is the same shape for one the
    table does list.

    Naming the assembly is what lets its ``chrom.sizes`` be found rather than passed, so
    an unlisted GTF has its chromosome names checked by default, exactly as a listed one
    does — and it is what says which reference these gene models are for. An assembly that
    is not prepared yet has no ``chrom.sizes`` to check against, and the record then says
    the names went unchecked — and that it was for want of that file, not because anyone
    stood the check down — rather than claiming they passed.

    Parameters
    ----------
    assembly : str
        The assembly the annotation belongs to, e.g. ``"hg38"``. Never inferred from the
        GTF: it says which reference these gene models are for (ADR-0003).
    gtf : str or pathlib.Path
        Path to the source GTF, plain or ``.gz``.
    name : str
        The **Registered name** to address it by, unique within the assembly.
    force : bool, default False
        Register again from scratch, repairing a directory that raises.
    cache_dir : str or pathlib.Path, optional
        Override which **assembly** directory the annotation is filed under. Defaults to
        :func:`assembly_data_dir(assembly) <genome.io.download.assembly_data_dir>`.
    check_chromosomes : bool, default True
        Check the GTF's chromosome names against the assembly's. Pass ``False`` to
        register a GTF whose mismatch you have inspected and accept; the record
        then says the check was stood down, rather than merely that it did not run.
    disable_infer_genes : bool, default True
        Do not reconstruct ``gene`` features from exon lines.
    disable_infer_transcripts : bool, default True
        Do not reconstruct ``transcript`` features from exon lines.

    Returns
    -------
    RegisteredAnnotation
        The completion record the run wrote, with the ``assembly`` and the ``directory``
        it lives in, exactly as :func:`register_annotation` returns them. The
        ``source_url`` is the path the GTF was taken from.

    Raises
    ------
    FileNotFoundError
        If ``gtf`` is not a file.
    ChromosomeMismatchError
        If the GTF names sequences the assembly does not carry.
    genome.io.completion.RegistrationError
        If the directory holds a build that cannot be trusted as finished, or (with
        ``force``) if the run somehow left no record behind.

    Examples
    --------
    >>> register_gtf(                     # doctest: +SKIP
    ...     "sacCer3", "custom.gtf.gz", "custom"
    ... )
    RegisteredAnnotation(assembly='sacCer3', directory=PosixPath('...'), record=...)
    """
    source = Path(gtf)
    annotation = _registry_for(assembly, cache_dir).register_path(
        source,
        name,
        force=force,
        check_chromosomes=check_chromosomes,
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )
    return _registered_annotation(
        annotation,
        assembly=assembly,
        repair=_path_repair_command(assembly, shlex.quote(str(source)), name),
    )


def _chromosome_check_details(known: frozenset[str] | None, *, requested: bool) -> dict[str, Any]:
    """Return what the record says about the chromosome check: whether it ran, and why not.

    ``known`` of ``None`` is a check that did not run, and the two reasons it can be are
    exactly what a record carrying only the bool could not tell apart. A check that ran
    and did not raise passed, so nothing is recorded beside the ``True``.
    """
    if known is not None:
        return {"chromosomes_checked": True, "chromosomes_unchecked_because": None}
    return {
        "chromosomes_checked": False,
        "chromosomes_unchecked_because": (
            UNCHECKED_NO_CHROM_SIZES if requested else UNCHECKED_CALLER_OVERRIDE
        ),
    }


def _registered_annotation(
    annotation: GtfAnnotation, *, assembly: str, repair: str
) -> RegisteredAnnotation:
    """Return the :class:`RegisteredAnnotation` a just-finished run left.

    The record is what registering *produced*, so a registration that reports success
    and leaves none is a contradiction rather than a missing file, and raises naming the
    command that does the job again. The annotation carries the directory to look in, so
    nothing here re-derives where it landed.
    """
    directory = annotation.gtf.parent
    record = read_record(directory)
    if record is None:
        raise RegistrationError(
            f"{annotation.name} was registered for {assembly} in {directory} but no "
            f"{RECORD_NAME} is there, so nothing can vouch for it. Register it again "
            f"with `{repair}`."
        )
    return RegisteredAnnotation(assembly=assembly, directory=directory, record=record)


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
    fetched = fetch.fetch_url(
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
    source_url: str | None,
    sha256: str,
    details: dict[str, Any],
    disable_infer_genes: bool,
    disable_infer_transcripts: bool,
) -> GtfAnnotation:
    """Build the database beside the placed GTF, then write the record that ends the job.

    The record is written last, once both files exist, and the working area goes only
    after it — so an interrupted run leaves its download in place and repairs without
    fetching the GTF again.

    The library that built the database says its own version, and it is recorded in
    ``details`` rather than in ``tool_versions`` for the reason
    :data:`~genome.io.annotation.database.LIBRARY_VERSION_KEY` gives.
    """
    library_version = build_database(
        annotation.gtf,
        annotation.db,
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )

    directory = annotation.gtf.parent
    record = build_record(
        directory,
        kind="annotation",
        name=annotation.name,
        files=[annotation.gtf, annotation.db],
        source_url=source_url,
        sha256=sha256,
        details={**details, LIBRARY_VERSION_KEY: library_version},
    )
    write_record(directory, record)
    clear_work_dir(directory)
    return annotation
