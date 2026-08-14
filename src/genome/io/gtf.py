"""Annotation registration — place a GTF, build its database, and record that it finished.

I/O boundary module. A reference assembly may carry several gene annotations (GENCODE,
RefSeq, WormBase, …). Each is registered under a **Registered name** and lives in its own
directory beside the assembly's sequence files::

    <LIULAB_DATA>/genome/<assembly>/gtf/<name>/
        <name>.gtf          # the annotation, kept decompressed
        <name>.db           # the gffutils SQLite database built from it
        .completion.json    # the record saying all of that finished
        .work/              # the disposable working area a fetch downloads into

:class:`AnnotationRegistry` is the way in. Bound once to one assembly — its name, its
**Assembly dir** and its ``chrom.sizes`` — it holds every annotation that assembly has and
answers everything about them: which are registered, which are broken, which the table
offers, which is the **Default annotation**, where one's GTF is, and the two acts that add
one. Everything that needs the four-way state asks a registry rather than assembling it
again: a :class:`~genome.genome.Genome` holds one for its lifetime, and each
assembly-addressed function here builds one for the length of the call.

There are two ways to add an annotation, each in two forms. By **name**:
:meth:`AnnotationRegistry.register` takes the name the curated annotation table lists for
this assembly, fetches that row's URL, checks the unpacked GTF against the sha256 the row
pins (ADR-0006), builds the database and writes the record; :func:`fetch_annotation` is it
addressed by directory. By **path**: :meth:`AnnotationRegistry.register_path` is the escape
hatch for a GTF no row lists — the caller says where the file is, and it is placed, built
and recorded the same way; :func:`register_gtf` is it addressed by directory, which is the
one form that knows no assembly name and so names a Python call where the rest name a
command. :func:`register_annotation` and :func:`register_annotation_by_path` are the two
that answer with the record rather than the paths, and those are what the CLI drives.

A third way in has exactly one caller. :func:`register_merged_gtf` writes the **Merged
annotation** a **Chimera** build derives from its components' own annotations, inside the
act that writes the chimera's FASTA (ADR-0008). Nothing is fetched, so its record pins no
source and no table row describes it; and nothing on disk is ever adopted, so the build
that owns it writes it every time it runs. Owning it cuts both ways: the merged name is
derived from what contributed, so a rebuild whose contributors changed writes a *different*
name, and :func:`discard_merged_annotation` takes the one that build no longer owns.

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

**Every annotation directory is registered, broken, or not begun**, and the middle one
has its own listing: :func:`list_broken_annotations`. Registering an annotation raises
over a directory it cannot trust, but listing must not — one annotation nobody can vouch
for cannot be allowed to stop a **Genome** opening or hide the annotations beside it — so
the two lists are reported side by side and each broken one carries the command that
repairs it.

**What the lab offers and what this machine holds are different questions.** The first is
the annotation table's to answer (:func:`~genome.metadata.list_annotation_metadata`), the
second this disk's (:func:`list_annotations`); :meth:`AnnotationRegistry.status` sets one
against the other, and :func:`default_annotation` is the one rule that picks the **Default
annotation** out of both. Those three scans plus that rule are what a registry is: they are
read together, once, and every later question is answered from the answer.

Examples
--------
>>> from pathlib import Path
>>> from genome.io.gtf import annotation_dir
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
from typing import IO, Any

import gffutils
import pooch

from genome.chimera import suffixed
from genome.io import fetch
from genome.io.completion import (
    RECORD_NAME,
    CompletionRecord,
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
from genome.io.registration import ANNOTATIONS_SUBDIR, AssemblyDir, assembly_repair_command
from genome.io.utils import ChecksumMismatchError, _gunzip, sha256_file
from genome.metadata import AnnotationMetadata, list_annotation_metadata, lookup_annotation

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
_CALLER_OVERRIDE = "caller-override"

#: …and when the check was asked for but had nothing to run against, the assembly having
#: no ``chrom.sizes`` yet. Registering the assembly is what makes the check possible, so
#: this is the one of the two states where saying so is useful advice rather than noise.
_NO_CHROM_SIZES = "no-chrom-sizes"

#: What each state of the chromosome check reads as, one sentence apiece — including the
#: one where it ran and passed, since a surface that says nothing about it reads as a pass.
#: Keyed by ``details["chromosomes_unchecked_because"]``; ``None`` is the check that ran.
_CHECK_SUMMARIES = {
    None: "chromosomes checked — every name the GTF uses is one the assembly carries",
    _NO_CHROM_SIZES: (
        "chromosomes not checked — nothing to check against; register the assembly first "
        "to verify them"
    ),
    _CALLER_OVERRIDE: (
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
            f"check_chromosomes=False — --no-check-chromosomes from a shell — to register "
            f"this one anyway."
        )


class AnnotationNotRegisteredError(KeyError):
    """No annotation of that name is registered here, so there is no path to hand back.

    Routinely *not* a mistake. An assembly's **Default annotation** comes from the
    curated table, and on a fresh machine the table's choice is exactly what nobody has
    registered yet — so a :class:`~genome.genome.Genome` opens with that default named
    and only asking for its path raises, naming the command that closes the gap. The
    other way in is a name nothing knows, and the message then says what the table does
    offer and how to register a GTF it does not list.

    A third way in is a directory that is there and cannot be trusted, which is not
    registered either. The next action is then neither of the above — registering it
    plainly would itself raise and demand ``--force`` — so ``broken`` carries what
    :func:`list_broken_annotations` found and the message quotes its repair, which is a
    command that runs as it stands.

    A :class:`KeyError`, because that is what asking a registry for a name it does not
    hold has always been.

    Parameters
    ----------
    assembly : str
        The assembly the annotation was asked for.
    name : str
        The **Registered name** that is not registered.
    registered : iterable of str
        The names that are registered on this machine.
    offered : iterable of str
        The names the annotation table offers for this assembly.
    broken : BrokenAnnotation, optional
        The broken registration filed under ``name``, when that is why there is no path
        to hand back.

    Attributes
    ----------
    assembly : str
        The assembly asked about.
    name : str
        The name that is not registered.
    registered : tuple of str
        The registered names, as they were passed in.
    offered : tuple of str
        The offered names, as they were passed in.
    broken : BrokenAnnotation or None
        The broken registration, or ``None`` when nothing of that name is on disk.

    Examples
    --------
    >>> raise AnnotationNotRegisteredError("hg38", "gencode_v50", [], ["gencode_v50"])
    Traceback (most recent call last):
    genome.io.gtf.AnnotationNotRegisteredError: "no annotation 'gencode_v50' ...
    """

    def __init__(
        self,
        assembly: str,
        name: str,
        registered: Iterable[str],
        offered: Iterable[str],
        *,
        broken: BrokenAnnotation | None = None,
    ) -> None:
        self.assembly = assembly
        self.name = name
        self.registered: tuple[str, ...] = tuple(registered)
        self.offered: tuple[str, ...] = tuple(offered)
        self.broken = broken
        if broken is not None:
            next_step = f"A broken registration for it is on disk: {broken.problem}"
        elif name in self.offered:
            next_step = (
                f"The annotation table offers it for {assembly!r}, so register it with "
                f"`{_register_command(assembly, name)}`."
            )
        else:
            next_step = (
                f"The table does not offer {name!r} for {assembly!r} either — it offers: "
                f"{_elide(self.offered) or '(none)'}. Register one of those by name, or a GTF "
                f"no row lists by path with Genome.register_gtf(<path>, {name!r})."
            )
        super().__init__(
            f"no annotation {name!r} is registered for {assembly!r}. Registered here: "
            f"{_elide(self.registered) or '(none)'}. {next_step}"
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
class BrokenAnnotation:
    """An annotation directory that is there and cannot be trusted as finished.

    What :func:`list_annotations` leaves out, said out loud. It is not a
    :class:`GtfAnnotation` and carries no file paths, because the whole point is that
    nothing vouches for the files: what it carries instead is why it cannot be trusted
    and the one command that makes it trustworthy again.

    Attributes
    ----------
    name : str
        The **Registered name** its directory is filed under.
    directory : pathlib.Path
        The annotation directory, whatever state it is in.
    problem : str
        What is wrong, in full — which files disagree or which are there with no record
        — ending in the ``repair`` below. This is
        :func:`~genome.io.completion.check_registration`'s own message, so re-registering
        the annotation says exactly what listing it says.
    repair : str
        The command that registers it again from scratch.

    Examples
    --------
    >>> from pathlib import Path
    >>> broken = BrokenAnnotation(
    ...     name="mine",
    ...     directory=Path("/data/genome/hg38/gtf/mine"),
    ...     problem="... holds files but no .completion.json ...",
    ...     repair="genome register-gtf hg38 /tmp/mine.gtf mine --force",
    ... )
    >>> broken.repair
    'genome register-gtf hg38 /tmp/mine.gtf mine --force'
    """

    name: str
    directory: Path
    problem: str
    repair: str


@dataclass(frozen=True)
class RegisteredAnnotation:
    """What registering one annotation produced: its record, and where that landed.

    :func:`register_annotation`'s answer and :func:`register_annotation_by_path`'s — what
    ``genome register-annotation`` and ``genome register-gtf`` print, and what their
    ``--json`` serializes. A :class:`GtfAnnotation` says where an annotation's two files
    are; this says what the run that wrote them did, which is the **Completion marker**
    itself, carried whole. Every question a surface then asks — the digest, the source,
    the files claimed, whether the chromosome names were actually checked — is answered
    from that one record rather than by reading the directory again.

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
        """Return this registration as the payload ``--json`` serializes.

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


@dataclass(frozen=True)
class AnnotationStatusRow:
    """One annotation, in whichever of its states it is: offered, registered, broken.

    One shape for all of them, so a reader never has to ask which fields a row has — a
    name the table does not list carries the table's columns as ``None``, and one nothing
    is wrong with carries the broken columns as ``None``. :attr:`registered` and
    :attr:`broken` are never both true: a registration nothing vouches for is not one.

    Attributes
    ----------
    name : str
        The **Registered name** this row is about.
    offered : bool
        Whether the annotation table lists it for this assembly.
    registered : bool
        Whether a record here vouches for it.
    broken : bool
        Whether its directory is here and cannot be trusted.
    default : bool
        The table's own default flag, ``False`` for a name no row lists.
    provider : str or None
        Who publishes it, from the table's row; ``None`` for an unlisted one.
    version : str or None
        The provider's release identifier; ``None`` for an unlisted one.
    url : str or None
        Where the table says its GTF is fetched from; ``None`` for an unlisted one.
    sha256 : str or None
        The digest the table pins; ``None`` when it pins none, and for an unlisted one.
    path : str or None
        The registered GTF's path, or ``None`` when it is not registered here.
    problem : str or None
        What is wrong, when :attr:`broken`; ``None`` otherwise.
    repair : str or None
        The command that registers it again from scratch, when :attr:`broken`.

    Examples
    --------
    >>> row = AnnotationStatusRow(
    ...     name="gencode_v50",
    ...     offered=True,
    ...     registered=False,
    ...     broken=False,
    ...     default=True,
    ...     provider="GENCODE",
    ...     version="v50",
    ...     url="https://example.org/gencode_v50.gtf.gz",
    ...     sha256=None,
    ...     path=None,
    ...     problem=None,
    ...     repair=None,
    ... )
    >>> row.as_json()["offered"]
    True
    """

    name: str
    offered: bool
    registered: bool
    broken: bool
    default: bool
    provider: str | None
    version: str | None
    url: str | None
    sha256: str | None
    path: str | None
    problem: str | None
    repair: str | None

    def as_json(self) -> dict[str, Any]:
        """Return this row as ``--json`` serializes it: every attribute above, in order.

        Returns
        -------
        dict
            The row's fields, under their own names.
        """
        return asdict(self)


@dataclass(frozen=True)
class AnnotationStatus:
    """What one assembly's table offers, set against what is registered on this machine.

    :meth:`AnnotationRegistry.status`'s answer, and what ``genome annotations`` prints.
    Two questions joined for one reader, with a third riding along because this is where
    anyone would look for it: a directory that cannot be trusted is ``broken`` rather than
    registered, and reporting one is the point — nothing here raises.

    Attributes
    ----------
    assembly : str
        The **Assembly** reported on.
    directory : pathlib.Path
        Its **Assembly dir**, whether or not anything is there.
    default_annotation : str or None
        The **Default annotation**'s name, or ``None`` when nothing decides one. It may
        name one nobody has registered here, which is a fresh machine's ordinary state.
    annotations : tuple of AnnotationStatusRow
        One row per name: the offered ones in table order, then anything on this disk
        that no row lists.

    Examples
    --------
    >>> from pathlib import Path
    >>> status = AnnotationStatus(
    ...     assembly="hg38",
    ...     directory=Path("/data/genome/hg38"),
    ...     default_annotation=None,
    ...     annotations=(),
    ... )
    >>> status.default_row is None
    True
    >>> status.as_json()["directory"]
    '/data/genome/hg38'
    """

    assembly: str
    directory: Path
    default_annotation: str | None
    annotations: tuple[AnnotationStatusRow, ...]

    @property
    def default_row(self) -> AnnotationStatusRow | None:
        """The **Default annotation**'s own row, or ``None`` when no row is about it.

        ``None`` covers both of the ways that happens, and a caller wanting to tell them
        apart reads :attr:`default_annotation` beside this: nothing decided a default, or
        one is decided and the table lists it under a name this disk knows nothing about.
        """
        return next((row for row in self.annotations if row.name == self.default_annotation), None)

    def as_json(self) -> dict[str, Any]:
        """Return this report as the payload ``--json`` serializes.

        Returns
        -------
        dict
            ``assembly``, the ``directory`` as text, the ``default_annotation`` name, and
            ``annotations`` as a list of :meth:`AnnotationStatusRow.as_json` rows.
        """
        return {
            "assembly": self.assembly,
            "directory": str(self.directory),
            "default_annotation": self.default_annotation,
            "annotations": [row.as_json() for row in self.annotations],
        }


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


def _register_command(assembly: str, name: str) -> str:
    """Return the command that registers ``name`` for ``assembly``."""
    return f"genome register-annotation {assembly} {name}"


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
    return f"{_register_command(assembly, name)} --force"


def _path_repair_call(source: Path, name: str) -> str:
    """Return the Python call that registers an unlisted GTF again from scratch.

    Given a directory and no assembly name there is no command to name: the caller who
    passed the path is the only one who knows which assembly it is for. So the repair
    names the call instead. Addressed by assembly name, it is a command a shell can run
    — see :func:`_path_repair_command`.
    """
    return f"register_gtf(<assembly dir>, {str(source)!r}, {name!r}, force=True)"


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


def list_annotations(assembly_dir: Path) -> dict[str, GtfAnnotation]:
    """Return the annotations registered under ``assembly_dir``, keyed by name.

    Registered means *a record is there and agrees with what is on disk* — never that a
    database file exists, which is true of a build killed half-way through as well as of
    a finished one. Anything else in the ``gtf/`` subtree is left out rather than raised
    over: listing is a question about this machine, and one unfinished annotation must
    not stop a genome from opening. Left out is not lost —
    :func:`list_broken_annotations` is where those go, and between the two every
    directory under ``gtf/`` is accounted for.

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


def list_broken_annotations(assembly_dir: Path, assembly: str) -> dict[str, BrokenAnnotation]:
    """Return the annotations under ``assembly_dir`` that cannot be trusted, keyed by name.

    The complement of :func:`list_annotations`. Registering an annotation over a
    directory like this raises (ADR-0007) and that is right, but a caller who never
    registers anything would otherwise never hear of it: a half-built annotation read as
    one nobody had fetched, and one no table row lists did not appear at all. So this
    reports rather than raises, and each entry carries the command that repairs it.

    A directory that is absent, or empty but for its working area, is a registration
    nobody has begun and is not broken — the same rule the registration path follows.

    ``assembly`` is needed for more than the message: whether the annotation table lists
    the name decides which command repairs it, since a listed one is re-fetched by name
    and an unlisted one has to be handed its GTF again.

    Parameters
    ----------
    assembly_dir : pathlib.Path
        The assembly directory whose ``gtf/`` subtree is inspected.
    assembly : str
        The assembly those annotations belong to, e.g. ``"hg38"``.

    Returns
    -------
    dict of str to BrokenAnnotation
        Registered name to what is wrong with it, in directory-name order. Empty when
        every annotation there is finished, which is the ordinary case.

    Examples
    --------
    >>> from pathlib import Path
    >>> list_broken_annotations(Path("/tmp/definitely-not-an-assembly"), "hg38")
    {}
    """
    root = _annotations_root(assembly_dir)
    if not root.is_dir():
        return {}
    offered = {record.name for record in list_annotation_metadata(assembly)}
    found: dict[str, BrokenAnnotation] = {}
    for directory in sorted(root.iterdir()):
        if not directory.is_dir():
            continue
        repair = _annotation_repair(directory, assembly=assembly, offered=offered)
        try:
            check_registration(directory, repair=repair)
        except RegistrationError as err:
            found[directory.name] = BrokenAnnotation(
                name=directory.name, directory=directory, problem=str(err), repair=repair
            )
    return found


def default_annotation(
    offered: Iterable[AnnotationMetadata],
    registered: Iterable[str],
    *,
    explicit: str | None = None,
) -> str | None:
    """Return the name of the **Default annotation**, or ``None`` when there is none.

    The whole rule, in one place, because two callers ask it: a
    :class:`~genome.genome.Genome` being opened, and :func:`annotation_status` reporting
    on an assembly nobody has opened. In order:

    1. an explicit choice, which is the caller overruling everything below it;
    2. the row the annotation table flags for this assembly, so everyone in the lab
       reaches for the same one without discussing it;
    3. the sole registered annotation, when exactly one is registered;
    4. otherwise none — a caller who did not choose between several is asked rather
       than guessed at.

    Only the first three lines of the table are consulted, never the disk: the name this
    returns may be one nothing has registered yet, which is the normal state of a fresh
    machine and is *not* an error. Where that name has to exist is
    :attr:`Genome.default_gtf_path <genome.genome.Genome.default_gtf_path>`.

    Parameters
    ----------
    offered : iterable of genome.metadata.AnnotationMetadata
        What the table offers for the assembly, in table order. The first flagged row
        wins, so a table that flags two names is read as naming the earlier one.
    registered : iterable of str
        The **Registered name**s on this machine.
    explicit : str, optional
        A name the caller chose, which wins over everything else. It is returned as
        given and is not checked against either list.

    Returns
    -------
    str or None
        The default annotation's name, or ``None`` when nothing decides one.

    Examples
    --------
    >>> from genome.metadata import AnnotationMetadata
    >>> row = AnnotationMetadata(
    ...     "hg38", "gencode_v50", "GENCODE", "v50", "https://example.org/g.gtf.gz", default=True
    ... )
    >>> default_annotation([row], [])                     # nothing registered yet
    'gencode_v50'
    >>> default_annotation([row], ["refseq_2023"], explicit="refseq_2023")
    'refseq_2023'
    >>> default_annotation([], ["refseq_2023"])           # no flag: the sole one stands
    'refseq_2023'
    >>> default_annotation([], ["refseq_2023", "mine"]) is None
    True
    """
    if explicit is not None:
        return explicit
    flagged = next((record.name for record in offered if record.default), None)
    if flagged is not None:
        return flagged
    names = list(registered)
    return names[0] if len(names) == 1 else None


class AnnotationRegistry:
    """One **Assembly**'s annotations, the state each is in, and the acts that add one.

    An annotation directory is *registered*, *broken*, *offered but not begun*, or nothing
    at all, and every useful question about one is a question about that four-way state:
    what may a caller name, what may it be handed the path of, which is the **Default
    annotation**, what does a surface print, what does a name nobody registered earn as an
    error. This settles all four **once**, at construction, and answers from that — so the
    state is assembled in one place rather than wherever it is needed.

    Bound to one assembly and carried, never re-derived: the **Assembly dir** comes in as
    an :class:`~genome.io.registration.AssemblyDir`, so a registry cannot file an
    annotation somewhere other than where the caller that built it is looking, and the
    ``chrom.sizes`` every GTF is checked against comes in beside it rather than being
    guessed from the layout.

    Reading is cheap and safe: nothing here is created, fetched or built by asking, an
    assembly with no directory at all answers emptily, and one broken annotation is
    reported rather than raised over. Only :meth:`register` and :meth:`register_path`
    write, and both fold what they wrote back in, so the four states stay current without
    reading the disk again.

    Parameters
    ----------
    assembly_dir : genome.io.registration.AssemblyDir
        The assembly this registry is for, and where its ``gtf/`` subtree lives.
    chrom_sizes : str or pathlib.Path, optional
        The assembly's ``chrom.sizes``, whose names every registered GTF's must be among.
        Defaults to the one the layout names, which is what an assembly prepared in place
        has; a caller that prepared it elsewhere passes the file it actually wrote. A path
        that is not there is *nothing to check against* rather than an error — an
        annotation may be registered before its assembly is.
    default : str, optional
        A **Default annotation** the caller chose, which wins over the table's flag and
        need not be registered. See :func:`default_annotation` for the whole rule.

    Attributes
    ----------
    assembly : str
        The assembly every annotation here belongs to.

    Examples
    --------
    >>> registry = AnnotationRegistry.locate("sacCer3", "/tmp/definitely-not-an-assembly")
    >>> registry.registered
    []
    >>> registry.default                       # the table's flag, registered or not
    'ensgene_v101'
    """

    def __init__(
        self,
        assembly_dir: AssemblyDir,
        *,
        chrom_sizes: str | Path | None = None,
        default: str | None = None,
    ) -> None:
        self._dir = assembly_dir
        self.assembly: str = assembly_dir.assembly
        self._chrom_sizes: Path = (
            assembly_dir.genome_files.chrom_sizes if chrom_sizes is None else Path(chrom_sizes)
        )
        self._registered: dict[str, GtfAnnotation] = list_annotations(assembly_dir.path)
        self._broken: dict[str, BrokenAnnotation] = list_broken_annotations(
            assembly_dir.path, self.assembly
        )
        self._offered: list[AnnotationMetadata] = list_annotation_metadata(self.assembly)
        self._default: str | None = default_annotation(
            self._offered, self._registered, explicit=default
        )

    @classmethod
    def locate(
        cls, assembly: str, cache_dir: str | Path | None = None, *, default: str | None = None
    ) -> AnnotationRegistry:
        """Return the registry for ``assembly``, wherever the layout says its files live.

        The assembly-addressed way in, and the only one the CLI has: a name and at most a
        directory override. :meth:`~genome.io.registration.AssemblyDir.locate` is where
        that override rule lives, and the ``chrom.sizes`` is the one that layout names.

        Parameters
        ----------
        assembly : str
            The assembly to open the registry of, e.g. ``"hg38"``.
        cache_dir : str or pathlib.Path, optional
            An explicit **Assembly dir**, overriding the **Data dir** layout.
        default : str, optional
            A **Default annotation** the caller chose.

        Returns
        -------
        AnnotationRegistry
            Its registry. Nothing is created and nothing is fetched.

        Examples
        --------
        >>> AnnotationRegistry.locate("hg38", "/tmp/definitely-not-an-assembly").registered
        []
        """
        return cls(AssemblyDir.locate(assembly, cache_dir), default=default)

    @property
    def registered(self) -> list[str]:
        """The **Registered name**s on this machine, in directory-name order.

        What is here, as against :attr:`offered`, which is what the lab supports, and
        :attr:`broken`, which is what is here and cannot be trusted.
        """
        return list(self._registered)

    @property
    def broken(self) -> list[BrokenAnnotation]:
        """The annotation directories here that cannot be trusted as finished.

        What :attr:`registered` leaves out, and between the two every directory under
        ``gtf/`` is accounted for. Each entry says what is wrong and names the one command
        that repairs it.
        """
        return list(self._broken.values())

    @property
    def offered(self) -> list[AnnotationMetadata]:
        """The annotation table's rows for this assembly, in table order.

        What the lab supports, whether or not anyone has registered it. Empty for an
        assembly the table offers nothing for, which is legal: it is a cross-reference
        rather than an allow-list (ADR-0003).
        """
        return list(self._offered)

    @property
    def default(self) -> str | None:
        """Name of the **Default annotation**, or ``None`` when nothing decides one.

        :func:`default_annotation`'s answer for this assembly, settled when the registry
        was built. It may name an annotation nobody has registered here — the normal state
        of a fresh machine — so it is :meth:`path` that says whether one exists. A default
        already decided is never displaced by a later registration.
        """
        return self._default

    def path(self, name: str) -> Path:
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
        AnnotationNotRegisteredError
            If nothing of that name is registered here. The four-way state decides what
            the message says next: the command that registers ``name`` when the table
            offers it, the path-based way in when it does not, and — for a directory of
            that name that is there and broken — the command that registers it again from
            scratch, so what is named is a command that runs rather than one that raises
            in turn.

        Examples
        --------
        >>> registry = AnnotationRegistry.locate("sacCer3", "/data/genome/sacCer3")
        >>> registry.path("ensgene_v101")              # doctest: +SKIP
        PosixPath('/data/genome/sacCer3/gtf/ensgene_v101/ensgene_v101.gtf')
        """
        annotation = self._registered.get(name)
        if annotation is None:
            raise AnnotationNotRegisteredError(
                self.assembly,
                name,
                self._registered,
                [record.name for record in self._offered],
                broken=self._broken.get(name),
            )
        return annotation.gtf

    def register(
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
        """Register the annotation the table lists for this assembly as ``name``.

        Naming an annotation is enough: where its GTF comes from and which digest it must
        match are the curated table's to know. The row's URL is fetched into the working
        area, the **unpacked** GTF is checked against the sha256 the row pins (ADR-0006) —
        so a GTF that is not the pinned one never reaches the annotation directory — the
        gffutils database is built, and the record is written last.

        Its chromosome names are checked too, against this registry's ``chrom.sizes`` and
        while the GTF is still in the working area: every name the GTF uses must be one the
        assembly carries, so an Ensembl-spelled GTF registered against a UCSC-spelled
        assembly fails in seconds rather than after the minutes the database build takes.
        The reverse is not required — an assembly may carry scaffolds the annotation never
        mentions. An assembly with no ``chrom.sizes`` yet has nothing to check against, and
        the record says so in ``details["chromosomes_checked"]`` — with
        ``details["chromosomes_unchecked_because"]`` saying whether that was for want of a
        ``chrom.sizes`` or because the caller stood the check down.

        An annotation that already has a valid record is returned silently: nothing is
        fetched, nothing is rebuilt and nothing is warned about. A directory that cannot be
        trusted — files with no record, or a record that disagrees with disk — **raises**,
        naming ``genome register-annotation <assembly> <name> --force`` (ADR-0007). That is
        what ``force=True`` is: it skips the question, keeps a GTF whose digest can be shown
        to be the pinned one, and fetches the source again when it cannot.

        Parameters
        ----------
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
            register an annotation whose mismatch you have inspected and accept; the record
            then says the check was stood down, rather than merely that it did not run.
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
            If the table lists no annotation ``name`` for this assembly; the message lists
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
        >>> AnnotationRegistry.locate("sacCer3").register(       # doctest: +SKIP
        ...     "ensgene_v101"
        ... )
        GtfAnnotation(name='ensgene_v101', ...)
        """
        row = metadata if metadata is not None else lookup_annotation(self.assembly, name)
        if row is None:
            offered = ", ".join(record.name for record in self._offered) or "(none)"
            raise ValueError(
                f"no annotation named {name!r} is listed for {self.assembly!r}. Listed for it: "
                f"{offered}. An annotation the table does not list is registered by path "
                f"instead — `{_register_gtf_command(self.assembly, '<path>', name)}`, or "
                f"Genome.register_gtf(<path>, {name!r}) from Python."
            )

        annotation = _annotation_files(self._dir.path, name)
        repair = _repair_command(self.assembly, name)
        if _already_registered(annotation.gtf.parent, force=force, repair=repair):
            return self._adopt(annotation)

        known = _assembly_chromosomes(self._chrom_sizes) if check_chromosomes else None
        digest = _proven_gtf(annotation.gtf, row.sha256)
        if digest is None:
            digest = _fetch_gtf(annotation, row, progressbar=progressbar, known=known)
        else:
            # Kept from a previous run, so it is already placed; there is nothing to undo.
            _reject_unknown_chromosomes(annotation.gtf, known, name=name)

        return self._adopt(
            _build_and_record(
                annotation,
                source_url=row.url,
                sha256=digest,
                details={
                    "provider": row.provider,
                    "version": row.version,
                    **_chromosome_check_details(known, requested=check_chromosomes),
                },
                disable_infer_genes=disable_infer_genes,
                disable_infer_transcripts=disable_infer_transcripts,
            )
        )

    def register_path(
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
        :meth:`register` is the way in for one it does. A gzipped (``.gz``) source is
        decompressed into the registered ``<name>.gtf``; a plain GTF is copied as-is. The
        digest recorded is of the placed GTF, since an unlisted annotation has no pinned
        digest to compare against.

        Its chromosome names are checked against this registry's ``chrom.sizes`` before
        anything is created, so a GTF that does not line up leaves the annotation directory
        exactly as it was found. Knowing the assembly is what buys that: the by-directory
        :func:`register_gtf` has to be handed the file or check nothing.

        Registering something already registered returns it silently, and a directory that
        cannot be trusted raises naming ``genome register-gtf <assembly> <gtf> <name>
        --force``, exactly as :meth:`register` does for a listed one.

        Parameters
        ----------
        gtf : str or pathlib.Path
            Path to the source GTF, plain or ``.gz``.
        name : str
            The **Registered name** to address it by, unique within the assembly.
        force : bool, default False
            Register again from scratch — the repair for a directory that raises.
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
            If the GTF names sequences the assembly does not carry.
        genome.io.completion.RegistrationError
            If the annotation's directory cannot be trusted as finished.

        Examples
        --------
        >>> AnnotationRegistry.locate("sacCer3").register_path(  # doctest: +SKIP
        ...     "custom.gtf.gz", "custom"
        ... )
        GtfAnnotation(name='custom', ...)
        """
        source = Path(gtf)
        return self._adopt(
            _register_gtf(
                self._dir.path,
                source,
                name,
                repair=_path_repair_command(self.assembly, shlex.quote(str(source)), name),
                force=force,
                chrom_sizes=self._chrom_sizes,
                check_chromosomes=check_chromosomes,
                disable_infer_genes=disable_infer_genes,
                disable_infer_transcripts=disable_infer_transcripts,
            )
        )

    def status(self) -> AnnotationStatus:
        """Report what this assembly's table offers against what is registered here.

        Two questions with two answers, joined for one reader: the table's rows say what
        the lab supports, the disk says what is on this machine, and every row carries
        which of the two it is. The command behind it is ``genome annotations``.

        A third answer rides along, because this is where a reader would look for it: a
        directory that cannot be trusted is ``broken`` rather than ``registered``. Nothing
        raises — reporting a broken annotation is the point, and one of them must not cost
        the rest.

        Returns
        -------
        AnnotationStatus
            The assembly, its directory, the **Default annotation**'s name, and one
            :class:`AnnotationStatusRow` per name — the offered ones in table order,
            followed by anything on this disk that no row lists.

        Examples
        --------
        >>> here = AnnotationRegistry.locate("sacCer3", "/tmp/definitely-not-an-assembly")
        >>> here.status().default_annotation
        'ensgene_v101'
        """
        rows: list[AnnotationStatusRow] = [
            _status_row(
                record.name, table_row=record, registered=self._registered, broken=self._broken
            )
            for record in self._offered
        ]
        listed = {record.name for record in self._offered}
        rows.extend(
            _status_row(name, table_row=None, registered=self._registered, broken=self._broken)
            for name in sorted((self._registered.keys() | self._broken.keys()) - listed)
        )
        return AnnotationStatus(
            assembly=self.assembly,
            directory=self._dir.path,
            default_annotation=self._default,
            annotations=tuple(rows),
        )

    def _adopt(self, annotation: GtfAnnotation) -> GtfAnnotation:
        """Fold a just-registered annotation into the four states, adopting it if alone.

        The sole-registered clause of the default rule, applied the moment it becomes
        true. A default already decided — the caller's choice, or the table's flag — is
        never displaced by one being registered. Registering over a broken directory is
        what repairs it, so the name stops being reported as broken here rather than only
        the next time the disk is read.
        """
        self._registered[annotation.name] = annotation
        self._broken.pop(annotation.name, None)
        if self._default is None and len(self._registered) == 1:
            self._default = annotation.name
        return annotation


def annotation_status(assembly: str, *, cache_dir: str | Path | None = None) -> AnnotationStatus:
    """Report what ``assembly``'s table offers against what is registered on this machine.

    :meth:`AnnotationRegistry.status` for an assembly named rather than opened, which is
    what ``genome annotations`` runs. Nothing is prepared, fetched, built or created to
    answer it — an assembly with no directory at all is the case it most needs to serve.

    Parameters
    ----------
    assembly : str
        The assembly to report on, e.g. ``"hg38"``.
    cache_dir : str or pathlib.Path, optional
        Override which assembly directory is inspected. Defaults to
        :func:`assembly_data_dir(assembly) <genome.io.download.assembly_data_dir>`.

    Returns
    -------
    AnnotationStatus
        The report :meth:`AnnotationRegistry.status` describes.

    Examples
    --------
    >>> annotation_status("sacCer3").default_annotation
    'ensgene_v101'
    """
    return AnnotationRegistry.locate(assembly, cache_dir).status()


def _status_row(
    name: str,
    *,
    table_row: AnnotationMetadata | None,
    registered: dict[str, GtfAnnotation],
    broken: dict[str, BrokenAnnotation],
) -> AnnotationStatusRow:
    """Return one :class:`AnnotationStatus` row, whichever of the three states it is in."""
    annotation = registered.get(name)
    broken_annotation = broken.get(name)
    return AnnotationStatusRow(
        name=name,
        offered=table_row is not None,
        registered=annotation is not None,
        broken=broken_annotation is not None,
        default=table_row.default if table_row is not None else False,
        provider=table_row.provider if table_row is not None else None,
        version=table_row.version if table_row is not None else None,
        url=table_row.url if table_row is not None else None,
        sha256=table_row.sha256 if table_row is not None else None,
        path=str(annotation.gtf) if annotation is not None else None,
        problem=broken_annotation.problem if broken_annotation is not None else None,
        repair=broken_annotation.repair if broken_annotation is not None else None,
    )


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

    :meth:`AnnotationRegistry.register_path` given a directory and no assembly name —
    which is the one thing a registry cannot be built from, and see it for what registering
    does. Having no assembly name costs this form exactly two things, and they are the
    whole difference. Its ``chrom.sizes`` cannot be derived, so it is passed or nothing is
    checked; and the repair a broken directory names can only be the Python call, since
    nobody but the caller knows which assembly that directory is for.

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
        register a GTF whose mismatch you have inspected and accept; the record
        then says the check was stood down, rather than merely that it did not run.
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
    return _register_gtf(
        assembly_dir,
        source,
        name,
        repair=_path_repair_call(source, name),
        force=force,
        chrom_sizes=chrom_sizes,
        check_chromosomes=check_chromosomes,
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )


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
    """Place, build and record the GTF at ``source``, as :func:`register_gtf` describes.

    ``repair`` is the only thing its two callers differ over: what a broken annotation
    directory tells the caller to run. Given an assembly name that is a command; given
    only a directory it can only be the Python call.
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

    :meth:`AnnotationRegistry.register` for an assembly given as a name and a directory
    rather than opened — see it for what registering does, what is checked and when this
    raises. The two agree about where the ``chrom.sizes`` is, since an assembly prepared
    in its own directory keeps it at ``<assembly_dir>/<assembly>.chrom.sizes``.

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
        A complete annotation record to use *instead of* the curated table's row.
    check_chromosomes : bool, default True
        Check the GTF's chromosome names against the assembly's.
    disable_infer_genes : bool, default True
        Do not reconstruct ``gene`` features from exon lines.
    disable_infer_transcripts : bool, default True
        Do not reconstruct ``transcript`` features from exon lines.

    Returns
    -------
    GtfAnnotation
        The registered annotation's name and its two file paths.

    Examples
    --------
    >>> from pathlib import Path
    >>> fetch_annotation(                                # doctest: +SKIP
    ...     Path("/data/genome/sacCer3"), "sacCer3", "ensgene_v101"
    ... )
    GtfAnnotation(name='ensgene_v101', ...)
    """
    return AnnotationRegistry(AssemblyDir(assembly=assembly, path=assembly_dir)).register(
        name,
        force=force,
        progressbar=progressbar,
        metadata=metadata,
        check_chromosomes=check_chromosomes,
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
    disable_infer_genes: bool = True,
    disable_infer_transcripts: bool = True,
) -> RegisteredAnnotation:
    """Register ``name`` for ``assembly`` and return the record of what that did.

    :meth:`AnnotationRegistry.register` addressed by assembly name, and answering with the
    record rather than the paths — the call ``genome register-annotation`` makes, and the
    one a script makes when it wants to serialize what happened. An annotation that is
    already registered is returned from its record without fetching anything.
    :func:`register_annotation_by_path` is the same shape for a GTF the table does not
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
    annotation = AnnotationRegistry.locate(assembly, cache_dir).register(
        name,
        force=force,
        progressbar=progressbar,
        metadata=metadata,
        check_chromosomes=check_chromosomes,
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )
    return _registration_payload(
        annotation, assembly=assembly, repair=_repair_command(assembly, name)
    )


def register_annotation_by_path(
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

    :meth:`AnnotationRegistry.register_path` addressed by assembly name, and answering
    with the record rather than the paths — the call ``genome register-gtf`` makes, and
    the way a script registers an annotation the curated table does not list and then
    serializes what happened. :func:`register_annotation` is the same shape for one the
    table does list.

    Naming the assembly buys the one thing the by-directory form cannot do for itself:
    its ``chrom.sizes`` is found rather than passed, so an unlisted GTF has its
    chromosome names checked by default, exactly as a listed one does. An assembly that
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
    >>> register_annotation_by_path(                     # doctest: +SKIP
    ...     "sacCer3", "custom.gtf.gz", "custom"
    ... )
    RegisteredAnnotation(assembly='sacCer3', directory=PosixPath('...'), record=...)
    """
    source = Path(gtf)
    annotation = AnnotationRegistry.locate(assembly, cache_dir).register_path(
        source,
        name,
        force=force,
        check_chromosomes=check_chromosomes,
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )
    return _registration_payload(
        annotation,
        assembly=assembly,
        repair=_path_repair_command(assembly, shlex.quote(str(source)), name),
    )


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
        "chromosomes_unchecked_because": _NO_CHROM_SIZES if requested else _CALLER_OVERRIDE,
    }


def _registration_payload(
    annotation: GtfAnnotation, *, assembly: str, repair: str
) -> RegisteredAnnotation:
    """Return the record a just-registered annotation left, as the value both ways in answer with.

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
