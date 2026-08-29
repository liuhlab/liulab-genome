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

There are two ways to add an annotation. By **name**: :meth:`AnnotationRegistry.register`
takes the name the curated annotation table lists for this assembly, fetches that row's
URL, checks the unpacked GTF against the sha256 the row pins (ADR-0006), builds the
database and writes the record. By **path**: :meth:`AnnotationRegistry.register_path` is
the escape hatch for a GTF no row lists — the caller says where the file is, and it is
placed, built and recorded the same way. :func:`register_annotation` and
:func:`register_gtf` are those same two addressed by assembly name and answering with the
record rather than the paths, one apiece, matching ``genome register-annotation`` and
``genome register-gtf`` exactly; both build a registry for the length of the call, so they
add no second code path. What they answer *with* —
:class:`~genome.io.results.RegisteredAnnotation`, and the two shapes
:meth:`AnnotationRegistry.status` reports in — is :mod:`genome.io.results`, so this module
changes when registering changes and not when the shape of an answer does.

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

**An annotation can also say which genes are in a category.**
:meth:`AnnotationRegistry.gene_list` and :meth:`AnnotationRegistry.gene_lists` answer from
the **Curated gene list** :mod:`genome.gene_list` ships for that annotation, addressed by
the same **Registered name** everything else here is. A **Merged annotation** answers per
contributor, read off its record, so its genes stay attributable to the component they came
from. Neither ever answers with an empty collection: an annotation that ships no list and
one whose list does not declare the category asked for raise errors of their own, since a
caller acts differently on those two facts and a silent zero is what that surface exists to
prevent.

**And which of its own gene ids a caller's stems name.**
:meth:`AnnotationRegistry.resolve_gene_ids` matches a **Gene id stem** — a gene id with
its version dropped — against the stem of every gene id in the **Annotation database**,
which makes it the first thing in this package to open the database registering an
annotation has always built. It answers with *every* gene id a stem names and never picks
one, and the stems that named nothing come back on the answer rather than being dropped.
It is general: nothing about it knows what the stems it is handed are a list *of*, and a
caller holding a few thousand of them resolves the lot in one pass.

**Its first caller is a published census.** :meth:`AnnotationRegistry.tf_gene_list`
resolves the **TF gene table** :mod:`genome.tf.gene` ships for this assembly's species
into the annotation's own gene ids, so the answer joins to a counts matrix with nothing
left to normalise. The species is read from the assembly's own metadata row and never
passed in, so one species' transcription factors cannot be asked for while holding
another species' assembly (ADR-0003);
nothing here decides what a transcription factor is, and the census's provenance travels
on the answer. Two absences raise rather than answering emptily, as the gene categories'
pair does: an assembly whose species has no census, and one nothing names a species for.

**Its second caller asks the other half of the same question.**
:meth:`AnnotationRegistry.tf_cofactor_list` resolves the **Cofactor table**
:mod:`genome.tf.cofactor` ships for the species the same way, and answers in the same
shape, so a caller who has read one answer has read both. It is a second caller of one
resolver and not a second crossing: the two differ in which table is read and in what a
row of it says, and in nothing else. The species is the assembly's own here too, and the
absences are the same pair — with worm answering here while the census half raises for it,
because a publisher assessed worm cofactors and none has released a worm TF census.

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
from collections.abc import Container, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any

import gffutils
import pooch

from genome.chimera import suffixed
from genome.gene_list import (
    CuratedGeneList,
    GeneCategoryNotDeclaredError,
    NoGeneCategoriesError,
    curated_annotations,
    curated_gene_list,
)
from genome.io import fetch
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
from genome.io.registration import ANNOTATIONS_SUBDIR, AssemblyDir, assembly_repair_command
from genome.io.results import (
    UNCHECKED_CALLER_OVERRIDE,
    UNCHECKED_NO_CHROM_SIZES,
    AnnotationStatus,
    AnnotationStatusRow,
    GeneList,
    GeneListSource,
    RegisteredAnnotation,
    ResolvedGeneIds,
    TFCofactor,
    TFCofactorList,
    TFGene,
    TFGeneList,
    annotation_register_command,
)
from genome.io.utils import ChecksumMismatchError, _gunzip, sha256_file
from genome.metadata import (
    AnnotationMetadata,
    assembly_metadata,
    list_annotation_metadata,
    lookup_annotation,
)
from genome.tf.cofactor import UNIFORM_COLUMNS as COFACTOR_UNIFORM_COLUMNS
from genome.tf.cofactor import (
    CofactorTable,
    cofactor_metadata,
    cofactor_species,
    cofactor_table,
)
from genome.tf.gene import (
    TRUE_CELL,
    UNIFORM_COLUMNS,
    TFGeneTable,
    census_metadata,
    census_species,
    species_slug,
    tf_gene_table,
)

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

#: What is asked of an **Annotation database** to resolve **Gene id stem**s: the id of
#: every gene feature, ascending, and nothing else. That id is the table's primary key —
#: so the ordering is an index walk — and for a GTF it is the ``gene_id`` gffutils keys
#: gene features by, which is the annotation's own spelling of its gene ids. It goes
#: through gffutils' own ``execute`` and is read a row at a time off the cursor: building
#: a feature object per row would parse an attribute blob per gene for a value already in
#: hand, and a GENCODE annotation has some 78,000 of them.
_GENE_IDS_QUERY = "SELECT id FROM features WHERE featuretype = 'gene' ORDER BY id"

#: What separates a gene id from its version — ``ENSG00000123456.7`` — and therefore what
#: a **Gene id stem** is everything before. An id carrying none is its own stem.
_VERSION_SEPARATOR = "."


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
                f"`{annotation_register_command(assembly, name)}`."
            )
        else:
            next_step = (
                f"The table does not offer {name!r} for {assembly!r} either — it offers: "
                f"{_elide(self.offered) or '(none)'}. Register one of those by name, or a GTF "
                f"no row lists by path with "
                f"genome.annotations.register_path(<path>, {name!r})."
            )
        super().__init__(
            f"no annotation {name!r} is registered for {assembly!r}. Registered here: "
            f"{_elide(self.registered) or '(none)'}. {next_step}"
        )


class NoGeneFeaturesError(LookupError):
    """An annotation's database holds no gene at all, so no gene id can be resolved.

    The absence a caller must never read as *this annotation carries none of my genes*.
    A GTF that declares only exons registers as exons alone —
    :meth:`AnnotationRegistry.register` leaves **Feature inference** off, and rightly, since
    it is gffutils' slow path and the publishers who matter declare their genes — so an
    annotation like that would answer every stem with *not found* while looking perfectly
    healthy. It says so instead, and names the argument that rebuilds it with the genes in.

    A :class:`LookupError`, as the other absences on this surface are.

    Parameters
    ----------
    annotation : str
        The **Registered name** that was asked about.
    assembly : str
        The **Assembly** it is registered for.

    Attributes
    ----------
    annotation : str
        The name asked about.
    assembly : str
        The assembly it is registered for.

    Examples
    --------
    >>> try:
    ...     raise NoGeneFeaturesError("mine", "tiny")
    ... except LookupError as error:
    ...     print("--infer-genes" in str(error))
    True
    """

    def __init__(self, annotation: str, assembly: str) -> None:
        self.annotation = annotation
        self.assembly = assembly
        super().__init__(
            f"the annotation {annotation!r} registered for {assembly!r} has no gene features "
            f"in its database, so there is nothing to resolve gene ids against. This is not "
            f"a gene it happens to lack: its GTF declares no gene lines at all, and "
            f"reconstructing them from the exons is off by default. Register it again with "
            f"that turned on — --infer-genes from a shell, disable_infer_genes=False from "
            f"Python — or register an annotation whose GTF declares its genes."
        )


class NoTFCensusError(LookupError):
    """No census has been published for this **Assembly**'s species, so none can answer.

    The first of the two absences a **TF gene list** has, and the one a caller must never
    read as *this species has no transcription factors*: nobody has published a census
    for it, which is a fact about the literature and not about the genome. Worm, yeast and
    *E. coli* land here today. The message names the species that do have one, since asking
    about one of those is the thing a caller can do instead.

    A :class:`LookupError`, so it may be caught together with :class:`UnknownSpeciesError`
    and still told apart — exactly as the **Curated gene list**'s two absences are.

    Parameters
    ----------
    assembly : str
        The **Assembly** asked about.
    species : str
        The species its metadata row names, which is what no census ships for.
    censused : iterable of str
        The species a census does ship for.

    Attributes
    ----------
    assembly : str
        The assembly asked about.
    species : str
        Its species.
    censused : tuple of str
        The species that do have a census.

    Examples
    --------
    >>> try:
    ...     raise NoTFCensusError("ce11", "Caenorhabditis elegans", ["Homo sapiens"])
    ... except LookupError as error:
    ...     print("Homo sapiens" in str(error))
    True
    """

    def __init__(self, assembly: str, species: str, censused: Iterable[str]) -> None:
        self.assembly = assembly
        self.species = species
        self.censused: tuple[str, ...] = tuple(censused)
        super().__init__(
            f"no TF census ships for {species!r}, which is the species the assembly table "
            f"names for {assembly!r} — so which of its genes are transcription factors is "
            f"unanswered here rather than answered with none. Censuses ship for: "
            f"{_elide(self.censused) or '(none)'}. Answering for {species!r} means shipping a "
            f"census for it, which scripts/build_tf_census.py writes."
        )


class NoCofactorTableError(LookupError):
    """Nobody has listed this **Assembly**'s species' cofactors, so no table can answer.

    :class:`NoTFCensusError`'s counterpart on the cofactor half, and the same kind of
    fact: nobody has published a cofactor list for this species, which is about the
    literature and not about the genome, and must never be read as *this species has no
    transcription cofactors*. Yeast and *E. coli* land here today.

    **Worm does not**, although :class:`NoTFCensusError` raises for the same assembly. A
    publisher assessed worm cofactors and none has released a worm TF census, so the two
    halves answer differently for one species — the publishers' shape, and not a defect.

    A :class:`LookupError`, so it may be caught together with :class:`UnknownSpeciesError`
    and still told apart, exactly as the census pair is.

    Parameters
    ----------
    assembly : str
        The **Assembly** asked about.
    species : str
        The species its metadata row names, which is what no table ships for.
    shipped_for : iterable of str
        The species a **Cofactor table** does ship for.

    Attributes
    ----------
    assembly : str
        The assembly asked about.
    species : str
        Its species.
    shipped_for : tuple of str
        The species that do have a cofactor table.

    Examples
    --------
    >>> try:
    ...     raise NoCofactorTableError("sacCer3", "Saccharomyces cerevisiae", ["Mus musculus"])
    ... except LookupError as error:
    ...     print("Mus musculus" in str(error))
    True
    """

    def __init__(self, assembly: str, species: str, shipped_for: Iterable[str]) -> None:
        self.assembly = assembly
        self.species = species
        self.shipped_for: tuple[str, ...] = tuple(shipped_for)
        super().__init__(
            f"no cofactor table ships for {species!r}, which is the species the assembly "
            f"table names for {assembly!r} — so which of its genes are transcription "
            f"cofactors is unanswered here rather than answered with none. Cofactor tables "
            f"ship for: {_elide(self.shipped_for) or '(none)'}. Answering for {species!r} "
            f"means shipping a table for it, which scripts/build_tf_cofactor.py writes."
        )


class UnknownSpeciesError(LookupError):
    """Nothing says what species this **Assembly** is, so no shipped table can be chosen.

    The second absence both gene-keyed halves have, and a different fact from
    :class:`NoTFCensusError` and :class:`NoCofactorTableError`: the question was not *has
    anyone published for this species* but *which species is this*, and nothing answered
    it. Two ways in, and neither is a mistake — a **Chimera** is more than one species by
    construction and nothing published answers for one, and an assembly the curated table
    does not list carries no species at all, which is the ordinary state of a free-form
    local key.

    The species is read from the assembly's own metadata and never passed in, so this is
    what a caller gets instead of quietly being handed another species' answer.

    One class for both halves because the fact is one fact: only what a caller can ask
    about instead differs, and ``shipped_table`` is what says which of them was asked for.

    Parameters
    ----------
    assembly : str
        The **Assembly** whose species nothing names.
    shipped_for : iterable of str
        The species the table asked for does ship for.
    shipped_table : str
        What could not be chosen, named in the message — ``"TF census"`` or ``"cofactor
        table"``. Keyword-only and required, so a message never says census over an
        answer that was about cofactors.

    Attributes
    ----------
    assembly : str
        The assembly asked about.
    shipped_for : tuple of str
        The species the table asked for does ship for.
    shipped_table : str
        What could not be chosen.

    Examples
    --------
    >>> try:
    ...     raise UnknownSpeciesError(
    ...         "ce11_ecHT115", ["Homo sapiens"], shipped_table="TF census"
    ...     )
    ... except LookupError as error:
    ...     print("Homo sapiens" in str(error))
    True
    """

    def __init__(self, assembly: str, shipped_for: Iterable[str], *, shipped_table: str) -> None:
        self.assembly = assembly
        self.shipped_table = shipped_table
        self.shipped_for: tuple[str, ...] = tuple(shipped_for)
        super().__init__(
            f"nothing says what species {assembly!r} is, so no {shipped_table} can be chosen "
            f"for it. The species comes from the assembly's own metadata row and is never "
            f"passed in, which is what makes asking about one species while holding another's "
            f"assembly impossible. A chimera is more than one species and no {shipped_table} "
            f"answers for one; anything else needs a row naming its species in the assembly "
            f"metadata table. One ships for: {_elide(self.shipped_for) or '(none)'}."
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
class _Contributor:
    """One annotation whose **Curated gene list** may answer for another, and for what.

    A plain annotation contributes itself, against the assembly it is registered for. A
    **Merged annotation** contributes one of these per entry of its completion record's
    ``merged_from``, each against its own **Component** — which is the assembly whose
    genes those are, and therefore the one its curated list has to name.
    """

    component: str | None
    assembly: str
    annotation: str


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
        return self._annotation(name).gtf

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
                f"genome.annotations.register_path(<path>, {name!r}) from Python."
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
        exactly as it was found. Knowing the assembly is what buys that: the file is found
        rather than passed, so an unlisted GTF is held to the same check a listed one gets.

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

    def gene_list(self, category: str, name: str | None = None) -> GeneList:
        """Return the genes one registered annotation puts in ``category``.

        The genes come from the **Curated gene list** shipped for that annotation and
        never from the GTF's own biotype attribute, which is spelled two ways across four
        publishers, carries three taxonomies that do not agree, and is absent altogether
        from some annotations (ADR-0011). Nothing here knows a category vocabulary: which
        categories exist is what the curated list declares.

        The annotation must be **registered here** — it is resolved through :meth:`path`,
        so an unregistered name earns the error that names the command registering it.
        The curated list is then held to the assembly it was curated against, since a name
        is unique only within its assembly and a list found by name alone is not yet known
        to be about this reference.

        For a **Merged annotation** the record's ``merged_from`` says who contributed, and
        each contributor's own curated list answers for its own **Component**: the result
        carries one source per contributor that declares the category, so a caller counting
        worm ribosomal RNA can drop the *E. coli* entry. A contributor that does not
        declare it is simply absent — a bacterium has no mitochondria, and that is not a
        failure.

        **There is no empty answer.** An annotation nothing ships a list for, and one whose
        list does not declare this category, are different facts and each raises an error
        of its own.

        Parameters
        ----------
        category : str
            The **Gene category** to ask for, as the curated list spells it — ``"rRNA"``,
            ``"Mt_rRNA"``.
        name : str, optional
            The **Registered name** to ask about. Omitted, this assembly's **Default
            annotation** answers.

        Returns
        -------
        genome.io.results.GeneList
            The category, its gene ids, and one
            :class:`~genome.io.results.GeneListSource` per contributing curated list.

        Raises
        ------
        ValueError
            If ``name`` is omitted and no **Default annotation** is decided; the message
            names the argument that chooses one.
        AnnotationNotRegisteredError
            If nothing of that name is registered here.
        genome.gene_list.NoGeneCategoriesError
            If no curated list ships for that annotation — nothing can be asked of it,
            which is not the same answer as its having no genes in this category.
        genome.gene_list.GeneCategoryNotDeclaredError
            If it declares categories and not this one; the message lists the ones it does.
        genome.gene_list.GeneListAssemblyMismatchError
            If the curated list found under that name was curated against another
            assembly, in which case it must not answer here.

        Examples
        --------
        >>> registry = AnnotationRegistry.locate("ce11")      # doctest: +SKIP
        >>> registry.gene_list("rRNA").gene_ids[:2]           # doctest: +SKIP
        ['WBGene00004512', 'WBGene00004513']
        >>> [source.component for source in registry.gene_list("rRNA").sources]  # doctest: +SKIP
        [None]
        """
        resolved, contributors = self._gene_list_contributors(name)
        return self._category_answer(resolved, category, contributors)

    def gene_lists(self, name: str | None = None) -> tuple[GeneList, ...]:
        """Return every **Gene category** one registered annotation declares.

        :meth:`gene_list` for all of them at once, in the order the curated lists spell
        them — and for a **Merged annotation**, each contributor's own order, contributors
        first-listed first. Everything :meth:`gene_list` says about resolution, the
        assembly guard and attribution holds here.

        **Never an empty tuple.** An annotation that declares nothing raises rather than
        answering emptily, which is the whole distinction this surface exists to keep: a
        caller that got ``()`` could not tell *no categories are declared* from *every
        category is empty*, and no declared category is ever empty.

        Parameters
        ----------
        name : str, optional
            The **Registered name** to ask about. Omitted, this assembly's **Default
            annotation** answers.

        Returns
        -------
        tuple of genome.io.results.GeneList
            One entry per declared category, in declaration order. Never empty.

        Raises
        ------
        ValueError
            If ``name`` is omitted and no **Default annotation** is decided.
        AnnotationNotRegisteredError
            If nothing of that name is registered here.
        genome.gene_list.NoGeneCategoriesError
            If no curated list ships for that annotation.
        genome.gene_list.GeneListAssemblyMismatchError
            If a curated list found by name was curated against another assembly.

        Examples
        --------
        >>> registry = AnnotationRegistry.locate("hg38")      # doctest: +SKIP
        >>> [answer.category for answer in registry.gene_lists()]   # doctest: +SKIP
        ['rRNA', 'rRNA_pseudogene', 'Mt_rRNA']
        """
        resolved, contributors = self._gene_list_contributors(name)
        return tuple(
            self._category_answer(resolved, category, contributors)
            for category in _declared_categories(contributors)
        )

    def resolve_gene_ids(self, stems: Iterable[str], name: str | None = None) -> ResolvedGeneIds:
        """Return the gene ids one registered annotation carries for each **Gene id stem**.

        A stem is a gene id with its version dropped — ``ENSG00000123456`` for
        ``ENSG00000123456.7`` — which is how every published table keyed by gene arrives,
        and never how a GENCODE **Annotation** spells the same gene. This is the crossing:
        every gene id in the **Annotation database** is reduced to its own stem, and a stem
        answers with every gene id that reduced to it. An id carrying no version is its own
        stem, so an annotation whose ids were never versioned — WormBase's, SGD's — resolves
        each of its genes to itself and is untouched by an Ensembl-shaped assumption.

        **Every id, and never a chosen one.** One stem naming two gene ids is not a
        malformed annotation: ``gencode_v50lift37`` has nine such stems, eight of them
        pseudoautosomal genes carrying a ``_PAR_Y`` copy, and a resolver taking the first
        would hand back the X copy of a Y gene without saying it had chosen. So the answer
        is a mapping to *all* of them, ascending.

        **Nothing is dropped.** Stems this annotation carries no gene for come back in
        :attr:`~genome.io.results.ResolvedGeneIds.unresolved`, so a caller resolving a few
        thousand at once can see which of them this annotation does not have rather than
        counting the answer and wondering.

        The annotation must be **registered here**, and a **Merged annotation** is read
        exactly as any other: it has one database of its own, holding both components' gene
        features under the components' own gene ids — a merge rewrites seqnames and never a
        ``gene_id`` — so a stem naming a gene in each component answers with both, which is
        the same rule as the pseudoautosomal one and needs no attribution to apply it.

        One indexed pass over the database's gene features answers the whole call, however
        many stems it was handed; nothing reads the GTF, and no annotation is held in
        memory.

        Parameters
        ----------
        stems : iterable of str
            The **Gene id stem**s to resolve, in the order they should come back. Repeats
            are asked once. Pass them all at once: the cost is the pass, not the stem.
        name : str, optional
            The **Registered name** to resolve against. Omitted, this assembly's **Default
            annotation** answers.

        Returns
        -------
        genome.io.results.ResolvedGeneIds
            The stems that named gene ids, mapped to every id each names, and the stems
            that named none.

        Raises
        ------
        ValueError
            If ``name`` is omitted and no **Default annotation** is decided; the message
            names the argument that chooses one.
        AnnotationNotRegisteredError
            If nothing of that name is registered here.
        NoGeneFeaturesError
            If its database holds no gene at all — every stem would otherwise resolve to
            nothing, which is a different fact and must not be reported as this one.

        Examples
        --------
        >>> registry = AnnotationRegistry.locate("hg19")             # doctest: +SKIP
        >>> answer = registry.resolve_gene_ids(                      # doctest: +SKIP
        ...     ["ENSG00000182378", "ENSG00000141510"], "gencode_v50lift37"
        ... )
        >>> answer.resolved["ENSG00000182378"]                       # doctest: +SKIP
        ('ENSG00000182378.14', 'ENSG00000182378.14_PAR_Y')
        """
        resolved_name = self._named(name)
        annotation = self._annotation(resolved_name)
        # Asked once each and answered in the order they arrived, so a caller can read its
        # own list against the answer.
        asked = tuple(dict.fromkeys(stems))
        if not asked:
            return ResolvedGeneIds(
                assembly=self.assembly, annotation=resolved_name, resolved={}, unresolved=()
            )
        found, any_genes = _gene_ids_by_stem(annotation.db, frozenset(asked))
        if not any_genes:
            raise NoGeneFeaturesError(resolved_name, self.assembly)
        return ResolvedGeneIds(
            assembly=self.assembly,
            annotation=resolved_name,
            resolved={stem: tuple(found[stem]) for stem in asked if stem in found},
            unresolved=tuple(stem for stem in asked if stem not in found),
        )

    def tf_gene_list(
        self, name: str | None = None, *, include_rejected: bool = False
    ) -> TFGeneList:
        """Return the genes a census judges transcription factors, in this annotation's ids.

        The **TF gene table** :mod:`genome.tf.gene` ships for this assembly's species, met
        with one registered annotation: every **Gene id stem** the census is keyed by is
        resolved through :meth:`resolve_gene_ids` into the gene ids this annotation
        actually spells, so the answer joins to a counts matrix with nothing left for the
        caller to normalise. A stem naming two gene ids answers with both, and the stems
        this annotation carries no gene for ride back on
        :attr:`~genome.io.results.TFGeneList.unresolved` rather than being dropped.

        **The species is the assembly's own** — the curated metadata table's, read here and
        never passed in, so asking for human transcription factors while holding a mouse
        assembly is not expressible (ADR-0003).

        **Nothing here decides what a transcription factor is.** The verdict is the
        census's, and the answer carries the publisher, version and PubMed id that reached
        it. Assessed-positive by default, because the common case is not 2,765 rows to
        filter down to 1,639; ``include_rejected`` widens to the genes the census assessed
        and turned down, which is the whole census and still not every gene there is — one
        it never assessed is absent from both answers, and that is a third fact.

        Wanting only ``Known motif``, or wanting ``Inferred motif`` included, is a re-filter
        on the **TF assessment** each gene already carries in
        :attr:`~genome.io.results.TFGene.judgements` rather than a second flag here.

        Parameters
        ----------
        name : str, optional
            The **Registered name** to answer in the gene ids of. Omitted, this assembly's
            **Default annotation** answers.
        include_rejected : bool, default False
            Carry the genes the census assessed and judged *not* to be transcription
            factors as well, each saying so in
            :attr:`~genome.io.results.TFGene.is_tf`. A census that records no rejections
            answers the same either way.

        Returns
        -------
        genome.io.results.TFGeneList
            The genes, the census's provenance, and the stems that resolved to nothing.

        Raises
        ------
        UnknownSpeciesError
            If nothing names this assembly's species — a chimera, or an assembly no
            curated row lists.
        NoTFCensusError
            If no census ships for that species; the message names the ones that do.
        ValueError
            If ``name`` is omitted and no **Default annotation** is decided.
        AnnotationNotRegisteredError
            If nothing of that name is registered here.
        NoGeneFeaturesError
            If its database holds no gene at all.

        Examples
        --------
        >>> registry = AnnotationRegistry.locate("hg38")              # doctest: +SKIP
        >>> answer = registry.tf_gene_list("gencode_v50")             # doctest: +SKIP
        >>> answer.genes[0].symbol, answer.genes[0].dbd_family        # doctest: +SKIP
        ('TFAP2A', 'AP-2')
        >>> answer.provenance.publisher                               # doctest: +SKIP
        'Lambert et al. 2018'
        """
        species = assembly_metadata(self.assembly).species
        if species is None:
            raise UnknownSpeciesError(self.assembly, _censused_species(), shipped_table="TF census")
        census = tf_gene_table(species)
        if census is None:
            raise NoTFCensusError(self.assembly, species, _censused_species())
        resolved = self.resolve_gene_ids(
            census.gene_id_stems if include_rejected else census.assessed_positive, name
        )
        return TFGeneList(
            assembly=self.assembly,
            annotation=resolved.annotation,
            species=species,
            provenance=census.provenance,
            genes=_tf_genes(census, resolved),
            unresolved=resolved.unresolved,
        )

    def tf_cofactor_list(self, name: str | None = None) -> TFCofactorList:
        """Return the genes a publisher lists as cofactors, in this annotation's ids.

        The **Cofactor table** :mod:`genome.tf.cofactor` ships for this assembly's
        species, met with one registered annotation, exactly as :meth:`tf_gene_list` meets
        a census: every **Gene id stem** the table is keyed by is resolved through
        :meth:`resolve_gene_ids` into the gene ids this annotation actually spells, so the
        answer joins to a counts matrix with nothing left for the caller to normalise. A
        stem naming two gene ids answers with both, and the stems this annotation carries
        no gene for ride back on
        :attr:`~genome.io.results.TFCofactorList.unresolved` rather than being dropped.

        **The species is the assembly's own** — the curated metadata table's, read here
        and never passed in, so asking for human cofactors while holding a mouse assembly
        is not expressible (ADR-0003).

        **Membership is this package's and classification is each publisher's.** A table
        built from two publishers is a union nobody else published (ADR-0016), which is
        why each entry says which publisher listed the gene in
        :attr:`~genome.io.results.TFCofactor.source` and keeps every publisher's own
        vocabulary under that publisher's namespaced column name in
        :attr:`~genome.io.results.TFCofactor.classifications`. A ``source`` of ``both``
        asserts agreement on membership only, never on classification.

        There is no widening flag here, as there is on :meth:`tf_gene_list`: no publisher
        shipping today releases a rejected set, so the table's listed genes are the whole
        table. A source that did record rejections would ship them, and they would be left
        out here rather than needing a second argument.

        Parameters
        ----------
        name : str, optional
            The **Registered name** to answer in the gene ids of. Omitted, this assembly's
            **Default annotation** answers.

        Returns
        -------
        genome.io.results.TFCofactorList
            The cofactors, the publishers' provenance, and the stems that resolved to
            nothing.

        Raises
        ------
        UnknownSpeciesError
            If nothing names this assembly's species — a chimera, or an assembly no
            curated row lists.
        NoCofactorTableError
            If no cofactor table ships for that species; the message names the ones that
            do. Worm has a table although no TF census covers it, so this and
            :class:`NoTFCensusError` do not raise for the same set of assemblies.
        ValueError
            If ``name`` is omitted and no **Default annotation** is decided.
        AnnotationNotRegisteredError
            If nothing of that name is registered here.
        NoGeneFeaturesError
            If its database holds no gene at all.

        Examples
        --------
        >>> registry = AnnotationRegistry.locate("mm39")                    # doctest: +SKIP
        >>> answer = registry.tf_cofactor_list("gencode_vM39")              # doctest: +SKIP
        >>> first = answer.cofactors[0]                                     # doctest: +SKIP
        >>> first.symbol, first.classifications["animaltfdb_category"]      # doctest: +SKIP
        ('Scmh1', 'Other Cofactors')
        >>> answer.provenance.sources[0].publisher                          # doctest: +SKIP
        'AnimalTFDB'
        """
        species = assembly_metadata(self.assembly).species
        if species is None:
            raise UnknownSpeciesError(
                self.assembly, _cofactor_species(), shipped_table="cofactor table"
            )
        table = cofactor_table(species)
        if table is None:
            raise NoCofactorTableError(self.assembly, species, _cofactor_species())
        resolved = self.resolve_gene_ids(table.cofactor_stems, name)
        return TFCofactorList(
            assembly=self.assembly,
            annotation=resolved.annotation,
            species=species,
            provenance=table.provenance,
            cofactors=_tf_cofactors(table, resolved),
            unresolved=resolved.unresolved,
        )

    def _annotation(self, name: str) -> GtfAnnotation:
        """Return the files of the annotation registered as ``name``, or raise for it.

        Every question about one registered annotation comes through here, so a name
        nothing registered earns one error with one next action wherever it was asked.
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
        return annotation

    def _gene_list_contributors(
        self, name: str | None
    ) -> tuple[str, tuple[tuple[_Contributor, CuratedGeneList], ...]]:
        """Resolve the annotation asked about and the curated lists that answer for it.

        Everything both gene-category questions share: which annotation is meant, that it
        is registered here, who contributed to it, whose curated list may speak for each
        contributor, and that each of those was curated against the right assembly. A
        contributor no list ships for is left out rather than raised over — one component
        of a merge carrying no curated list must not silence the others — and only when
        *no* contributor has one is there nothing to answer with.
        """
        resolved = self._named(name)
        contributors = _contributors_of(self.path(resolved).parent, resolved, self.assembly)
        answering: list[tuple[_Contributor, CuratedGeneList]] = []
        for contributor in contributors:
            listed = curated_gene_list(contributor.annotation)
            if listed is None:
                continue
            listed.check_assembly(contributor.assembly)
            answering.append((contributor, listed))
        if not answering:
            merged = [one.annotation for one in contributors if one.component is not None]
            raise NoGeneCategoriesError(
                resolved, self.assembly, curated_annotations(), contributors=merged
            )
        return resolved, tuple(answering)

    def _category_answer(
        self,
        annotation: str,
        category: str,
        contributors: Sequence[tuple[_Contributor, CuratedGeneList]],
    ) -> GeneList:
        """Assemble one category's answer, or say that nobody declared it."""
        sources: list[GeneListSource] = []
        for contributor, listed in contributors:
            declared = listed.categories.get(category)
            if declared is None:
                continue
            sources.append(
                GeneListSource(
                    component=contributor.component,
                    annotation=contributor.annotation,
                    description=declared.description,
                    source=declared.source,
                    gene_ids=declared.gene_ids,
                )
            )
        if not sources:
            raise GeneCategoryNotDeclaredError(
                annotation, self.assembly, category, _declared_categories(contributors)
            )
        return GeneList(
            assembly=self.assembly,
            annotation=annotation,
            category=category,
            sources=tuple(sources),
        )

    def _named(self, name: str | None) -> str:
        """Return the annotation a caller meant, which is the default when they named none."""
        if name is not None:
            return name
        if self._default is None:
            raise ValueError(
                f"no annotation was named and {self.assembly!r} has no default one to fall "
                f"back on, so there is nothing to ask about. Name one with the annotation "
                f"argument — gene_list(<category>, <name>) here, annotation=<name> from the "
                f"assembly-addressed functions, --annotation <name> from a shell — or decide "
                f"the assembly's default once with Genome({self.assembly!r}, "
                f"default_gtf=<name>). genome.annotations.registered says what is registered "
                f"here."
            )
        return self._default

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


def gene_list(
    assembly: str,
    category: str,
    *,
    annotation: str | None = None,
    cache_dir: str | Path | None = None,
) -> GeneList:
    """Return the genes ``assembly``'s annotation puts in ``category``.

    :meth:`AnnotationRegistry.gene_list` for an assembly named rather than opened, which
    is what ``genome gene-list`` runs. A registry is built for the length of the call, so
    there is no second code path. Nothing is prepared, fetched or built to answer it.

    Parameters
    ----------
    assembly : str
        The assembly to ask about, e.g. ``"ce11"``.
    category : str
        The **Gene category**, as the curated list spells it.
    annotation : str, optional
        The **Registered name** to ask about; the **Default annotation** when omitted.
    cache_dir : str or pathlib.Path, optional
        Override which assembly directory is inspected, as
        :func:`annotation_status` takes it.

    Returns
    -------
    genome.io.results.GeneList
        The answer :meth:`AnnotationRegistry.gene_list` describes.

    Raises
    ------
    ValueError
        If ``annotation`` is omitted and no **Default annotation** is decided.
    AnnotationNotRegisteredError
        If that annotation is not registered here.
    genome.gene_list.NoGeneCategoriesError
        If no curated gene list ships for it.
    genome.gene_list.GeneCategoryNotDeclaredError
        If it declares categories and not this one.

    Examples
    --------
    >>> gene_list("ce11", "rRNA").category                   # doctest: +SKIP
    'rRNA'
    """
    return AnnotationRegistry.locate(assembly, cache_dir).gene_list(category, annotation)


def gene_lists(
    assembly: str, *, annotation: str | None = None, cache_dir: str | Path | None = None
) -> tuple[GeneList, ...]:
    """Return every **Gene category** ``assembly``'s annotation declares.

    :meth:`AnnotationRegistry.gene_lists` addressed by assembly name — what ``genome
    gene-categories`` runs, built the same way :func:`gene_list` is. Never an empty tuple:
    an annotation that declares nothing raises instead.

    Parameters
    ----------
    assembly : str
        The assembly to ask about, e.g. ``"hg38"``.
    annotation : str, optional
        The **Registered name** to ask about; the **Default annotation** when omitted.
    cache_dir : str or pathlib.Path, optional
        Override which assembly directory is inspected.

    Returns
    -------
    tuple of genome.io.results.GeneList
        One entry per declared category, in declaration order. Never empty.

    Raises
    ------
    ValueError
        If ``annotation`` is omitted and no **Default annotation** is decided.
    AnnotationNotRegisteredError
        If that annotation is not registered here.
    genome.gene_list.NoGeneCategoriesError
        If no curated gene list ships for it.

    Examples
    --------
    >>> [answer.category for answer in gene_lists("hg38")]    # doctest: +SKIP
    ['rRNA', 'rRNA_pseudogene', 'Mt_rRNA']
    """
    return AnnotationRegistry.locate(assembly, cache_dir).gene_lists(annotation)


def tf_gene_list(
    assembly: str,
    *,
    annotation: str | None = None,
    include_rejected: bool = False,
    cache_dir: str | Path | None = None,
) -> TFGeneList:
    """Return the genes a published census judges transcription factors in ``assembly``.

    :meth:`AnnotationRegistry.tf_gene_list` for an assembly named rather than opened,
    built the way :func:`gene_list` is: a registry for the length of the call, so a shell
    surface over it adds no second code path. Nothing is prepared, fetched or built to
    answer it, and the census is read from inside the package.

    Parameters
    ----------
    assembly : str
        The assembly to ask about, e.g. ``"hg38"``. Its own metadata row names the
        species, which is what selects the census.
    annotation : str, optional
        The **Registered name** to answer in the gene ids of; the **Default annotation**
        when omitted.
    include_rejected : bool, default False
        Carry the genes the census assessed and turned down as well.
    cache_dir : str or pathlib.Path, optional
        Override which assembly directory is inspected, as :func:`gene_list` takes it.

    Returns
    -------
    genome.io.results.TFGeneList
        The answer :meth:`AnnotationRegistry.tf_gene_list` describes.

    Raises
    ------
    UnknownSpeciesError
        If nothing names the assembly's species.
    NoTFCensusError
        If no census ships for that species.
    ValueError
        If ``annotation`` is omitted and no **Default annotation** is decided.
    AnnotationNotRegisteredError
        If that annotation is not registered here.
    NoGeneFeaturesError
        If its database holds no gene at all.

    Examples
    --------
    >>> tf_gene_list("hg38").provenance.publisher            # doctest: +SKIP
    'Lambert et al. 2018'
    """
    return AnnotationRegistry.locate(assembly, cache_dir).tf_gene_list(
        annotation, include_rejected=include_rejected
    )


def tf_cofactor_list(
    assembly: str,
    *,
    annotation: str | None = None,
    cache_dir: str | Path | None = None,
) -> TFCofactorList:
    """Return the genes a publisher lists as transcription cofactors in ``assembly``.

    :meth:`AnnotationRegistry.tf_cofactor_list` for an assembly named rather than opened,
    built the way :func:`tf_gene_list` is: a registry for the length of the call, so a
    shell surface over it adds no second code path. Nothing is prepared, fetched or built
    to answer it, and the table is read from inside the package.

    Parameters
    ----------
    assembly : str
        The assembly to ask about, e.g. ``"mm39"``. Its own metadata row names the
        species, which is what selects the table.
    annotation : str, optional
        The **Registered name** to answer in the gene ids of; the **Default annotation**
        when omitted.
    cache_dir : str or pathlib.Path, optional
        Override which assembly directory is inspected, as :func:`tf_gene_list` takes it.

    Returns
    -------
    genome.io.results.TFCofactorList
        The answer :meth:`AnnotationRegistry.tf_cofactor_list` describes.

    Raises
    ------
    UnknownSpeciesError
        If nothing names the assembly's species.
    NoCofactorTableError
        If no cofactor table ships for that species.
    ValueError
        If ``annotation`` is omitted and no **Default annotation** is decided.
    AnnotationNotRegisteredError
        If that annotation is not registered here.
    NoGeneFeaturesError
        If its database holds no gene at all.

    Examples
    --------
    >>> tf_cofactor_list("mm39").species                     # doctest: +SKIP
    'Mus musculus'
    """
    return AnnotationRegistry.locate(assembly, cache_dir).tf_cofactor_list(annotation)


def _censused_species() -> tuple[str, ...]:
    """Return the species a census ships for, in the spelling an assembly's row uses.

    What both absences name as the thing a caller can ask about instead. The shipped files
    are what is enumerated, since one is what makes a species answerable; the provenance
    table beside them is only read for the publisher's own spelling, so a census shipping
    without a row is still named — badly, as its slug, which is the state
    :func:`~genome.tf.gene.census.tf_gene_table` raises over anyway.
    """
    named = {species_slug(record.species): record.species for record in census_metadata()}
    return tuple(named.get(slug, slug) for slug in census_species())


def _tf_genes(census: TFGeneTable, resolved: ResolvedGeneIds) -> tuple[TFGene, ...]:
    """Return one entry per resolved stem, carrying the census's own row for that gene.

    The census's four uniform columns become fields and everything after them stays under
    the publisher's own name, because beyond those four no two censuses carry the same
    columns and nothing here compares one publisher's with another's (ADR-0014). The order
    is the census's own row order, which is the order the stems were asked about.
    """
    publisher_columns = census.columns[len(UNIFORM_COLUMNS) :]
    rows = {row[0]: row for row in census.rows}
    genes: list[TFGene] = []
    for stem, gene_ids in resolved.resolved.items():
        cells = dict(zip(census.columns, rows[stem], strict=True))
        genes.append(
            TFGene(
                gene_id_stem=stem,
                gene_ids=gene_ids,
                symbol=cells["symbol"],
                is_tf=cells["is_tf"] == TRUE_CELL,
                dbd_family=cells["dbd_family"],
                judgements={name: cells[name] for name in publisher_columns},
            )
        )
    return tuple(genes)


def _cofactor_species() -> tuple[str, ...]:
    """Return the species a **Cofactor table** ships for, in an assembly row's spelling.

    :func:`_censused_species` for the cofactor half, and answering the same question for
    the same reason: it is what both absences name as the thing a caller can ask about
    instead. The shipped files are what is enumerated, since one is what makes a species
    answerable; the provenance table beside them is only read for the publisher's own
    spelling, so a table shipping without a row is still named — badly, as its slug,
    which is the state :func:`~genome.tf.cofactor.table.cofactor_table` raises over
    anyway.
    """
    named = {species_slug(record.species): record.species for record in cofactor_metadata()}
    return tuple(named.get(slug, slug) for slug in cofactor_species())


def _tf_cofactors(table: CofactorTable, resolved: ResolvedGeneIds) -> tuple[TFCofactor, ...]:
    """Return one entry per resolved stem, carrying the table's own row for that gene.

    :func:`_tf_genes` for the cofactor half. The table's four uniform columns become
    fields and everything after them stays under the publisher's own namespaced name,
    because no two publishers carry the same columns and nothing here compares one's
    vocabulary with another's (ADR-0014). The order is the table's own row order, which
    is the order the stems were asked about.
    """
    publisher_columns = table.columns[len(COFACTOR_UNIFORM_COLUMNS) :]
    rows = {row[0]: row for row in table.rows}
    cofactors: list[TFCofactor] = []
    for stem, gene_ids in resolved.resolved.items():
        cells = dict(zip(table.columns, rows[stem], strict=True))
        cofactors.append(
            TFCofactor(
                gene_id_stem=stem,
                gene_ids=gene_ids,
                symbol=cells["symbol"],
                is_cofactor=cells["is_cofactor"] == TRUE_CELL,
                # ``source`` is one of a closed vocabulary, checked as the file is read,
                # so ``or ""`` only narrows the type of a cell that is always text.
                source=cells["source"] or "",
                classifications={name: cells[name] for name in publisher_columns},
            )
        )
    return tuple(cofactors)


def _contributors_of(directory: Path, name: str, assembly: str) -> tuple[_Contributor, ...]:
    """Return the annotations that contributed to the one registered in ``directory``.

    The completion record is what knows, and it is asked rather than the name parsed: a
    **Merged annotation**'s name is the ``+``-join of what went into it, but splitting on
    ``+`` would guess at which component each half came from, and a caller cannot attribute
    genes to a component nobody wrote down. No ``merged_from`` marker means one
    contributor, the annotation itself, against the assembly it is registered for.
    """
    record = read_record(directory)
    merged = record.details.get(_MERGED_FROM_KEY) if record is not None else None
    if not merged:
        return (_Contributor(component=None, assembly=assembly, annotation=name),)
    return tuple(
        _Contributor(
            component=str(entry[_MERGED_COMPONENT_KEY]),
            assembly=str(entry[_MERGED_COMPONENT_KEY]),
            annotation=str(entry[_MERGED_ANNOTATION_KEY]),
        )
        for entry in merged
    )


def _declared_categories(
    contributors: Iterable[tuple[_Contributor, CuratedGeneList]],
) -> tuple[str, ...]:
    """Return every category the contributors declare between them, in declaration order.

    First-listed contributor first and each one's own order kept, with a category two of
    them declare appearing once, where the first put it — so a merged annotation reads as
    the union of its parts rather than as one part with the rest appended.
    """
    declared: dict[str, None] = {}
    for _contributor, listed in contributors:
        declared.update(dict.fromkeys(listed.categories))
    return tuple(declared)


def _gene_id_stem(gene_id: str) -> str:
    """Return ``gene_id`` with its version dropped — its **Gene id stem**.

    Everything before the first separator, and the whole id when it carries none, which is
    what makes an unversioned annotation's gene its own stem. GENCODE's pseudoautosomal
    ``ENSG00000182378.14_PAR_Y`` stems to ``ENSG00000182378`` alongside the copy it is of,
    which is the collision the caller is handed both halves of.
    """
    stem, separator, _version = gene_id.partition(_VERSION_SEPARATOR)
    return stem if separator else gene_id


def _gene_ids_by_stem(database: Path, wanted: Container[str]) -> tuple[dict[str, list[str]], bool]:
    """Return the gene ids in ``database`` under each wanted stem, and whether it has any.

    One pass, a row at a time off a cursor over the gene features alone: the database is
    queried rather than read, so a GENCODE-sized annotation costs an index walk of its gene
    rows and holds only what matched. The ids under a stem arrive ascending because the
    query does, so two machines answer in one order.

    The second half of the answer is *whether the database holds a gene at all*, because
    **no genes** and **no matching genes** are different facts that the mapping alone
    cannot tell apart, and the caller says something different about each.
    """
    found: dict[str, list[str]] = {}
    any_genes = False
    handle = gffutils.FeatureDB(str(database))
    try:
        for row in handle.execute(_GENE_IDS_QUERY):
            any_genes = True
            gene_id = str(row["id"])
            stem = _gene_id_stem(gene_id)
            if stem in wanted:
                found.setdefault(stem, []).append(gene_id)
    finally:
        # The registry hands back an answer and not a handle, so the SQLite connection
        # this opened closes here — including on the way out of an exception.
        handle.conn.close()
    return found, any_genes


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
    """Place, build and record the GTF at ``source``, as :meth:`register_path` describes.

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
    annotation = AnnotationRegistry.locate(assembly, cache_dir).register(
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

    :meth:`AnnotationRegistry.register_path` addressed by assembly name, and answering
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
    annotation = AnnotationRegistry.locate(assembly, cache_dir).register_path(
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
    """Return the :class:`~genome.io.results.RegisteredAnnotation` a just-finished run left.

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
