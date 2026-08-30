"""One **Assembly**'s annotations — the four states each is in, and what may be asked.

:class:`AnnotationRegistry` is here, and it is the way in to everything the
:mod:`annotation <genome.annotation>` package does. Bound once to one assembly — its
name, its **Assembly dir** and its ``chrom.sizes`` — it holds every annotation that
assembly has and answers everything about them: which are registered, which are broken,
which the table offers, which is the **Default annotation**, where one's GTF is, and the
two acts that add one. Everything that needs the four-way state asks a registry rather
than assembling it again: a :class:`~genome.assembly.genome.Genome` holds one for its lifetime, and
each assembly-addressed function here builds one for the length of the call.

**Every annotation directory is registered, broken, or not begun**, and the middle one has
its own listing: :func:`list_broken_annotations`. Registering an annotation raises over a
directory it cannot trust, but listing must not — one annotation nobody can vouch for
cannot be allowed to stop a **Genome** opening or hide the annotations beside it — so the
two lists are reported side by side and each broken one carries the command that repairs
it.

**What the lab offers and what this machine holds are different questions.** The first is
the annotation table's to answer (:func:`~genome.annotation.metadata.list_annotation_metadata`), the
second this disk's (:func:`list_annotations`); :meth:`AnnotationRegistry.status` sets one
against the other, and :func:`default_annotation` is the one rule that picks the **Default
annotation** out of both. Those three scans plus that rule are what a registry is: they are
read together, once, at construction, and every later question is answered from the answer.

**An annotation can also say which genes are in a category.**
:meth:`AnnotationRegistry.gene_list` and :meth:`AnnotationRegistry.gene_lists` answer from
the **Curated gene list** :mod:`genome.annotation.curated` ships for that annotation, addressed by
the same **Registered name** everything else here is. A **Merged annotation** answers per
contributor, read off its record, so its genes stay attributable to the component they came
from. Neither ever answers with an empty collection: an annotation that ships no list and
one whose list does not declare the category asked for raise errors of their own, since a
caller acts differently on those two facts and a silent zero is what that surface exists to
prevent.

The class stays **one class with one interface**, and its methods straddle three of the
package's four modules: the query surface is here, :meth:`AnnotationRegistry.register` and
:meth:`AnnotationRegistry.register_path` are placement
(:mod:`genome.annotation.registration`), and
:meth:`AnnotationRegistry.resolve_gene_ids` is stem resolution
(:mod:`genome.annotation.stems`). The two that write cannot leave — each folds what it
wrote back into the four states so a later read is current without touching the disk again,
and each is held to the ``chrom.sizes`` the registry was bound to — and neither can the
third, being a name resolution, one call and the shape of an answer. So the other modules
hold free functions, dataclasses and errors, and the class calls across them.

:func:`annotation_status`, :func:`gene_list` and :func:`gene_lists` are the by-assembly-name
forms the CLI calls, each building a registry for the length of the call.

Examples
--------
>>> from genome.annotation.registry import AnnotationRegistry
>>> AnnotationRegistry.locate("sacCer3", "/tmp/definitely-not-an-assembly").registered
[]
"""

from __future__ import annotations

import shlex
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from genome.annotation.curated import (
    CuratedGeneList,
    GeneCategoryNotDeclaredError,
    NoGeneCategoriesError,
    curated_annotations,
    curated_gene_list,
)
from genome.annotation.metadata import (
    AnnotationMetadata,
    list_annotation_metadata,
    lookup_annotation,
)
from genome.annotation.registration import (
    _MERGED_ANNOTATION_KEY,
    _MERGED_COMPONENT_KEY,
    _MERGED_FROM_KEY,
    GtfAnnotation,
    _already_registered,
    _annotation_files,
    _annotation_repair,
    _annotations_root,
    _assembly_chromosomes,
    _build_and_record,
    _chromosome_check_details,
    _elide,
    _fetch_gtf,
    _path_repair_command,
    _proven_gtf,
    _register_gtf,
    _register_gtf_command,
    _reject_unknown_chromosomes,
    _repair_command,
    annotation_register_command,
)
from genome.annotation.stems import NoGeneFeaturesError, ResolvedGeneIds, gene_ids_by_stem
from genome.assembly.registration import AssemblyDir
from genome.store.completion import (
    RegistrationError,
    check_registration,
    disagreements,
    read_record,
)


class AnnotationNotRegisteredError(KeyError):
    """No annotation of that name is registered here, so there is no path to hand back.

    Routinely *not* a mistake. An assembly's **Default annotation** comes from the
    curated table, and on a fresh machine the table's choice is exactly what nobody has
    registered yet — so a :class:`~genome.assembly.genome.Genome` opens with that default named
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
    genome.annotation.registry.AnnotationNotRegisteredError: "no annotation ...
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
        :func:`~genome.store.completion.check_registration`'s own message, so re-registering
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
    ...     repair="genome annotation register-gtf hg38 /tmp/mine.gtf mine --force",
    ... )
    >>> broken.repair
    'genome annotation register-gtf hg38 /tmp/mine.gtf mine --force'
    """

    name: str
    directory: Path
    problem: str
    repair: str


@dataclass(frozen=True)
class AnnotationStatusRow:
    """One annotation, in whichever of its states it is: offered, registered, broken.

    One shape for all of them, so a reader never has to ask which fields a row has — a
    name the table does not list carries the table's columns as ``None``, and one nothing
    is wrong with carries the broken columns as ``None``. :attr:`registered` and
    :attr:`broken` are never both true: a registration nothing vouches for is not one, and
    :attr:`state` is that invariant said in one word.

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
    >>> row.state
    'offered, not registered'
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

    @property
    def state(self) -> str:
        """Which of its states this row is in, in the words a surface prints.

        ``broken`` first, because the three fields it is read from are not independent:
        a broken annotation is not registered — no record vouches for it — so answering
        with the absence of one would be true and useless, and it is the state that needs
        acting on. A row nothing offers is one this disk holds and the table does not, so
        the only thing left to say about it is that.

        Returns
        -------
        str
            One of ``broken``, ``registered, not offered``, ``registered``, or
            ``offered, not registered``.
        """
        if self.broken:
            return "broken"
        if not self.offered:
            return "registered, not offered"
        return "registered" if self.registered else "offered, not registered"

    def as_json(self) -> dict[str, Any]:
        """Return this row as ``--json`` serializes it: every attribute above, in order.

        :attr:`state` is not among them: it is read from :attr:`broken`, :attr:`offered`
        and :attr:`registered`, which are all here, so writing it out too would be a
        second spelling of the same rule for a reader to disagree with.

        Returns
        -------
        dict
            The row's fields, under their own names.
        """
        return asdict(self)


@dataclass(frozen=True)
class AnnotationStatus:
    """What one assembly's table offers, set against what is registered on this machine.

    :meth:`AnnotationRegistry.status`'s answer, and what ``genome annotation list`` prints. Two
    questions joined for one reader, with a third riding along because this is where anyone
    would look for it: a directory that cannot be trusted is ``broken`` rather than
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
    >>> status.default_summary
    'default: (none)'
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

    @property
    def default_summary(self) -> str:
        """The closing sentence naming the **Default annotation**, and how to get it.

        Four answers, and three of them tell the reader what to do next: nothing decided a
        default; one is decided and registered, which needs no advice; one is decided and
        broken, which is repaired by the command its own row already carries; or one is
        decided and absent, which is the ordinary state of a fresh machine and is
        registered by :func:`annotation_register_command`. Both commands come off an
        interface rather than being assembled here, so the two halves of this sentence
        cannot drift apart.

        Returns
        -------
        str
            One line, beginning ``default: `` — the whole of what ``genome annotation list``
            prints last.
        """
        default = self.default_annotation
        if default is None:
            return "default: (none)"
        row = self.default_row
        if row is not None and row.broken:
            return f"default: {default} — broken here; repair it with `{row.repair}`"
        if row is not None and row.registered:
            return f"default: {default}"
        return (
            f"default: {default} — not registered here; register it with "
            f"`{annotation_register_command(self.assembly, default)}`"
        )

    def as_json(self) -> dict[str, Any]:
        """Return this report as ``--json`` serializes it.

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


@dataclass(frozen=True)
class GeneListSource:
    """One **Curated gene list** that contributed to an answer, and what it contributed.

    What makes a **Merged annotation**'s genes attributable: one of these per contributing
    annotation, so a caller counting worm ribosomal RNA can drop the *E. coli* entry
    rather than being handed one number it cannot take apart. An annotation that is not a
    merge has exactly one, whose ``component`` is ``None``.

    ``description`` and ``source`` travel with the ids rather than being looked up
    separately, because they are what says whether these ids mean what the caller's metric
    needs: two annotations spelling a category the same way need not have curated it the
    same way.

    Attributes
    ----------
    component : str or None
        The **Component** assembly whose genes these are, for a contributor to a **Merged
        annotation**; ``None`` for anything else.
    annotation : str
        The **Registered name** of the contributing annotation.
    description : str
        What membership in this category means for that annotation.
    source : str
        Where that membership came from, and the caveats on using it.
    gene_ids : tuple of str
        The gene ids it contributed, in the order its curated list lists them.

    Examples
    --------
    >>> source = GeneListSource(
    ...     component="ce11",
    ...     annotation="wormbase_ws298",
    ...     description="the mature ribosomal RNA genes",
    ...     source="WormBase WS298 gene_biotype",
    ...     gene_ids=("WBGene00004512", "WBGene00004513"),
    ... )
    >>> source.as_json()["component"]
    'ce11'
    >>> len(source.gene_ids)
    2
    """

    component: str | None
    annotation: str
    description: str
    source: str
    gene_ids: tuple[str, ...]

    def as_json(self) -> dict[str, Any]:
        """Return this contribution as ``--json`` serializes it: every attribute, in order.

        Returns
        -------
        dict
            The fields above, under their own names, with ``gene_ids`` as a list.
        """
        return asdict(self)


@dataclass(frozen=True)
class GeneList:
    """The genes one annotation puts in one **Gene category**, attributed to their sources.

    :meth:`AnnotationRegistry.gene_list`'s answer — what ``genome annotation gene-list``
    prints and what its ``--json`` serializes. There is no empty one: an annotation that declares no
    categories, and one that declares categories but not this one, each raise a
    :class:`LookupError` of their own rather than answering with nothing, so holding one of
    these means the category was really declared and really has genes in it.

    Attributes
    ----------
    assembly : str
        The **Assembly** asked about.
    annotation : str
        The **Registered name** asked about — the merged name for a **Merged annotation**,
        whose contributors are named in :attr:`sources`.
    category : str
        The **Gene category**, as the curated lists spell it.
    sources : tuple of GeneListSource
        One entry per contributing **Curated gene list**, in contributor order. Never
        empty. A contributor that does not declare this category is simply absent — a
        bacterium has no mitochondria, and that is not a failure to report.

    Examples
    --------
    >>> genes = GeneList(
    ...     assembly="ce11",
    ...     annotation="wormbase_ws298",
    ...     category="rRNA",
    ...     sources=(
    ...         GeneListSource(None, "wormbase_ws298", "rRNA genes", "WormBase", ("a", "b")),
    ...     ),
    ... )
    >>> genes.gene_ids
    ['a', 'b']
    >>> genes.as_json()["category"]
    'rRNA'
    """

    assembly: str
    annotation: str
    category: str
    sources: tuple[GeneListSource, ...]

    @property
    def gene_ids(self) -> list[str]:
        """Every source's gene ids, concatenated in source order — a fresh list each call.

        **Concatenated and not de-duplicated.** A merge rewrites only the seqname and
        never the ``gene_id``, so two components carrying the same id would be a real
        ambiguity in the merged annotation, and collapsing it here would hide exactly that
        — a caller summing over these would silently under-count one of the two. Where
        that matters, :attr:`sources` says which contributor each id came from.
        """
        return [gene_id for source in self.sources for gene_id in source.gene_ids]

    def as_json(self) -> dict[str, Any]:
        """Return this answer as ``--json`` serializes it.

        Returns
        -------
        dict
            ``assembly``, ``annotation``, ``category``, the concatenated ``gene_ids``, and
            ``sources`` as a list of :meth:`GeneListSource.as_json` entries.
            :attr:`gene_ids` is written out beside the sources it is read from rather than
            left to the reader: assembling it is where a reader would reach for a set and
            de-duplicate, which is the one thing this answer must not do.
        """
        return {
            "assembly": self.assembly,
            "annotation": self.annotation,
            "category": self.category,
            "gene_ids": self.gene_ids,
            "sources": [source.as_json() for source in self.sources],
        }


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
    :class:`~genome.assembly.genome.Genome` being opened, and :func:`annotation_status` reporting
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
    :attr:`Genome.default_gtf_path <genome.assembly.genome.Genome.default_gtf_path>`.

    Parameters
    ----------
    offered : iterable of genome.annotation.metadata.AnnotationMetadata
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
    >>> from genome.annotation.metadata import AnnotationMetadata
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
    an :class:`~genome.assembly.registration.AssemblyDir`, so a registry cannot file an
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
    assembly_dir : genome.assembly.registration.AssemblyDir
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
        directory override. :meth:`~genome.assembly.registration.AssemblyDir.locate` is where
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
        naming ``genome annotation register <assembly> <name> --force`` (ADR-0007). That is
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
        metadata : genome.annotation.metadata.AnnotationMetadata, optional
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
        genome.store.checksum.ChecksumMismatchError
            If the row pins a sha256 and the unpacked GTF is not it; the message names both
            digests.
        genome.store.completion.UnfinishedRegistrationError
            If the annotation's directory holds files but no record.
        genome.store.completion.RegistrationMismatchError
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
        cannot be trusted raises naming ``genome annotation register-gtf <assembly> <gtf> <name>
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
        genome.store.completion.RegistrationError
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
        which of the two it is. The command behind it is ``genome annotation list``.

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
        GeneList
            The category, its gene ids, and one
            :class:`GeneListSource` per contributing curated list.

        Raises
        ------
        ValueError
            If ``name`` is omitted and no **Default annotation** is decided; the message
            names the argument that chooses one.
        AnnotationNotRegisteredError
            If nothing of that name is registered here.
        genome.annotation.curated.NoGeneCategoriesError
            If no curated list ships for that annotation — nothing can be asked of it,
            which is not the same answer as its having no genes in this category.
        genome.annotation.curated.GeneCategoryNotDeclaredError
            If it declares categories and not this one; the message lists the ones it does.
        genome.annotation.curated.GeneListAssemblyMismatchError
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
        tuple of GeneList
            One entry per declared category, in declaration order. Never empty.

        Raises
        ------
        ValueError
            If ``name`` is omitted and no **Default annotation** is decided.
        AnnotationNotRegisteredError
            If nothing of that name is registered here.
        genome.annotation.curated.NoGeneCategoriesError
            If no curated list ships for that annotation.
        genome.annotation.curated.GeneListAssemblyMismatchError
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
        :attr:`~genome.annotation.stems.ResolvedGeneIds.unresolved`, so a caller resolving
        a few thousand at once
        can see which of them this annotation does not have rather than counting the
        answer and wondering.

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
        ResolvedGeneIds
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
        found, any_genes = gene_ids_by_stem(annotation.db, frozenset(asked))
        if not any_genes:
            raise NoGeneFeaturesError(resolved_name, self.assembly)
        return ResolvedGeneIds(
            assembly=self.assembly,
            annotation=resolved_name,
            resolved={stem: tuple(found[stem]) for stem in asked if stem in found},
            unresolved=tuple(stem for stem in asked if stem not in found),
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
    what ``genome annotation list`` runs. Nothing is prepared, fetched, built or created to
    answer it — an assembly with no directory at all is the case it most needs to serve.

    Parameters
    ----------
    assembly : str
        The assembly to report on, e.g. ``"hg38"``.
    cache_dir : str or pathlib.Path, optional
        Override which assembly directory is inspected. Defaults to
        :func:`assembly_data_dir(assembly) <genome.assembly.download.assembly_data_dir>`.

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
    is what ``genome annotation gene-list`` runs. A registry is built for the length of the call, so
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
    GeneList
        The answer :meth:`AnnotationRegistry.gene_list` describes.

    Raises
    ------
    ValueError
        If ``annotation`` is omitted and no **Default annotation** is decided.
    AnnotationNotRegisteredError
        If that annotation is not registered here.
    genome.annotation.curated.NoGeneCategoriesError
        If no curated gene list ships for it.
    genome.annotation.curated.GeneCategoryNotDeclaredError
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
    annotation gene-categories`` runs, built the same way :func:`gene_list` is. Never an
    empty tuple: an annotation that declares nothing raises instead.

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
    tuple of GeneList
        One entry per declared category, in declaration order. Never empty.

    Raises
    ------
    ValueError
        If ``annotation`` is omitted and no **Default annotation** is decided.
    AnnotationNotRegisteredError
        If that annotation is not registered here.
    genome.annotation.curated.NoGeneCategoriesError
        If no curated gene list ships for it.

    Examples
    --------
    >>> [answer.category for answer in gene_lists("hg38")]    # doctest: +SKIP
    ['rRNA', 'rRNA_pseudogene', 'Mt_rRNA']
    """
    return AnnotationRegistry.locate(assembly, cache_dir).gene_lists(annotation)


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
