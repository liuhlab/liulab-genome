"""Where an assembly's FASTA comes from — resolved before a build, and checked after one.

The **Source** as a value. Naming an assembly is the whole interface, so one name has to
answer *where do these bytes come from*, and there are three answers: a URL — the one the
assembly's **Assembly metadata** row pins, else the UCSC golden path derived from the name
— the local path or URL a caller handed over, or a recipe, *these components* (ADR-0008).
:func:`resolve_source` gives that answer once and the registration dispatches on which kind
came back, rather than each step asking the question again in its own words.

Deliberately the bottom of the stack: nothing here fetches, writes, or opens a
:class:`~genome.genome.Genome`. It reads the **Completion marker** already on disk, the
curated table, and the name — which is what lets :mod:`genome.io.download` import this at
the top of the file. The module that *builds* a chimera cannot be imported there, because
:mod:`genome.io.gtf` reaches back into the downloader for the package's one fetch step.

What a finished build wrote down about its own source lives here too. :class:`ChimeraDetails`
is the ``details`` shape a chimera's completion record carries — written by
:class:`~genome.io.chimera.ChimeraBuilder` and read by everything else, so the keys are
spelled in one module and nowhere else — and :func:`components_status` asks the one thing
that record cannot answer about itself: are the components it names still the ones whose
bytes it holds?

Examples
--------
>>> from genome.io.registration import AssemblyDir
>>> nowhere = AssemblyDir.locate("hg38", "/tmp/definitely-not-a-build")
>>> resolve_source(nowhere, metadata=None, golden_path_url="https://example.org/hg38.fa.gz")
FetchedSource(url='https://example.org/hg38.fa.gz', derived=True)
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genome.chimera import ChimeraNamingError, derive_name, split_name
from genome.io.completion import (
    CompletionRecord,
    RegistrationMismatchError,
    read_record,
    record_path,
)
from genome.io.registration import AssemblyDir, assembly_data_dir, assembly_repair_command
from genome.metadata import AssemblyMetadata, lookup_assembly

#: ``details`` key holding the run of underscores this chimera's chromosome names were
#: written with. Recorded rather than assumed: a component carrying a doubled underscore
#: of its own forces a longer run, and every later reader must split with the run the
#: build actually used.
_SEPARATOR_KEY = "separator"

#: ``details`` key holding one entry per component, in the sorted order the chimera's own
#: name spells them. A list of objects rather than a mapping, so an entry stays one
#: self-describing record as later builds add to it.
_COMPONENTS_KEY = "components"

#: Keys of one entry under :data:`_COMPONENTS_KEY`: the component assembly name, and the
#: sha256 that component's *own* completion record pinned when this chimera was built —
#: a record-to-record fact, with no bytes rehashed.
_COMPONENT_NAME_KEY = "name"
_COMPONENT_DIGEST_KEY = "sha256"

#: …and the two beside them that describe the **Merged annotation**: which of that
#: component's annotations went into it, and what that annotation's own record pinned.
#: ``null`` is a component that contributed nothing to a merge that did happen. A build
#: that registered no merged annotation at all omits both keys instead, exactly as
#: ``tool_versions`` omits a tool that never answered — an absent key is a fact nobody
#: gathered, which is not the same as a fact whose answer is *none*.
_COMPONENT_ANNOTATION_KEY = "annotation"
_COMPONENT_ANNOTATION_DIGEST_KEY = "annotation_sha256"

#: What joins the contributing annotations' names into the merged one's. Deliberately not
#: the ``_`` that joins component names into a chimera's own: one name, read left to
#: right, then says which level each join is at.
_ANNOTATION_JOIN = "+"

#: Every component of this chimera was compared against its own record and is still what
#: it was. The strong answer :func:`components_status` can give.
COMPONENTS_UNCHANGED = "unchanged"

#: A comparison could not be made — a digest is absent on one side or the other — so this
#: chimera is *unproven* rather than proven stale. Reported instead of left silent,
#: because a caller who cannot tell the two apart reads the weaker one as the stronger.
COMPONENTS_UNKNOWN = "unknown"


@dataclass(frozen=True)
class FetchedSource:
    """The assembly's FASTA is downloaded from a URL.

    One kind and not two, because a pinned URL and a derived one differ in exactly one
    consequence and it is carried here: validation is a property of the source (ADR-0003),
    so a URL the assembly's row named is fetched as it stands, while one derived from the
    assembly name has that name checked at UCSC first — a typo then fails on a ``HEAD``
    request rather than three minutes into a download.

    Attributes
    ----------
    url : str
        The URL the FASTA is fetched from.
    derived : bool
        Whether :attr:`url` was derived from the assembly name because no row pinned one.

    Examples
    --------
    >>> FetchedSource("https://example.org/hg38.fa.gz", derived=False).derived
    False
    """

    url: str
    derived: bool


@dataclass(frozen=True)
class SeededSource:
    """The caller handed the bytes over — a local path or a URL of their own.

    Never resolved, only declared: this is what ``Genome(path_or_url=...)`` and ``genome
    register <name> --source`` say, so nothing is inferred and UCSC is never consulted.
    The bytes are recorded rather than compared, and the assembly name degrades to a label
    for the directory they land in (ADR-0003).

    Attributes
    ----------
    location : str or pathlib.Path
        The local FASTA path or the http(s)/ftp/sftp URL that was named.

    Examples
    --------
    >>> SeededSource("/data/my_ref.fa").location
    '/data/my_ref.fa'
    """

    location: str | Path


@dataclass(frozen=True)
class ComponentSource:
    """A recipe rather than a location: this assembly's FASTA is *these components*.

    The third kind, and what makes a **Chimera** an assembly rather than a type of its own
    (ADR-0008). Nothing is fetched — the bytes are copied from components already prepared
    on this disk — so naming one is never a way to ask for its parts.

    Attributes
    ----------
    components : tuple of str
        The **Component** assembly names, in the sorted order the chimera's name spells
        them.

    Examples
    --------
    >>> ComponentSource(("ce11", "ecHT115")).components
    ('ce11', 'ecHT115')
    """

    components: tuple[str, ...]


#: An assembly's **Source**, in the three kinds one name can resolve to. The registration
#: dispatches on which of them came back from :func:`resolve_source`; a
#: :class:`SeededSource` reaches it from the caller instead, since a source somebody named
#: is not one anything has to work out.
Source = FetchedSource | SeededSource | ComponentSource


def fetched_source(metadata: AssemblyMetadata | None, golden_path_url: str) -> FetchedSource:
    """Return the URL this assembly's FASTA is downloaded from, and how that was decided.

    The metadata row's pinned source when it has one, and ``golden_path_url`` when it does
    not — the whole of *where does an ordinary assembly come from*, in one place, so that
    the URL a fetch uses and the question of whether the assembly name still has to be
    confirmed at UCSC cannot drift apart.

    Parameters
    ----------
    metadata : genome.metadata.AssemblyMetadata or None
        The assembly's curated row, or ``None`` for one the table does not list.
    golden_path_url : str
        The URL derived from the assembly name, used when nothing is pinned.

    Returns
    -------
    FetchedSource
        The URL to fetch, and whether it was derived rather than pinned.

    Examples
    --------
    >>> fetched_source(None, "https://example.org/hg38.fa.gz")
    FetchedSource(url='https://example.org/hg38.fa.gz', derived=True)
    """
    pinned = metadata.source_url if metadata is not None else None
    return FetchedSource(url=pinned if pinned else golden_path_url, derived=pinned is None)


def is_prepared(assembly: str) -> bool:
    """Return whether ``assembly`` is registered under the shared data root, by its record alone.

    By name and not by path, exactly as a chimera's recorded components are found again: a
    component is addressed by the key it was registered under. A directory holding files
    but no record is *not* prepared — nothing vouches for it.

    Parameters
    ----------
    assembly : str
        The assembly name to look for.

    Returns
    -------
    bool
        Whether a completion record is there.

    Examples
    --------
    >>> is_prepared("definitely-not-an-assembly")
    False
    """
    return read_record(assembly_data_dir(assembly)) is not None


def _could_be_a_component(assembly: str) -> bool:
    """Whether ``assembly`` names something a chimera could be built from.

    Prepared here, or listed in the curated table. The second clause is what separates a
    chimera's name from a free-form local key on a machine that holds neither: the shipped
    row for a chimera's components is why ``ce11_ecHT115`` reads as two assemblies while
    ``my_ref`` reads as one name somebody chose (ADR-0003).
    """
    return is_prepared(assembly) or lookup_assembly(assembly) is not None


def resolve_source(
    assembly_dir: AssemblyDir,
    *,
    metadata: AssemblyMetadata | None,
    golden_path_url: str,
) -> FetchedSource | ComponentSource:
    """Return where the assembly in ``assembly_dir`` gets its FASTA — a URL, or components.

    What a name means, decided once. Both ways in — :class:`~genome.genome.Genome` and
    :func:`~genome.io.download.register_assembly` — reach this through one registration
    step, which is what makes ``genome register <name>`` one command for every kind of
    **Source** and keeps the command line a thin client, since neither the resolution nor
    its refusals are written there. Four checks, in this order:

    1. **A record here.** An existing completion record already says what this assembly is,
       and is believed outright: a chimera's record rebuilds a chimera, and any other
       record — a plain ``hg38_mm10`` somebody seeded years ago — keeps the assembly
       whatever it was registered as. Only a record that is *lost* falls through to the
       name.
    2. **A source the caller named.** An explicitly seeded assembly is what the caller said
       it is, so it never reaches this function: ``--source`` is a :class:`SeededSource`
       handed straight to the registration. That route consults the record first as well —
       a registration already here is returned from it without the source being read again
       — so the first check stands whichever way in was taken, and the source answers only
       where there is nothing to believe.
    3. **The name.** It splits on ``_`` into two or more parts, and every one of them is
       either prepared here or listed in the curated table. That second clause is the whole
       separation between ``ce11_ecHT115`` and a free-form local key like ``my_ref`` on a
       machine holding neither (ADR-0003, ADR-0008).
    4. **Today's path**, for everything else — the FASTA is fetched from the pinned source
       or the derived golden-path URL, exactly as it always was.

    Rebuilding from scratch is not one of the checks and never overrules them: what an
    assembly *is* is not something a repair may change.

    Parameters
    ----------
    assembly_dir : genome.io.registration.AssemblyDir
        The **Assembly dir** in question. It carries the name being resolved and the
        directory whose record answers first.
    metadata : genome.metadata.AssemblyMetadata or None
        The assembly's curated row, or ``None`` for one the table does not list.
    golden_path_url : str
        The URL derived from the assembly name, used when nothing pins one.

    Returns
    -------
    FetchedSource or ComponentSource
        The components this assembly is concatenated from, or the URL its FASTA is
        downloaded from.

    Raises
    ------
    FileNotFoundError
        If the name spells a legal chimera's components in the wrong order. The message
        names the canonical spelling — a mis-ordered name being detectable is what leaves a
        ``--component`` flag with nothing left to buy.

    Examples
    --------
    >>> from genome.io.registration import AssemblyDir
    >>> nowhere = AssemblyDir.locate("my_ref", "/tmp/definitely-not-a-build")
    >>> resolve_source(nowhere, metadata=None, golden_path_url="https://example.org/x.fa.gz")
    FetchedSource(url='https://example.org/x.fa.gz', derived=True)
    """
    record = assembly_dir.read_record()
    if record is not None:
        details = ChimeraDetails.from_record(record)
        if details is None:
            return fetched_source(metadata, golden_path_url)
        return ComponentSource(tuple(details.components))
    try:
        candidates = split_name(assembly_dir.assembly)
    except ChimeraNamingError:
        return fetched_source(metadata, golden_path_url)
    if not all(_could_be_a_component(name) for name in candidates):
        return fetched_source(metadata, golden_path_url)
    canonical = derive_name(candidates)
    if canonical != assembly_dir.assembly:
        raise FileNotFoundError(
            f"nothing is registered as {assembly_dir.assembly!r} in {assembly_dir.path}, and "
            f"it is not how a chimera of {', '.join(candidates)} is spelled: a chimera's name "
            f"is its components sorted, so that one set of components means one "
            f"directory whatever order they are typed in (ADR-0008). Build it with "
            f"`genome register {canonical}`."
        )
    return ComponentSource(candidates)


def merged_annotation_name(annotations: Sequence[str]) -> str:
    """Return the **Registered name** a merge of ``annotations`` is filed under.

    The contributing annotations' names joined by ``+``, in the sorted-component order the
    chimera's own name spells — ``wormbase_ws298+refseq_rs_2025_06_26``. Derived, like
    everything else about a chimera, so it changes the moment any component's default
    annotation does and a database built from the old set can never be found under it.
    Written by the build and read back from the record with this one spelling, so the two
    cannot disagree.

    It needs no parse-back: what a merged annotation is made of is recovered from the
    components, and written down in its own record besides. And it is not asked to carry
    *which* components contributed — a chimera with a component that contributes nothing
    spells the same name a different subset would — which is why
    :func:`~genome.io.gtf.register_merged_gtf` adopts nothing from disk and writes the
    annotation every time it runs.

    Parameters
    ----------
    annotations : sequence of str
        The contributing annotations' registered names, in sorted-component order.

    Returns
    -------
    str
        The merged annotation's registered name.

    Examples
    --------
    >>> merged_annotation_name(["wormbase_ws298", "refseq_rs_2025_06_26"])
    'wormbase_ws298+refseq_rs_2025_06_26'
    """
    return _ANNOTATION_JOIN.join(annotations)


@dataclass(frozen=True)
class ComponentDetails:
    """What a chimera's completion record says about one of its components.

    One entry of the ``details`` shape :class:`~genome.io.chimera.ChimeraBuilder` writes,
    read back. Every field is a fact about the component *at the time this chimera was
    built*, taken from that component's own records rather than by rehashing anything —
    which is what lets a later pass notice a component re-registered underneath the chimera.

    Attributes
    ----------
    name : str
        The **Component** assembly name.
    sha256 : str or None
        The digest that component's completion record pinned, or ``None`` when it pinned
        none — *unknown*, rather than wrong.
    annotation : str or None
        The **Registered name** of the annotation it contributed to the **Merged
        annotation**, or ``None`` when it contributed none.
    annotation_sha256 : str or None
        That annotation's own recorded digest, or ``None``.

    Examples
    --------
    >>> ComponentDetails("ce11", "1a2b3c", "wormbase_ws298", "4d5e6f").annotation
    'wormbase_ws298'
    """

    name: str
    sha256: str | None
    annotation: str | None
    annotation_sha256: str | None

    def as_entry(self, *, merged: bool) -> dict[str, Any]:
        """Return this component's entry as a record writes it down.

        The inverse of the reading :meth:`ChimeraDetails.from_details` does, so the keys
        are spelled once. The two annotation keys are written only when a merged annotation
        was registered: a build that registered none says nothing about annotations at all
        rather than writing ``null`` beside every component.

        Parameters
        ----------
        merged : bool
            Whether this build registered a **Merged annotation**.

        Returns
        -------
        dict
            One entry of a completion record's ``details``, ready to serialize.

        Examples
        --------
        >>> ComponentDetails("ce11", "1a2b3c", None, None).as_entry(merged=False)
        {'name': 'ce11', 'sha256': '1a2b3c'}
        """
        entry: dict[str, Any] = {
            _COMPONENT_NAME_KEY: self.name,
            _COMPONENT_DIGEST_KEY: self.sha256,
        }
        if merged:
            entry[_COMPONENT_ANNOTATION_KEY] = self.annotation
            entry[_COMPONENT_ANNOTATION_DIGEST_KEY] = self.annotation_sha256
        return entry


@dataclass(frozen=True)
class ChimeraDetails:
    """What a chimera's completion record says about the build that produced it.

    The **Source** a finished chimera records, in the ``details`` shape
    :class:`~genome.io.chimera.ChimeraBuilder` writes — written and read here, so nothing
    else has to know its keys. It answers the only question that decides whether an
    assembly is a **Chimera** at runtime — the record, never the metadata row — and
    carries the facts a later pass needs: the separator its chromosome names were written
    with, and what each component and each contributed annotation was at build time.

    Attributes
    ----------
    separator : str
        The run of underscores this chimera's chromosome names carry.
    component_details : tuple of ComponentDetails
        One entry per **Component**, in the sorted order the chimera's name spells them.

    Examples
    --------
    >>> details = ChimeraDetails(
    ...     "__",
    ...     (
    ...         ComponentDetails("ce11", "1a2b3c", "wormbase_ws298", "4d5e6f"),
    ...         ComponentDetails("ecHT115", "7a8b9c", None, None),
    ...     ),
    ... )
    >>> details.components
    ['ce11', 'ecHT115']
    """

    separator: str
    component_details: tuple[ComponentDetails, ...]

    @property
    def components(self) -> list[str]:
        """The component assembly names, sorted — a fresh list each call."""
        return [entry.name for entry in self.component_details]

    @property
    def merged_annotation(self) -> str | None:
        """The **Registered name** of the **Merged annotation** this build wrote, or ``None``.

        Read back rather than looked up: the name is the contributing annotations' names
        joined by :func:`merged_annotation_name`, which is what spelled it when the merge
        was registered. ``None`` when no component contributed one, which is a build that
        registered no annotation at all rather than an empty one.

        Examples
        --------
        >>> ChimeraDetails(
        ...     "__",
        ...     (
        ...         ComponentDetails("ce11", "1a2b3c", "wormbase_ws298", "4d5e6f"),
        ...         ComponentDetails("ecHT115", "7a8b9c", None, None),
        ...     ),
        ... ).merged_annotation
        'wormbase_ws298'
        """
        contributed = [
            entry.annotation for entry in self.component_details if entry.annotation is not None
        ]
        return merged_annotation_name(contributed) if contributed else None

    def as_details(self, *, merged: bool) -> dict[str, Any]:
        """Return these details as a completion record writes them down.

        The inverse of :meth:`from_details`. Kept to facts a later pass cannot re-derive
        from the name alone: which spelling was used, what each component was when this
        chimera was built, and — when there is a merged annotation — which of each
        component's annotations went into it.

        Parameters
        ----------
        merged : bool
            Whether this build registered a **Merged annotation**.

        Returns
        -------
        dict
            The ``details`` of a chimera's completion record, ready to serialize.

        Examples
        --------
        >>> ChimeraDetails("__", (ComponentDetails("ce11", None, None, None),)).as_details(
        ...     merged=False
        ... )
        {'separator': '__', 'components': [{'name': 'ce11', 'sha256': None}]}
        """
        return {
            _SEPARATOR_KEY: self.separator,
            _COMPONENTS_KEY: [entry.as_entry(merged=merged) for entry in self.component_details],
        }

    @classmethod
    def from_record(cls, record: CompletionRecord | None) -> ChimeraDetails | None:
        """Read a completion record's chimera details, or ``None`` when it has none.

        :meth:`from_details` over the record's ``details``, and ``None`` for an absent
        record — the shape a caller holding a record rather than a payload wants.

        Parameters
        ----------
        record : genome.io.completion.CompletionRecord or None
            The record to read, as :func:`~genome.io.completion.read_record` returns it.

        Returns
        -------
        ChimeraDetails or None
            The details, or ``None`` when the record is not a chimera's.

        Examples
        --------
        >>> ChimeraDetails.from_record(None) is None
        True
        """
        return None if record is None else cls.from_details(record.details)

    @classmethod
    def from_details(cls, details: Mapping[str, Any]) -> ChimeraDetails | None:
        """Read a completion record's ``details``, or ``None`` when they are not a chimera's.

        ``None`` means *this is not a chimera* — the details of an ordinary downloaded or
        seeded assembly, or ones that do not carry the shape a chimera build writes. Those
        read alike on purpose: nothing but a build of this package's own writing may make
        an assembly answer as a chimera.

        The two annotation fields are optional, and a build that registered no merged
        annotation writes neither; both then read as ``None``.

        Taking the mapping rather than the record is what lets a caller that already holds
        one — the CLI, whose ``register`` payload *is* the record — answer from what it has
        instead of reading the same file again.

        Parameters
        ----------
        details : mapping of str to object
            A registration record's ``details``. Anything else it holds is ignored, and an
            empty mapping reads as *not a chimera*.

        Returns
        -------
        ChimeraDetails or None
            The details, or ``None`` when they are not a chimera build's.

        Examples
        --------
        >>> ChimeraDetails.from_details({}) is None
        True
        """
        separator = details.get(_SEPARATOR_KEY)
        entries = details.get(_COMPONENTS_KEY)
        if not isinstance(separator, str) or not isinstance(entries, list):
            return None
        components: list[ComponentDetails] = []
        for entry in entries:
            if not isinstance(entry, dict) or not isinstance(entry.get(_COMPONENT_NAME_KEY), str):
                return None
            components.append(
                ComponentDetails(
                    name=entry[_COMPONENT_NAME_KEY],
                    sha256=_text(entry.get(_COMPONENT_DIGEST_KEY)),
                    annotation=_text(entry.get(_COMPONENT_ANNOTATION_KEY)),
                    annotation_sha256=_text(entry.get(_COMPONENT_ANNOTATION_DIGEST_KEY)),
                )
            )
        return cls(separator=separator, component_details=tuple(components))


def read_chimera_details(directory: Path) -> ChimeraDetails | None:
    """Return the chimera details recorded in ``directory``, or ``None`` when it has none.

    :meth:`ChimeraDetails.from_record` over :func:`~genome.io.completion.read_record` —
    the one call that answers *is the assembly registered here a chimera?*, and the only
    thing :attr:`genome.genome.Genome.components` consults.

    Parameters
    ----------
    directory : pathlib.Path
        An **Assembly dir** — the directory a registration filled.

    Returns
    -------
    ChimeraDetails or None
        The details, or ``None`` for an assembly that is not a chimera (or is not
        registered at all).

    Examples
    --------
    >>> from pathlib import Path
    >>> read_chimera_details(Path("/tmp/definitely-not-a-build")) is None
    True
    """
    return ChimeraDetails.from_record(read_record(directory))


def components_status(assembly_dir: AssemblyDir) -> str | None:
    """Check every component of the assembly in ``assembly_dir``, and report the answer.

    The one failure a digest of a chimera's own bytes cannot show. Those bytes are a copy
    of its components', so a component registered again underneath it leaves the chimera
    intact, agreeing with its own record, and no longer a copy of anything that exists —
    silently stale sequence, and stale gene models one level down. Both are caught here,
    and both are **record against record**: the digests this chimera wrote down are
    compared against the ones the components' own records pin now, so this reads a
    handful of small JSON files and not one base of sequence.

    One entry point for two needs, because they are one comparison. A component proved to
    have changed **raises**, wherever the question was asked from — reopening a finished
    chimera, rebuilding one, verifying one. What is *returned* is what nothing raised over,
    and a caller that only wants the refusal ignores it: a comparison that could not be
    made is :data:`COMPONENTS_UNKNOWN` and is not a pass, and a surface silent about it
    would say exactly what it says when everything agreed.

    An assembly with no components recorded has nothing to compare and returns at once,
    which is what makes an ordinary assembly pay nothing rather than be asked about.
    Likewise an absent digest on either side, which means *unknown* rather than wrong: a
    component that pinned none, or an annotation registered before its digest was
    recorded, leaves that component unguarded rather than refused.

    Parameters
    ----------
    assembly_dir : genome.io.registration.AssemblyDir
        The **Assembly dir** the chimera was built in. It carries the assembly's name,
        quoted in the error along with the command that repairs it, and it is what each
        component is found beside.

    Returns
    -------
    str or None
        :data:`COMPONENTS_UNCHANGED` when every component was compared and agreed,
        :data:`COMPONENTS_UNKNOWN` when any comparison could not be made, and ``None`` for
        an assembly that is not a chimera — which has no components to be asked about.

    Raises
    ------
    genome.io.completion.RegistrationMismatchError
        If a component's FASTA, or the annotation it contributed to the **Merged
        annotation**, is not the one this chimera was built from. The message names both
        digests and ``genome register <assembly> --force``, which rebuilds both halves.

    Examples
    --------
    >>> from genome.io.registration import AssemblyDir
    >>> nowhere = AssemblyDir.locate("notAChimera", "/tmp/definitely-not-a-build")
    >>> components_status(nowhere) is None
    True
    """
    details = read_chimera_details(assembly_dir.path)
    if details is None:
        return None
    compared = [
        _component_was_compared(entry, chimera=assembly_dir) for entry in details.component_details
    ]
    return COMPONENTS_UNCHANGED if all(compared) else COMPONENTS_UNKNOWN


def _component_was_compared(entry: ComponentDetails, *, chimera: AssemblyDir) -> bool:
    """Return whether one component was compared at all — raising unless it is still itself.

    Two digests are looked at, the component's FASTA and then the annotation it
    contributed, and one absent on either side leaves nothing compared: the answer is
    ``False``, which is how a caller reporting the outcome says *unknown* rather than
    claiming a pass nobody checked. A comparison that was made and disagreed raises.

    The component is found **beside the chimera**, by the name its record carries: a
    chimera's record holds component names and no paths, which is what keeps a registered
    directory movable, so the name is resolved against where the chimera itself is rather
    than against whatever **Data dir** this process is pointed at. Under the ordinary
    layout those are the same directory.
    """
    component = chimera.sibling(entry.name)
    component_dir = component.path
    assembly = chimera.assembly
    directory = chimera.path
    record = read_record(component_dir)
    current = None if record is None else record.sha256
    if _disagree(entry.sha256, current):
        raise RegistrationMismatchError(
            f"the chimera {assembly} in {directory} was built from a different "
            f"{entry.name} than the one registered now: it recorded {entry.sha256} for "
            f"that component, and {record_path(component_dir)} now pins {current}. "
            f"{entry.name} was registered again afterwards, so the sequences it "
            f"contributed here are a copy of bytes that are no longer anywhere — which "
            f"nothing about this chimera's own files can show, since they are unchanged. "
            f"Build it again with `{assembly_repair_command(assembly)}`."
        )
    known = _both_known(entry.sha256, current)
    if entry.annotation is None:
        return known
    gtf_dir = component.annotation_dir(entry.annotation)
    annotation = read_record(gtf_dir)
    current_annotation = None if annotation is None else annotation.sha256
    if _disagree(entry.annotation_sha256, current_annotation):
        raise RegistrationMismatchError(
            f"the merged annotation of the chimera {assembly} in {directory} was merged "
            f"from a different {entry.annotation} of {entry.name}: it recorded "
            f"{entry.annotation_sha256} for that annotation, and {record_path(gtf_dir)} "
            f"now pins {current_annotation}. That annotation was registered again after "
            f"the merge, so the merged gene models are not the ones {entry.name} carries "
            f"— the failure the sequence digests cannot show, since no base of this "
            f"chimera's FASTA changed. Build it again with "
            f"`{assembly_repair_command(assembly)}`, which rewrites the annotation and "
            f"the FASTA together."
        )
    return known and _both_known(entry.annotation_sha256, current_annotation)


def _disagree(recorded: str | None, current: str | None) -> bool:
    """Whether two recorded digests are known to differ — unknown on either side is not."""
    return recorded is not None and current is not None and recorded != current


def _both_known(recorded: str | None, current: str | None) -> bool:
    """Whether two recorded digests were both there to be compared."""
    return recorded is not None and current is not None


def _text(value: Any) -> str | None:
    """Return ``value`` when it is a string, else ``None`` — a record field read loosely.

    A record is JSON somebody else's version may have written, so a field that is absent,
    null or the wrong type all read as *not known*, and none of them raises.
    """
    return value if isinstance(value, str) else None
