"""What a finished chimera wrote down about the components it is made of.

The ``details`` shape a **Chimera**'s completion record carries — written by
:class:`~genome.io.chimera.ChimeraBuilder` and read by everything else, so the keys are
spelled in one module and nowhere else — and :func:`components_status`, which asks the
one thing that record cannot answer about itself: are the components it names still the
ones whose bytes it holds?

Split from :mod:`genome.io.source`, which answers *what does this name resolve to*.
That is a question about a name; this is a question about a finished build, and the two
changed for different reasons under one roof. The seam between them is one call —
``source`` reads :meth:`ChimeraDetails.from_record` to tell a chimera's record from any
other, and nothing here reaches back.

Like ``source`` this is near the bottom of the stack: it reads records already on disk
and never fetches, writes or opens a :class:`~genome.genome.Genome`. The comparison is
**record against record** — the digests a chimera wrote down against the ones its
components' own records pin now — so it reads a handful of small JSON files and not one
base of sequence.

Examples
--------
>>> ChimeraDetails("__", (ComponentDetails("ce11", "1a2b3c", None, None),)).components
['ce11']
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genome.io.completion import (
    CompletionRecord,
    RegistrationError,
    RegistrationMismatchError,
    read_record,
    record_path,
)
from genome.io.registration import AssemblyDir, assembly_repair_command

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
        record — the shape a caller holding a record rather than a payload wants. The
        record names the assembly, so a broken one is refused in its own name without the
        caller having to supply it.

        Parameters
        ----------
        record : genome.io.completion.CompletionRecord or None
            The record to read, as :func:`~genome.io.completion.read_record` returns it.

        Returns
        -------
        ChimeraDetails or None
            The details, or ``None`` when the record is not a chimera's.

        Raises
        ------
        genome.io.completion.RegistrationError
            If the record claims to be a chimera's and cannot be read as one — see
            :meth:`from_details`.

        Examples
        --------
        >>> ChimeraDetails.from_record(None) is None
        True
        """
        return None if record is None else cls.from_details(record.details, assembly=record.name)

    @classmethod
    def from_details(cls, details: Mapping[str, Any], *, assembly: str) -> ChimeraDetails | None:
        """Read a completion record's ``details``, or ``None`` when they are not a chimera's.

        Two answers, and they are not the same answer. ``None`` means *this is not a
        chimera* — the details of an ordinary downloaded or seeded assembly, which say
        nothing about components at all. Details that *do* speak of components but cannot
        be read as a build of this package's own writing are a **broken registration**,
        and one raises rather than reading back as an ordinary assembly (ADR-0007): both
        :func:`components_status` and :func:`~genome.io.download.verify_assembly` decide
        what to check from this answer, so a chimera silently demoted to an ordinary
        assembly is one nothing ever compares against its components again.

        A chimera's record spells the separator and the components together, and each
        component entry names itself. Anything short of that — one key without the other,
        either of the wrong type, an entry that is not an object or does not name an
        assembly — is the broken case. The two annotation fields are the exception and
        stay optional: a build that registered no merged annotation writes neither, and
        both then read as ``None``.

        Taking the mapping rather than the record is what lets a caller that already holds
        one — the CLI, whose ``register`` payload *is* the record — answer from what it has
        instead of reading the same file again. ``assembly`` is what the refusal is
        addressed to, since a mapping does not know whose it is.

        Parameters
        ----------
        details : mapping of str to object
            A registration record's ``details``. Anything else it holds is ignored, and an
            empty mapping reads as *not a chimera*.
        assembly : str
            The assembly these details were recorded for, named in the refusal along with
            the command that repairs it.

        Returns
        -------
        ChimeraDetails or None
            The details, or ``None`` when they are not a chimera build's.

        Raises
        ------
        genome.io.completion.RegistrationError
            If ``details`` claim to be a chimera's and cannot be read as one. The message
            says which part did not read and quotes ``genome register <assembly>
            --force``, which is the repair.

        Examples
        --------
        >>> ChimeraDetails.from_details({}, assembly="hg38") is None
        True
        """
        separator = details.get(_SEPARATOR_KEY)
        entries = details.get(_COMPONENTS_KEY)
        if separator is None and entries is None:
            return None
        if not isinstance(separator, str):
            raise _broken(assembly, f"its {_SEPARATOR_KEY!r} is {_shown(separator)}")
        if not isinstance(entries, list):
            raise _broken(assembly, f"its {_COMPONENTS_KEY!r} is {_shown(entries)}")
        components: list[ComponentDetails] = []
        for position, entry in enumerate(entries):
            if not isinstance(entry, dict):
                raise _broken(
                    assembly,
                    f"{_COMPONENTS_KEY!r} entry {position} is {_shown(entry)} rather than "
                    f"an object",
                )
            name = entry.get(_COMPONENT_NAME_KEY)
            if not isinstance(name, str):
                raise _broken(
                    assembly,
                    f"{_COMPONENTS_KEY!r} entry {position} has no {_COMPONENT_NAME_KEY!r}, "
                    f"so it names no assembly",
                )
            components.append(
                ComponentDetails(
                    name=name,
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

    Raises
    ------
    genome.io.completion.RegistrationError
        If the record there claims to be a chimera's and cannot be read as one — see
        :meth:`ChimeraDetails.from_details`.

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
    genome.io.completion.RegistrationError
        If the record claims to be a chimera's and cannot be read as one — see
        :meth:`ChimeraDetails.from_details`.

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


def _broken(assembly: str, problem: str) -> RegistrationError:
    """Return the refusal a record that half-claims to be a chimera earns.

    ``problem`` completes "…and it cannot be read as one:", so it names the part that did
    not read rather than restating the shape.
    """
    return RegistrationError(
        f"the completion record for {assembly} says it was built from components, and it "
        f"cannot be read as one: {problem}. A chimera's record spells {_SEPARATOR_KEY!r} "
        f"and {_COMPONENTS_KEY!r} together, each component an object naming the assembly "
        f"it is, and only a build of this package's own writing produces that. Reading "
        f"this back as an ordinary assembly would leave a chimera nothing ever compares "
        f"against its components again, so it is refused instead. Build it again with "
        f"`{assembly_repair_command(assembly)}`."
    )


def _shown(value: Any) -> str:
    """Render a record field for a refusal: what it was, without quoting a whole genome."""
    return "absent" if value is None else f"{type(value).__name__} {value!r:.60}"


def _disagree(recorded: str | None, current: str | None) -> bool:
    """Whether two recorded digests are known to differ — unknown on either side is not."""
    return recorded is not None and current is not None and recorded != current


def _both_known(recorded: str | None, current: str | None) -> bool:
    """Whether two recorded digests were both there to be compared."""
    return recorded is not None and current is not None


def _text(value: Any) -> str | None:
    """Return ``value`` when it is a string, else ``None`` — a record field read loosely.

    A record is JSON somebody else's version may have written, so a field that is absent,
    null or the wrong type all read as *not known*, and none of them raises. Deliberately
    looser than the shape :meth:`ChimeraDetails.from_details` insists on: these are the
    optional fields, and a build that gathered none of them is not a broken one.
    """
    return value if isinstance(value, str) else None
