"""Writing a chimera's FASTA and its merged annotation, from its components.

I/O boundary module, and the third kind of **Source** an assembly can have — a recipe,
*these components*, rather than a URL (ADR-0008). :class:`ChimeraBuilder` is an
:class:`~genome.io.registration.AssemblyRegistration` whose FASTA is concatenated from
already-prepared component assemblies instead of fetched, so it shares every other step
of a registration — the working area, placing the FASTA, deriving the companions, the
completion record — with :class:`~genome.io.download.UCSCGenomeDownloader`, and reaches
no network at any point. The way in is :meth:`genome.genome.Genome.chimera`.

**Nothing here reimplements a tool.** The rule is to shell out, and the burden is
discharged on evidence rather than waived: no installed tool renames a FASTA header.
``samtools faidx`` cannot rename anything and rewraps what it writes, ``bedtools
getfasta`` emits each sequence as one 20 Mb line, and seqkit would be a new dependency
that also rewraps every byte. Rewrapping is not cosmetic here, so the concatenation is a
streaming pass in this module — the only thing in the build that is not a native binary.

**Sequence lines are copied verbatim.** The only line that changes is a header, and the
only part of it that changes is the name it gives its sequence — the first non-whitespace
token after ``>``, which is where ``samtools faidx`` and ``faToTwoBit`` read a name from,
whitespace after the ``>`` skipped. ``>I some description`` becomes ``>I__tinyCe some
description``. The suffix rides on that token because both STAR and chromap truncate a
FASTA header at the first whitespace (measured in
``docs/research/aligner-index-params-and-reference-names.md``). Leaving every
sequence byte alone is what makes components that disagree about wrap width and about
soft-masking correct by default rather than by special case: ``ce11`` wraps at 60 and is
21.95% lower case, ``ecHT115`` wraps at 80 and is not masked at all, and a chimera of
the two is heterogeneous in both respects because neither was normalized towards the
other.

The FASTA is written into the working area and hashed **in the same pass**, so a whole
genome is neither held in memory nor read a second time to produce the digest the record
carries. Once ``prepare_fasta`` has derived the companions, the built ``chrom.sizes`` is
compared against the one the components predict — every chromosome, in component-sorted
then file order, suffixed, with its own component's length — and only if that agrees is
the completion record written: ``source_url`` is ``None`` because nothing was fetched,
and ``sha256`` is the digest of the bytes this module just wrote, which is what a later
verification falls back to.

**What no digest of a chimera can see** is a component registered again underneath it:
the chimera's own bytes are untouched and still agree with its record, while the
component they were copied from is gone. That is what the digests in ``details`` are for,
and :func:`check_components_unchanged` is the comparison — record against record, with
nothing rehashed — made both when a finished chimera is handed back and when one is
verified. An absent digest on either side reads as *unknown* rather than as wrong.

**The merge is part of the build.** Between those two steps the **Merged annotation** is
registered, so there is no second surface to remember and ``genome register <chimera>
--force`` repairs the annotation and the FASTA together. Each component contributes its
own **Default annotation** and nothing is passed in: a component that has none contributes
nothing, and no contributors at all means no annotation rather than an empty one — while a
component naming an annotation nobody registered, or carrying several with no default,
**raises before a byte is written**, naming what closes the gap. The price, paid
knowingly, is that every chimera build now pays a ``gffutils`` database build.
:func:`~genome.io.gtf.register_merged_gtf` does the writing; what is decided here is which
annotation each component contributes and what the result is called.

Examples
--------
>>> from genome.io.chimera import ChimeraBuilder
>>> builder = ChimeraBuilder([worm, draft])             # doctest: +SKIP
>>> builder.assembly                                    # doctest: +SKIP
'tinyCe_tinyEc'
>>> builder.build_genome().chrom_sizes.name             # doctest: +SKIP
'tinyCe_tinyEc.chrom.sizes'
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from genome.chimera import (
    ChimeraNamingError,
    check_roundtrip,
    derive_name,
    derive_separator,
    suffixed,
)
from genome.io.completion import (
    CompletionRecord,
    RegistrationError,
    RegistrationMismatchError,
    read_record,
    record_path,
)
from genome.io.fasta import GenomeFiles, prepare_fasta, read_chrom_sizes
from genome.io.gtf import GtfAnnotation, MergeSource, annotation_dir, register_merged_gtf
from genome.io.registration import (
    AssemblyRegistration,
    assembly_data_dir,
    assembly_repair_command,
)

if TYPE_CHECKING:
    from genome.genome import Genome

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

#: A FASTA header split into the whitespace that follows ``>``, the sequence name after
#: it, and everything after that. The leading run is its own group because ``samtools
#: faidx`` and ``faToTwoBit`` both *skip* it and name the sequence from the first
#: non-whitespace token — measured, not assumed — so a header ``> chrA desc`` names
#: ``chrA``, and a pattern that read the token as empty would suffix the wrong thing.
#: ``DOTALL`` so the trailing newline lands in the remainder and is written back
#: unchanged; ASCII ``\s``/``\S`` because a bytes pattern is ASCII by definition.
_HEADER_RE = re.compile(rb">(\s*)(\S*)(.*)", re.DOTALL)


class AmbiguousDefaultAnnotationError(ValueError):
    """A component carries several annotations and nothing says which one a chimera takes.

    A **Component** contributes its own **Default annotation**, and a component with
    several registered and none flagged by the annotation table has no default at all —
    which is the ordinary, deliberate answer to *pick one for me* everywhere else in the
    package, and the one place it cannot stand. Guessing would put a set of gene models
    into a merged annotation nobody chose, under a name that would look identical to the
    one the caller meant.

    The message names ``default_gtf=``, which is how a caller says which, and it says so
    of that component rather than of the chimera: the fix is one argument on one
    component's constructor.

    A :class:`ValueError` because the component handed in is not one this build can use.

    Examples
    --------
    >>> raise AmbiguousDefaultAnnotationError("ce11 carries 2 annotations ...")
    Traceback (most recent call last):
    genome.io.chimera.AmbiguousDefaultAnnotationError: ce11 carries 2 annotations ...
    """


@dataclass(frozen=True)
class ComponentDetails:
    """What a chimera's completion record says about one of its components.

    One entry of the ``details`` shape :class:`ChimeraBuilder` writes, read back. Every
    field is a fact about the component *at the time this chimera was built*, taken from
    that component's own records rather than by rehashing anything — which is what lets a
    later pass notice a component re-registered underneath the chimera.

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


@dataclass(frozen=True)
class ChimeraDetails:
    """What a chimera's completion record says about the build that produced it.

    The reader of the ``details`` shape :class:`ChimeraBuilder` writes, so that nothing
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

    @classmethod
    def from_record(cls, record: CompletionRecord | None) -> ChimeraDetails | None:
        """Read a completion record's chimera details, or ``None`` when it has none.

        ``None`` means *this is not a chimera* — the record of an ordinary downloaded or
        seeded assembly, an absent record, or one whose ``details`` do not carry the
        shape a chimera build writes. Those read alike on purpose: nothing but a build of
        this package's own writing may make an assembly answer as a chimera.

        The two annotation fields are optional, and a record written by a build that
        registered no merged annotation carries neither; both then read as ``None``.

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
        if record is None:
            return None
        separator = record.details.get(_SEPARATOR_KEY)
        entries = record.details.get(_COMPONENTS_KEY)
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


def check_components_unchanged(directory: Path, assembly: str) -> None:
    """Raise unless every component of the assembly in ``directory`` is what it was.

    The one failure a digest of a chimera's own bytes cannot show. Those bytes are a copy
    of its components', so a component registered again underneath it leaves the chimera
    intact, agreeing with its own record, and no longer a copy of anything that exists —
    silently stale sequence, and stale gene models one level down. Both are caught here,
    and both are **record against record**: the digests this chimera wrote down are
    compared against the ones the components' own records pin now, so this reads a
    handful of small JSON files and not one base of sequence.

    An assembly with no components recorded has nothing to compare and returns at once,
    which is what makes an ordinary assembly pay nothing rather than be asked about.
    Likewise an absent digest on either side, which means *unknown* rather than wrong: a
    component that pinned none, or an annotation registered before its digest was
    recorded, leaves that component unguarded rather than refused.

    Parameters
    ----------
    directory : pathlib.Path
        The **Assembly dir** the chimera was built in.
    assembly : str
        Its assembly name, quoted in the error along with the command that repairs it.

    Raises
    ------
    genome.io.completion.RegistrationMismatchError
        If a component's FASTA, or the annotation it contributed to the **Merged
        annotation**, is not the one this chimera was built from. The message names both
        digests and ``genome register <assembly> --force``, which rebuilds both halves.

    Examples
    --------
    >>> from pathlib import Path
    >>> check_components_unchanged(Path("/tmp/definitely-not-a-build"), "notAChimera")
    """
    details = read_chimera_details(directory)
    if details is None:
        return
    for entry in details.component_details:
        _check_component_unchanged(entry, assembly=assembly, directory=directory)


def _check_component_unchanged(entry: ComponentDetails, *, assembly: str, directory: Path) -> None:
    """Raise unless one component — its FASTA, then its annotation — is still itself.

    The component is found by name under the shared data root, exactly as an aligner index
    finds the assembly it was built from: a chimera's record carries component names and
    no paths, which is what keeps a registered directory movable.
    """
    component_dir = assembly_data_dir(entry.name)
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
    if entry.annotation is None:
        return
    gtf_dir = annotation_dir(component_dir, entry.annotation)
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


def _disagree(recorded: str | None, current: str | None) -> bool:
    """Whether two recorded digests are known to differ — unknown on either side is not."""
    return recorded is not None and current is not None and recorded != current


@dataclass(frozen=True)
class _Contribution:
    """What one component puts into the **Merged annotation**, settled before any writing.

    ``None`` in place of one of these is a component that contributes nothing, which is
    an ordinary shape rather than a failure: an assembly may carry sequences no annotation
    mentions, and the chromosome check is one-directional for exactly that reason.
    """

    annotation: str
    gtf: Path
    sha256: str | None


class ChimeraBuilder(AssemblyRegistration):
    """Concatenate prepared assemblies into one, and register the result.

    An :class:`~genome.io.registration.AssemblyRegistration` whose FASTA is written here
    rather than fetched: the base owns the assembly's directory and the steps that finish
    a registration in it, and what is added is the concatenation — which components, in
    which order, with which chromosome names.

    Its name is not given but derived from the component set (ADR-0008), so one set of
    components means one directory whatever order they were passed in. Its chromosome
    names are suffixed unconditionally with the component they came from (ADR-0009),
    under a separator derived from those components and recorded, and every one of those
    names is proved to read back **before** a byte is written.

    Parameters
    ----------
    components : sequence of genome.genome.Genome
        Two or more prepared component assemblies, in any order, each appearing once. A
        :class:`~genome.genome.Genome` is prepared by construction, which is what makes
        most of the invariant a signature rather than a check; what is checked here is
        what a type cannot say — that there are at least two of them, that none repeats,
        and that none is itself a chimera.
    cache_dir : str or pathlib.Path, optional
        Override the directory the chimera is built in. Defaults to
        :func:`assembly_data_dir(<derived name>)
        <genome.io.registration.assembly_data_dir>`.

    Attributes
    ----------
    assembly : str
        The derived chimera name — the component names sorted and joined by ``_``.
    cache_dir : pathlib.Path
        The **Assembly dir** this build fills.
    components : tuple of genome.genome.Genome
        The components, sorted by assembly name — the order their sequences are written
        in, and the order the derived name spells them.
    separator : str
        The run of underscores every chromosome name written here carries.

    Raises
    ------
    genome.chimera.ChimeraNamingError
        If fewer than two components are given, a component repeats, a component's name
        is not alphanumeric, or a component is itself a chimera.

    Examples
    --------
    >>> builder = ChimeraBuilder([worm, draft])            # doctest: +SKIP
    >>> builder.assembly, builder.separator                # doctest: +SKIP
    ('tinyCe_tinyEc', '__')
    >>> files = builder.build_genome()                     # doctest: +SKIP
    """

    def __init__(self, components: Sequence[Genome], cache_dir: str | Path | None = None) -> None:
        _check_not_nested(components)
        names = [component.assembly for component in components]
        super().__init__(derive_name(names), cache_dir)
        self.components: tuple[Genome, ...] = tuple(
            sorted(components, key=lambda component: component.assembly)
        )
        chromosomes = {
            component.assembly: [str(name) for name in component.chromosomes]
            for component in self.components
        }
        self.separator: str = derive_separator(chromosomes)
        # Before a byte is written, as the contract intends: a name that could not be
        # read back is a reference nobody can attribute, and it costs 94 string
        # operations to know for the two assemblies the lab ships.
        check_roundtrip(chromosomes, self.separator)

    def build_genome(self, *, overwrite: bool = False) -> GenomeFiles:
        """Write the chimera's FASTA and its merged annotation, and record that it finished.

        The whole build: settle which annotation each component contributes, concatenate
        the components' bytes into the working area while hashing them, place the result
        as ``<assembly>.fa``, derive the ``.fai``, ``.2bit`` and ``chrom.sizes`` with
        :func:`~genome.io.fasta.prepare_fasta`, check the sequence names and lengths that
        came back against the ones the components predict, register the **Merged
        annotation**, and write the completion record last.

        A chimera whose record says it finished is returned from that record without
        rewriting anything — the annotation included, since a finished chimera already has
        the one its build wrote — but only once its components are shown to still be the
        ones it was built from (:func:`check_components_unchanged`). A directory that
        cannot be trusted **raises** instead of being rebuilt, naming ``genome register
        <name> --force`` (ADR-0007). That command is what ``overwrite=True`` is, and it
        repairs both halves in one pass.

        Parameters
        ----------
        overwrite : bool, default False
            Build again from scratch: the completion record is not consulted and every
            derived file is rebuilt even when it looks fresh. The repair for a directory
            that raises.

        Returns
        -------
        genome.io.fasta.GenomeFiles
            Paths to the written FASTA and its three derived files.

        Raises
        ------
        genome.chimera.ChimeraNamingError
            If a component's FASTA carries a header that names no sequence, so there is
            nothing for the suffix to ride on.
        genome.io.gtf.AnnotationNotRegisteredError
            If a component's **Default annotation** names something nobody registered on
            this machine. Raised before anything is written, and the message names the
            command that registers it.
        AmbiguousDefaultAnnotationError
            If a component carries several annotations and none is its default. Likewise
            before anything is written.
        genome.io.completion.UnfinishedRegistrationError
            If the chimera's directory holds files but no record.
        genome.io.completion.RegistrationMismatchError
            If its record disagrees with what is on disk, or if a component was
            registered again since this chimera was built, so the copy of that
            component's bytes held here is of bytes that no longer exist.
        genome.io.completion.RegistrationError
            If the built ``chrom.sizes`` is not the concatenation the components
            predict. Nothing vouches for the directory in that case — the record is
            written after this check, never before.
        genome.io.gtf.ChromosomeMismatchError
            If a merged seqname is not one the built FASTA carries — the two halves of
            the build disagreeing, which nothing else would catch.
        genome.external.ToolNotFoundError
            If ``samtools``, ``faToTwoBit`` or ``twoBitInfo`` are not on ``PATH``.
        RuntimeError
            If any native preparation tool exits non-zero.

        Examples
        --------
        >>> ChimeraBuilder([worm, draft]).build_genome()    # doctest: +SKIP
        GenomeFiles(fasta=PosixPath('.../tinyCe_tinyEc.fa'), ...)
        """
        registered = self._completed_genome(overwrite=overwrite, repair=self._repair_command())
        if registered is not None:
            # The record vouches for this directory's own files and can say nothing about
            # the components those files were copied from, so the one question it cannot
            # answer is asked here — before a stale chimera is handed back as finished.
            check_components_unchanged(self.cache_dir, self.assembly)
            return registered
        # Settled here rather than beside the merge, so that a component naming an
        # annotation this machine has not registered costs nothing: the two refusals below
        # are a caller's mistake, and finding out after a whole genome had been written
        # would leave a directory nothing vouches for behind. It is *after* the early
        # return above on purpose — reopening a finished chimera must not depend on its
        # components' annotations still being registered.
        contributions = self._contributions()
        work = self._work_dir
        work.mkdir(parents=True, exist_ok=True)
        staged = work / f"{self.assembly}.fa"
        digest = self._concatenate(staged)
        files = prepare_fasta(self._place_fasta(staged), overwrite=overwrite)
        self._check_built_names(files)
        # The merged annotation lands here: the chimera's own chrom.sizes now exists for
        # the merge to be checked against, and this is above the record, whose details
        # describe what the merge did.
        merged = self._merge_annotation(files, contributions)
        self._record_completion(
            files,
            source_url=None,
            sha256=digest,
            details=self._details(contributions, merged=merged is not None),
        )
        return files

    def _contributions(self) -> dict[str, _Contribution | None]:
        """Return what each component contributes to the merged annotation, in order.

        Keyed by component assembly name in sorted order, ``None`` for a component that
        contributes nothing. Nothing is written and nothing is read but records, so this
        is the cheap half of the build and the right place for it to refuse.
        """
        return {component.assembly: self._contribution(component) for component in self.components}

    def _contribution(self, component: Genome) -> _Contribution | None:
        """Return one component's contribution, or ``None`` when it has none.

        The intention, spelled out. A component's **Default annotation** is whatever that
        component already decided — an explicit ``default_gtf=``, the annotation table's
        flag, or the sole registered one — so a chimera needs no annotation argument of
        its own. Three answers rather than two: the component has one and it is
        registered here; it has none, and contributes nothing; or it has one in name only,
        which is a cold machine and raises with the command that fixes it.
        """
        if component.default_gtf is None:
            if component.annotations:
                raise AmbiguousDefaultAnnotationError(
                    f"component {component.assembly!r} has {len(component.annotations)} "
                    f"annotations registered and no default among them: "
                    f"{', '.join(component.annotations)}. A component contributes its own "
                    f"default annotation to a chimera's merged one, so this build cannot "
                    f"tell which set of gene models you meant. Open that component with "
                    f"Genome({component.assembly!r}, default_gtf=<name>) — naming one of "
                    f"the above — and build the chimera from it."
                )
            return None
        # Raises AnnotationNotRegisteredError, naming what registers it, when the default
        # is a name and not yet a file: the same refusal a cold machine gets for a
        # component it has never prepared, one level down.
        gtf = component.get_gtf_path(component.default_gtf)
        record = read_record(gtf.parent)
        return _Contribution(
            annotation=component.default_gtf,
            gtf=gtf,
            sha256=None if record is None else record.sha256,
        )

    def _merge_annotation(
        self, files: GenomeFiles, contributions: dict[str, _Contribution | None]
    ) -> GtfAnnotation | None:
        """Register the merged annotation, or return ``None`` when nothing contributed.

        No contributors means **no annotation registered**, rather than an empty one that
        would answer every query with nothing while looking perfectly healthy — the same
        reason a GTF whose chromosomes are all unknown is refused rather than kept.
        """
        sources = [
            MergeSource(component=component, annotation=entry.annotation, gtf=entry.gtf)
            for component, entry in contributions.items()
            if entry is not None
        ]
        if not sources:
            return None
        return register_merged_gtf(
            self.cache_dir,
            _merged_name(sources),
            sources,
            separator=self.separator,
            chrom_sizes=files.chrom_sizes,
        )

    def _concatenate(self, destination: Path) -> str:
        """Write every component's FASTA to ``destination`` and return the sha256 of it.

        One pass, a line at a time, hashing what is written as it is written — so a whole
        genome is never held in memory and never read again to produce the digest the
        record carries. Sequence lines go through untouched; a header line has its first
        token suffixed and nothing else about it changed. A component whose last line
        carries no newline gets one, since otherwise the next component's header would
        land on the end of its sequence.
        """
        digest = hashlib.sha256()
        with destination.open("wb") as output:
            for component in self.components:
                terminated = True
                with component.fasta_path.open("rb") as source:
                    for line in source:
                        written = (
                            _extend_header(line, component.assembly, self.separator)
                            if line.startswith(b">")
                            else line
                        )
                        output.write(written)
                        digest.update(written)
                        terminated = line.endswith(b"\n")
                if not terminated:
                    output.write(b"\n")
                    digest.update(b"\n")
        return digest.hexdigest()

    def _expected_chrom_sizes(self) -> list[tuple[str, int]]:
        """Return the ``(name, length)`` pairs this chimera's components predict, in order.

        Every component's chromosomes in the order that component declares them,
        suffixed, with the component's own lengths — components taken in sorted order,
        which is the order :meth:`_concatenate` writes them in.
        """
        return [
            (suffixed(str(name), component.assembly, self.separator), int(length))
            for component in self.components
            for name, length in component.chrom_sizes.items()
        ]

    def _check_built_names(self, files: GenomeFiles) -> None:
        """Raise unless the built ``chrom.sizes`` is what the components predict.

        The one thing between a concatenation and a record that vouches for it. It
        compares names *and* lengths *and* order, reading the file the native tools
        actually produced rather than the bytes this module thinks it wrote — so a header
        that did not get suffixed, a sequence that did not survive, or a component that
        was read short is caught here instead of by whatever is built on it later.
        """
        expected = self._expected_chrom_sizes()
        built = [
            (str(name), int(length)) for name, length in read_chrom_sizes(files.chrom_sizes).items()
        ]
        if built == expected:
            return
        raise RegistrationError(
            f"the {self.assembly} FASTA just written does not carry the sequences its "
            f"components predict: {_first_difference(built, expected)}. No "
            f"{files.chrom_sizes.name} that disagrees is ever recorded, so nothing "
            f"vouches for {self.cache_dir} and it reads as an interrupted build. Build it "
            f"again with `{self._repair_command()}` — and if it fails the same way, "
            f"report both lists, since bytes copied verbatim under a derived name are "
            f"meant to make this impossible."
        )

    def _details(
        self, contributions: dict[str, _Contribution | None], *, merged: bool
    ) -> dict[str, Any]:
        """Return what this build records about itself beyond the files it wrote.

        The separator its chromosome names carry, and one entry per component. Kept to
        facts a later pass cannot re-derive from the name alone: which spelling was used,
        what each component was when this chimera was built, and — when there is a merged
        annotation — which of each component's annotations went into it.
        """
        return {
            _SEPARATOR_KEY: self.separator,
            _COMPONENTS_KEY: [
                self._component_details(component, contributions[component.assembly], merged=merged)
                for component in self.components
            ],
        }

    def _component_details(
        self, component: Genome, contribution: _Contribution | None, *, merged: bool
    ) -> dict[str, Any]:
        """Return one component's entry: what it was, and what it contributed.

        Record to record throughout — the component's completion record already holds the
        sha256 of its FASTA and its annotation's holds that annotation's, so nothing is
        rehashed here. ``None`` for a digest a record pinned none of, which reads as
        unknown rather than as wrong. The two annotation keys are written only when a
        merged annotation was registered; a build that registered none says nothing about
        annotations at all rather than writing ``null`` beside every component.
        """
        # The component's own **Assembly dir**, which is where its FASTA sits. Asked of
        # the component rather than derived from its name, since a component may be
        # registered somewhere other than the shared data root.
        record = read_record(component.fasta_path.parent)
        entry: dict[str, Any] = {
            _COMPONENT_NAME_KEY: component.assembly,
            _COMPONENT_DIGEST_KEY: None if record is None else record.sha256,
        }
        if merged:
            entry[_COMPONENT_ANNOTATION_KEY] = (
                None if contribution is None else contribution.annotation
            )
            entry[_COMPONENT_ANNOTATION_DIGEST_KEY] = (
                None if contribution is None else contribution.sha256
            )
        return entry


def _merged_name(sources: Sequence[MergeSource]) -> str:
    """Return the **Registered name** the merge of ``sources`` is filed under.

    The contributing annotations' names joined by ``+``, in the sorted-component order the
    chimera's own name spells — ``wormbase_ws298+refseq_rs_2025_06_26``. Derived, like
    everything else about a chimera, so it changes the moment any component's default
    annotation does and a database built from the old set can never be found under it.

    It needs no parse-back: what a merged annotation is made of is recovered from the
    components, and written down in its own record besides. And it is not asked to carry
    *which* components contributed — a chimera with a component that contributes nothing
    spells the same name a different subset would — which is why
    :func:`~genome.io.gtf.register_merged_gtf` adopts nothing from disk and writes the
    annotation every time it runs.
    """
    return _ANNOTATION_JOIN.join(source.annotation for source in sources)


def _text(value: Any) -> str | None:
    """Return ``value`` when it is a string, else ``None`` — a record field read loosely.

    A record is JSON somebody else's version may have written, so a field that is absent,
    null or the wrong type all read as *not known*, and none of them raises.
    """
    return value if isinstance(value, str) else None


def _check_not_nested(components: Sequence[Genome]) -> None:
    """Raise if any of ``components`` is itself a chimera.

    The half of the no-nesting rule that :func:`genome.chimera.derive_name` cannot make:
    it refuses a name spelled like a chimera's, while whether a *prepared* assembly is
    one is answered by the record on its disk, which the naming contract cannot see.
    """
    for component in components:
        nested = component.components
        if nested is not None:
            raise ChimeraNamingError(
                f"component {component.assembly!r} is itself a chimera, of "
                f"{', '.join(nested)}; a component is always a canonical assembly, so "
                f"nesting is forbidden by the model rather than deferred (ADR-0008). "
                f"List those components alongside the others instead."
            )


def _extend_header(line: bytes, component: str, separator: str) -> bytes:
    """Return ``line`` with the name it gives its sequence suffixed by ``component``.

    ``b">I some description"`` becomes ``b">I__tinyCe some description"``. The name is
    the first non-whitespace token after ``>``, which is the sequence name **as the
    native tools read it**: ``samtools faidx`` and ``faToTwoBit`` both skip whitespace
    that follows ``>``, so ``b"> desc"`` declares a sequence called ``desc`` and that is
    what gets suffixed. Every other byte of the line — the skipped whitespace, the
    description, and the line ending, whatever it is — is written back untouched, because
    none of it is what a tool reads: STAR and chromap truncate a header at the first
    whitespace.

    A header that names nothing raises rather than writing a name made only of separator
    and component, which the naming contract cannot read back and which ``samtools
    faidx`` would file under the empty name.
    """
    match = _HEADER_RE.match(line)
    if match is None:  # pragma: no cover - the caller only passes lines starting with '>'
        return line
    leading, name, remainder = match.groups()
    if not name:
        raise ChimeraNamingError(
            f"the {component} FASTA carries a header that names no sequence: {line!r}. A "
            f"sequence is named by the first non-whitespace token after '>' — that is "
            f"what samtools faidx and faToTwoBit read — so there is nothing here to "
            f"suffix with {component!r} and nothing a chimera could attribute. Name every "
            f"sequence in that FASTA and register {component!r} again."
        )
    spelled = suffixed(name.decode("utf-8", "surrogateescape"), component, separator)
    return b">" + leading + spelled.encode("utf-8", "surrogateescape") + remainder


def _first_difference(built: list[tuple[str, int]], expected: list[tuple[str, int]]) -> str:
    """Describe how two ``(name, length)`` lists disagree, in one line."""
    for index, (found, wanted) in enumerate(zip(built, expected, strict=False)):
        if found != wanted:
            return (
                f"sequence {index + 1} is {found[0]!r} at {found[1]} bp where "
                f"{wanted[0]!r} at {wanted[1]} bp was expected"
            )
    return f"it carries {len(built)} sequences where {len(expected)} were expected"
