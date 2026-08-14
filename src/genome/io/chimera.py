"""Writing a chimera's FASTA and its merged annotation, from its components.

I/O boundary module, and what a **Source** of the third kind — a recipe, *these
components*, rather than a URL (ADR-0008) — is turned into. :class:`ChimeraBuilder` is an
:class:`~genome.io.registration.AssemblyRegistration` whose FASTA is concatenated from
already-prepared component assemblies instead of fetched, so it shares every other step
of a registration — the working area, placing the FASTA, deriving the companions, the
completion record — with :class:`~genome.io.download.UCSCGenomeDownloader`, and reaches
no network at any point. Two ways in: :meth:`genome.genome.Genome.chimera` for components
a caller already holds open, and :func:`build_chimera` for a name that resolved to a
component set, which opens them itself.

What a chimera's **Source** *is* — the recipe read back off its completion record, and
whether the components it names are still the ones it copied — is
:mod:`genome.io.source`, which this module writes through and re-exports from. That split
is what lets the downloader answer *is this name a chimera?* without importing anything
that builds one.

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
and :func:`components_status` is the comparison — record against record, with nothing
rehashed — made both when a finished chimera is handed back and when one is verified. An
absent digest on either side reads as *unknown* rather than as wrong.

**The merge is part of the build.** Between those two steps the **Merged annotation** is
registered, so there is no second surface to remember and ``genome register <chimera>
--force`` repairs the annotation and the FASTA together. Each component contributes its
own **Default annotation** and nothing is passed in: a component that has none contributes
nothing, and no contributors at all means no annotation rather than an empty one — while a
component naming an annotation nobody registered, or carrying several with no default,
**raises before a byte is written**, naming what closes the gap. The price, paid
knowingly, is that every chimera build now pays a ``gffutils`` database build.
:func:`~genome.io.gtf.register_merged_gtf` does the writing; what is decided here is which
annotation each component contributes and what the result is called. The name is the
contributing annotations' names joined, so a rebuild whose contributors changed writes a
different one — and the build that owns the merged annotation owns the superseded one too,
which :func:`~genome.io.gtf.discard_merged_annotation` removes, since two derived
annotations side by side leave the chimera with no default at all.

Examples
--------
>>> from genome.io.chimera import ChimeraBuilder
>>> builder = ChimeraBuilder([worm, draft])             # doctest: +SKIP
>>> builder.assembly                                    # doctest: +SKIP
'tinyCe_tinyEc'
>>> builder.build_genome().chrom_sizes.name             # doctest: +SKIP
'tinyCe_tinyEc.chrom.sizes'

Those need two prepared assemblies. What that build wrote down needs nothing, and reads
back with no disk between:

>>> from genome.io.chimera import ChimeraDetails, ComponentDetails
>>> ChimeraDetails("__", (ComponentDetails("tinyCe", None, "genes", None),)).components
['tinyCe']
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Sequence
from contextlib import ExitStack
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
from genome.io.completion import RegistrationError, read_record
from genome.io.fasta import GenomeFiles, prepare_fasta, read_chrom_sizes
from genome.io.gtf import (
    GtfAnnotation,
    MergeSource,
    discard_merged_annotation,
    register_merged_gtf,
)
from genome.io.registration import AssemblyDir, AssemblyRegistration, liulab_data_dir

# The **Source** half, which is where these are written and read. They stay importable
# from this module, which is where they used to live and where a caller reaching for
# "what is this chimera made of" still looks first.
from genome.io.source import COMPONENTS_UNCHANGED as COMPONENTS_UNCHANGED
from genome.io.source import COMPONENTS_UNKNOWN as COMPONENTS_UNKNOWN
from genome.io.source import ChimeraDetails as ChimeraDetails
from genome.io.source import ComponentDetails as ComponentDetails
from genome.io.source import components_status, is_prepared, merged_annotation_name
from genome.io.source import read_chimera_details as read_chimera_details

if TYPE_CHECKING:
    from genome.genome import Genome

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


def build_chimera(
    assembly_dir: AssemblyDir, components: Sequence[str], *, overwrite: bool = False
) -> GenomeFiles:
    """Concatenate ``components`` into the assembly at ``assembly_dir``, opening each to do it.

    What a :class:`~genome.io.source.ComponentSource` is turned into: the registration
    resolved a name to a component set and this is the only thing that can act on one, since
    a component is read as a whole :class:`~genome.genome.Genome`.

    The gate on a name that resolves to a chimera is here. Every component must already be
    prepared under the shared data root: a machine that is missing one **raises**, naming it
    and the command that prepares it, rather than turning one mistyped string into a
    whole-genome download per part — which is the cost of guessing wrong, and the reason the
    name alone never fetches anything.

    Each component is opened and closed again around the build, so a build that raises
    leaves no ``.2bit`` handle behind.

    Parameters
    ----------
    assembly_dir : genome.io.registration.AssemblyDir
        Where the chimera is built. Its name must be the one ``components`` derive, which
        is what resolving the name guaranteed.
    components : sequence of str
        The **Component** assembly names, in the sorted order the chimera's name spells
        them.
    overwrite : bool, default False
        Build again from scratch — the repair for a directory that raises.

    Returns
    -------
    genome.io.fasta.GenomeFiles
        Paths to the written FASTA and its three derived files.

    Raises
    ------
    FileNotFoundError
        If any component is not registered on this machine; the message names each missing
        one and the command that prepares it.

    See Also
    --------
    ChimeraBuilder.build_genome : everything this then runs, and what else it raises.

    Examples
    --------
    >>> from genome.io.registration import AssemblyDir
    >>> build_chimera(AssemblyDir.locate("ce11_ecHT115"), ("ce11", "ecHT115"))  # doctest: +SKIP
    GenomeFiles(fasta=PosixPath('.../ce11_ecHT115.fa'), ...)
    """
    # Deferred, and not for a layering slip: a component is a whole prepared genome, so
    # opening one is the top of the stack reaching back down — `genome.genome` imports this
    # module to build a chimera from components a caller already holds.
    from genome.genome import Genome

    missing = [name for name in components if not is_prepared(name)]
    if missing:
        listed = ", ".join(missing)
        commands = ", ".join(f"`genome register {name}`" for name in missing)
        raise FileNotFoundError(
            f"{assembly_dir.assembly} is a chimera of {', '.join(components)}, and a chimera "
            f"copies the bytes of components that are already prepared: {listed} is not "
            f"registered under {liulab_data_dir()}. Nothing was downloaded, because "
            f"naming a chimera is not a way to ask for its parts. Prepare what is "
            f"missing with {commands}, then run this again."
        )
    with ExitStack() as opened:
        genomes = [opened.enter_context(Genome(name, progressbar=False)) for name in components]
        files = ChimeraBuilder(genomes, assembly_dir.path).build_genome(overwrite=overwrite)
    return files


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

    Every one of those wants two prepared assemblies. The naming contract wants none: it
    is checked on the arguments alone, before this even has a directory to fill.

    >>> try:
    ...     ChimeraBuilder([])
    ... except ChimeraNamingError:
    ...     print("a chimera of nothing is refused before anything is read")
    a chimera of nothing is refused before anything is read
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
        annotation**, and write the completion record last. A build whose contributors
        changed since the last one merges under a different name, and the annotation the
        previous record names is removed rather than left registered beside the new one.

        A chimera whose record says it finished is returned from that record without
        rewriting anything — the annotation included, since a finished chimera already has
        the one its build wrote — but only once its components are shown to still be the
        ones it was built from (:func:`components_status`). A directory that cannot be
        trusted **raises** instead of being rebuilt, naming ``genome register <name>
        --force`` (ADR-0007). That command is what ``overwrite=True`` is, and it repairs
        both halves in one pass.

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
            # Its refusal is what is wanted; the answer it returns is for a surface to
            # print, and there is none here.
            components_status(self.dir)
            return registered
        # Settled here rather than beside the merge, so that a component naming an
        # annotation this machine has not registered costs nothing: the two refusals below
        # are a caller's mistake, and finding out after a whole genome had been written
        # would leave a directory nothing vouches for behind. It is *after* the early
        # return above on purpose — reopening a finished chimera must not depend on its
        # components' annotations still being registered.
        contributions = self._contributions()
        # Read while the previous build's record is still the one on disk: the name it
        # merged under is written down there and nowhere else, and nothing can re-derive
        # it once the contributors have changed.
        superseded = self._recorded_merged_annotation()
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
        self._discard_superseded(superseded, merged)
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

    def _recorded_merged_annotation(self) -> str | None:
        """Return the **Merged annotation** the last build here registered, or ``None``.

        :attr:`ChimeraDetails.merged_annotation` of the record this build is about to
        overwrite. ``None`` for a directory no chimera build has finished in, and for one
        whose build merged nothing.
        """
        details = read_chimera_details(self.cache_dir)
        return None if details is None else details.merged_annotation

    def _discard_superseded(self, previous: str | None, merged: GtfAnnotation | None) -> None:
        """Remove the merged annotation the previous build registered and this one replaced.

        The build owns the merged annotation, so it owns the stale one too. The name is
        derived from what contributed, so a contributing set that changed across a rebuild
        changes it — a component whose default annotation is a different one now, or one
        that contributes nothing where it used to. The old registration would otherwise
        stay beside the new one, leaving the chimera carrying two derived annotations with
        neither flagged and therefore **no default at all**, which is a chimera that
        arrived annotated coming back from a legitimate repair with none.

        Surgical in both directions: only the name the previous record itself names, and
        only when that directory's own record shows a merge wrote it — an annotation a
        caller registered by hand is never a build's to remove.
        """
        current = merged.name if merged is not None else None
        if previous is None or previous == current:
            return
        discard_merged_annotation(self.cache_dir, previous)

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

        The separator its chromosome names carry, and one entry per component, in the
        shape :class:`~genome.io.source.ChimeraDetails` reads back — which is where the
        keys are spelled, so the writing and the reading cannot drift.
        """
        return ChimeraDetails(
            separator=self.separator,
            component_details=tuple(
                self._component_details(component, contributions[component.assembly])
                for component in self.components
            ),
        ).as_details(merged=merged)

    def _component_details(
        self, component: Genome, contribution: _Contribution | None
    ) -> ComponentDetails:
        """Return one component's facts: what it was, and what it contributed.

        Record to record throughout — the component's completion record already holds the
        sha256 of its FASTA and its annotation's holds that annotation's, so nothing is
        rehashed here. ``None`` for a digest a record pinned none of, which reads as
        unknown rather than as wrong.
        """
        # The component's own **Assembly dir**, which is where its FASTA sits. Asked of
        # the component rather than derived from its name, since a component may be
        # registered somewhere other than the shared data root.
        record = read_record(component.fasta_path.parent)
        return ComponentDetails(
            name=component.assembly,
            sha256=None if record is None else record.sha256,
            annotation=None if contribution is None else contribution.annotation,
            annotation_sha256=None if contribution is None else contribution.sha256,
        )


def _merged_name(sources: Sequence[MergeSource]) -> str:
    """Return the **Registered name** the merge of ``sources`` is filed under.

    :func:`~genome.io.source.merged_annotation_name` over what contributed, in the
    sorted-component order the chimera's own name spells. The join lives there because a
    record read back has to recover the same name from the same parts.
    """
    return merged_annotation_name([source.annotation for source in sources])


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
