"""Writing a chimera's FASTA: its components' bytes, one token per header extended.

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
only part of it that changes is the first whitespace-delimited token: ``>I some
description`` becomes ``>I__tinyCe some description``. The suffix rides on that token
because both STAR and chromap truncate a FASTA header at the first whitespace (measured
in ``docs/research/aligner-index-params-and-reference-names.md``). Leaving every
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
from genome.io.completion import CompletionRecord, RegistrationError, read_record
from genome.io.fasta import GenomeFiles, prepare_fasta, read_chrom_sizes
from genome.io.registration import AssemblyRegistration

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

#: A FASTA header split into its first whitespace-delimited token and everything after
#: it. ``DOTALL`` so the trailing newline lands in the remainder and is written back
#: unchanged; ASCII ``\S`` because a bytes pattern is ASCII by definition.
_HEADER_RE = re.compile(rb">(\S*)(.*)", re.DOTALL)


@dataclass(frozen=True)
class ChimeraDetails:
    """What a chimera's completion record says about the build that produced it.

    The reader of the ``details`` shape :class:`ChimeraBuilder` writes, so that nothing
    else has to know its keys. It answers the only question that decides whether an
    assembly is a **Chimera** at runtime — the record, never the metadata row — and
    carries the two facts a later pass needs: the separator its chromosome names were
    written with, and what each component hashed to at build time.

    Attributes
    ----------
    separator : str
        The run of underscores this chimera's chromosome names carry.
    component_digests : dict of str to (str or None)
        Each **Component** assembly name, in the sorted order the chimera's name spells
        them, mapped to the sha256 its own completion record pinned when this chimera was
        built. ``None`` for a component whose record pinned none, which reads as
        *unknown* rather than as *wrong*.

    Examples
    --------
    >>> details = ChimeraDetails("__", {"ce11": "1a2b3c", "ecHT115": "4d5e6f"})
    >>> details.components
    ['ce11', 'ecHT115']
    """

    separator: str
    component_digests: dict[str, str | None]

    @property
    def components(self) -> list[str]:
        """The component assembly names, sorted — a fresh list each call."""
        return list(self.component_digests)

    @classmethod
    def from_record(cls, record: CompletionRecord | None) -> ChimeraDetails | None:
        """Read a completion record's chimera details, or ``None`` when it has none.

        ``None`` means *this is not a chimera* — the record of an ordinary downloaded or
        seeded assembly, an absent record, or one whose ``details`` do not carry the
        shape a chimera build writes. Those read alike on purpose: nothing but a build of
        this package's own writing may make an assembly answer as a chimera.

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
        digests: dict[str, str | None] = {}
        for entry in entries:
            if not isinstance(entry, dict) or not isinstance(entry.get(_COMPONENT_NAME_KEY), str):
                return None
            digest = entry.get(_COMPONENT_DIGEST_KEY)
            digests[entry[_COMPONENT_NAME_KEY]] = digest if isinstance(digest, str) else None
        return cls(separator=separator, component_digests=digests)


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
        """Write the chimera's FASTA, prepare it, and record that it finished.

        The whole build: concatenate the components' bytes into the working area while
        hashing them, place the result as ``<assembly>.fa``, derive the ``.fai``,
        ``.2bit`` and ``chrom.sizes`` with :func:`~genome.io.fasta.prepare_fasta`, check
        the sequence names and lengths that came back against the ones the components
        predict, and write the completion record last.

        A chimera whose record says it finished is returned from that record without
        rewriting anything, exactly as a downloaded assembly is; a directory that cannot
        be trusted **raises** instead of being rebuilt, naming
        ``genome register <name> --force`` (ADR-0007). That command is what
        ``overwrite=True`` is.

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
        genome.io.completion.UnfinishedRegistrationError
            If the chimera's directory holds files but no record.
        genome.io.completion.RegistrationMismatchError
            If its record disagrees with what is on disk.
        genome.io.completion.RegistrationError
            If the built ``chrom.sizes`` is not the concatenation the components
            predict. Nothing vouches for the directory in that case — the record is
            written after this check, never before.
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
            return registered
        work = self._work_dir
        work.mkdir(parents=True, exist_ok=True)
        staged = work / f"{self.assembly}.fa"
        digest = self._concatenate(staged)
        files = prepare_fasta(self._place_fasta(staged), overwrite=overwrite)
        self._check_built_names(files)
        # Where the merged annotation lands: every component is here, prepared, with its
        # own default annotation, and the chimera's own chrom.sizes now exists to check a
        # merged GTF against. It belongs above the record, whose details describe it.
        self._record_completion(files, source_url=None, sha256=digest, details=self._details())
        return files

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

    def _details(self) -> dict[str, Any]:
        """Return what this build records about itself beyond the files it wrote.

        The separator its chromosome names carry, and one entry per component. Kept to
        facts a later pass cannot re-derive from the name alone: which spelling was used,
        and what each component was when this chimera was built.
        """
        return {
            _SEPARATOR_KEY: self.separator,
            _COMPONENTS_KEY: [self._component_details(c) for c in self.components],
        }

    def _component_details(self, component: Genome) -> dict[str, Any]:
        """Return one component's entry: its name, and the digest its own record pins.

        Record to record — the component's completion record already holds the sha256 of
        its FASTA, so nothing is rehashed here. ``None`` when that component's record
        pins none, which reads as unknown rather than as wrong.
        """
        # The component's own **Assembly dir**, which is where its FASTA sits. Asked of
        # the component rather than derived from its name, since a component may be
        # registered somewhere other than the shared data root.
        record = read_record(component.fasta_path.parent)
        return {
            _COMPONENT_NAME_KEY: component.assembly,
            _COMPONENT_DIGEST_KEY: None if record is None else record.sha256,
        }


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
    """Return ``line`` with its first whitespace-delimited token suffixed by ``component``.

    ``b">I some description"`` becomes ``b">I__tinyCe some description"``. Everything
    after the first token — the description, and the line ending, whatever it is — is
    written back byte for byte, because it is not what any tool reads: STAR and chromap
    both truncate a header at the first whitespace.
    """
    match = _HEADER_RE.match(line)
    if match is None:  # pragma: no cover - the caller only passes lines starting with '>'
        return line
    name, remainder = match.groups()
    spelled = suffixed(name.decode("utf-8", "surrogateescape"), component, separator)
    return b">" + spelled.encode("utf-8", "surrogateescape") + remainder


def _first_difference(built: list[tuple[str, int]], expected: list[tuple[str, int]]) -> str:
    """Describe how two ``(name, length)`` lists disagree, in one line."""
    for index, (found, wanted) in enumerate(zip(built, expected, strict=False)):
        if found != wanted:
            return (
                f"sequence {index + 1} is {found[0]!r} at {found[1]} bp where "
                f"{wanted[0]!r} at {wanted[1]} bp was expected"
            )
    return f"it carries {len(built)} sequences where {len(expected)} were expected"
