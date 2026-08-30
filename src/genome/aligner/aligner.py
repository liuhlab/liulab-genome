"""General aligner abstraction for building genome indexes.

An :class:`Aligner` binds one **External tool** (STAR, chromap, …) to a
:class:`~genome.assembly.genome.Genome` and knows how to build that aligner's index for it.
The base class owns the cross-aligner plumbing — the on-disk layout at
``<assembly dir>/index/<name>/`` and the completion record that says a build
finished — while each concrete subclass supplies the aligner-specific command and
its exposed parameters via :meth:`Aligner.index`. Locating the binary, asking its
version, running it and saying what installs it are none of this module's business:
they belong to :mod:`genome.external`, which every tool in the package goes through.

An index writes the same record as every other build in this package (see
:mod:`genome.store.completion`), carrying the exact command, every tuning knob that
determined the build, the aligner version and the FASTA consumed, so a directory
can be explained months later. Reading that record is the only thing that ever
answers "is this index finished?"; no caller consults an index file's mere
existence.

That record also pins the **digest of the assembly it was built from**, copied
from the assembly's own record one directory up, so an assembly re-registered
underneath an index stops reading as a finished index. The comparison is record
against record and reads no sequence bytes. Both the layout and that record are
asked of the bound genome's :attr:`~genome.assembly.genome.Genome.assembly_dir`, never
re-derived from the **Data dir**: an index belongs inside the assembly it
indexes, and only the genome knows where it was opened.

Only index construction is implemented here; mapping/alignment is out of scope.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable, Mapping, Sequence
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING, Any

from genome.external import ExternalTool, InstalledTool
from genome.store.completion import (
    RECORD_NAME,
    WORK_DIR_NAME,
    CompletionRecord,
    RegistrationMismatchError,
    build_record,
    check_registration,
    record_path,
    write_record,
)

if TYPE_CHECKING:
    from genome.assembly.genome import Genome

#: ``details`` key under which an index pins the digest of the assembly it was built
#: from. Named for what it covers — the assembly's bytes, not this directory's, which
#: the record's own ``files`` sizes already hold it to — because ``details`` is
#: free-form and is read by a human months later.
_ASSEMBLY_DIGEST_KEY = "assembly_sha256"


class IndexNotBuiltError(RuntimeError):
    """No index has ever been built where one was asked for.

    Distinct from the two states :mod:`genome.store.completion` calls broken: nothing
    on disk is damaged or half-written, there is simply no index there yet. The
    message names the call that builds one, so the gap is self-explaining.
    """


class Aligner(ABC):
    """Base class for an external aligner that can build a genome index.

    Subclasses set the class attributes :attr:`name` (the lowercase identifier
    used in the index path), :attr:`binary` (the executable on ``PATH``) and
    :attr:`_flag_separator` (how its long options join words), and implement
    :meth:`index`, :attr:`_artifact` and :attr:`_build_arguments`. An :meth:`index`
    hands :meth:`_build` a way to compose its command line; :meth:`_build` owns the
    sequence every build shares and composes only when a build is going to run.

    An aligner is *given* its **External tool** rather than making one, and constructing
    it runs nothing at all: the binary is located, and its version asked for, the first
    time either is needed. A binary that is not installed raises
    :class:`~genome.external.ToolNotFoundError` carrying the install instructions, at the
    point the build would have started rather than at construction.

    Parameters
    ----------
    genome : genome.assembly.genome.Genome
        The genome whose reference FASTA will be indexed.
    tool : genome.external.ExternalTool, optional
        The tool to drive. Defaults to :attr:`binary` as installed on this machine; pass
        one to bind the build to a particular executable, or to a recording stand-in.
    """

    #: Lowercase identifier used in the index directory path (e.g. ``"star"``).
    name: str
    #: Executable name expected on ``PATH`` (e.g. ``"STAR"``). What installs it, and
    #: where its own documentation is, belong to :mod:`genome.external` and not here.
    binary: str
    #: What this aligner's long options put between words, which is the whole of how
    #: one tool's flags differ from another's: ``"_"`` leaves a Python keyword alone
    #: (STAR's ``genomeSAindexNbases``), ``"-"`` hyphenates it (chromap's
    #: ``min-frag-length``). Read by :meth:`_flags`, the one renderer.
    _flag_separator: str

    def __init__(self, genome: Genome, *, tool: ExternalTool | None = None) -> None:
        self._genome = genome
        self._tool: ExternalTool = InstalledTool(self.binary) if tool is None else tool

    # -- identity / layout ---------------------------------------------------

    @property
    def assembly(self) -> str:
        """The assembly name of the bound genome."""
        return self._genome.assembly

    @property
    def version(self) -> str:
        """The installed aligner version, asked of the binary on first use."""
        return self._tool.version

    def install_instructions(self) -> str:
        """Return how to install this aligner — its tool's own text.

        Every aligner owes one, and it is the same text the
        :class:`~genome.external.ToolNotFoundError` carries, so a caller reading either
        gets the command to run.

        Returns
        -------
        str
            What to run, naming the bioconda package.

        Examples
        --------
        >>> from genome.aligner.star import STAR
        >>> print(STAR(genome, gtf="gencode_v50").install_instructions())  # doctest: +SKIP
        STAR is not installed. Add it to the project environment with:
            pixi add star            # channels: conda-forge, bioconda
        ...
        """
        return self._tool.install_instructions()

    @property
    def index_dir(self) -> Path:
        """Directory holding this aligner's index for the assembly.

        ``<assembly dir>/index/<name>/`` — inside the **Assembly dir** the bound genome
        was opened in, asked of that genome rather than re-derived from the **Data dir**.
        The two agree under the ordinary layout and part company for a genome opened
        somewhere of its own, where re-deriving would put the index beside a different
        assembly's files and read a completion record that is not there.
        """
        return self._genome.assembly_dir.index_dir(self.name)

    @property
    def index_path(self) -> Path:
        """The built index file or prefix this aligner consumes.

        The exact flavour (a directory, a file, or a path prefix) is decided by
        the subclass via :attr:`_artifact`. Reading this property asserts that
        the build finished, and the completion record in :attr:`index_dir` is
        what says so — never the presence of an index file, which is exactly what
        that record exists to distrust.

        Returns
        -------
        pathlib.Path
            The path to hand to the aligner's own command line.

        Raises
        ------
        IndexNotBuiltError
            If nothing has been built here yet.
        genome.store.completion.UnfinishedRegistrationError
            If the directory holds index files but no record — a build that was
            interrupted before it finished.
        genome.store.completion.RegistrationMismatchError
            If the record disagrees with what is on disk, naming every file that
            differs and how, or if the assembly was re-registered after this index
            was built, naming both digests.
        """
        if self._registration() is None:
            raise IndexNotBuiltError(
                f"no {self.name} index for {self.assembly!r} at {self.index_dir}: nothing "
                f"has been built there yet. Build it with `{self._build_command()}`."
            )
        return self._artifact

    # -- subclass contract ---------------------------------------------------

    @property
    @abstractmethod
    def _artifact(self) -> Path:
        """The index file or prefix this aligner consumes (flavour-specific)."""

    @property
    @abstractmethod
    def _build_arguments(self) -> str:
        """The arguments identifying this index, as written in the build call.

        Rendered into ``Genome.build_<name>_index(<here>)`` by
        :meth:`_build_command`; empty when an index needs no selector.
        """

    @abstractmethod
    def index(self, *, overwrite: bool = False, **kwargs: Any) -> Path:
        """Build the genome index and return :attr:`index_path`."""

    # -- shared helpers for subclasses ---------------------------------------

    def _build(
        self,
        compose: Callable[[], tuple[Sequence[str], dict[str, Any]]],
        *,
        overwrite: bool,
    ) -> Path:
        """Run one index build end to end and return :attr:`index_path`.

        The sequence every aligner shares, so that what an :meth:`index` still writes
        for itself is only what genuinely differs: its command line, its tuning knobs,
        and the sizing it computes on the way.

        A finished index is handed back untouched unless ``overwrite`` forces a
        rebuild, and ``compose`` is reached only past that branch. That is why it is
        a callable rather than a pair already built: composing STAR's command line
        resolves the annotation that named the index directory, and reusing an index
        must not come to depend on that annotation still being registered. Past the
        branch, the earlier record is dropped *before* anything else, so a rebuild
        that dies — in the tool or in ``compose`` — leaves a directory that reads as
        interrupted rather than as finished.

        Parameters
        ----------
        compose : callable
            Called with no arguments, once, and only when a build is going to run.
            Returns the pair ``(args, parameters)``.

            ``args`` is the complete argument list for :attr:`binary` — the flags this
            aligner spells itself, followed by :meth:`_flags` over the caller's extra
            keywords.

            ``parameters`` is every tuning knob that determined this build,
            caller-supplied and package-computed alike, whether or not ``args`` also
            spells it under the tool's own flag name. This is what the completion
            record keeps, and it is never what ``args`` is rendered from: the question
            the field answers is *what was this index built with*, and a knob this
            package computed — STAR's ``genomeSAindexNbases`` and
            ``genomeChrBinNbits`` — is recoverable from a dict key and not from
            grepping a command string.

            Index records written before that contract was settled hold whichever of
            the two meanings their aligner had then, and are not migrated: nothing
            branches on a record's ``details``, which is read only to explain a
            directory to whoever is looking at it.
        overwrite : bool
            Rebuild even over a finished index, and over a directory that cannot be
            trusted rather than raising on it.

        Returns
        -------
        pathlib.Path
            The built index file or prefix — :attr:`index_path`.

        Raises
        ------
        genome.store.completion.RegistrationError
            If the index directory cannot be trusted and ``overwrite`` is false.
        genome.external.ToolNotFoundError
            If the binary is not installed; the message names what installs it.
        RuntimeError
            If the tool exits non-zero.
        """
        if not overwrite and self._registration() is not None:
            return self.index_path

        self._begin_build()

        args, parameters = compose()
        self._run(args)
        self._record_completion(command=[self.binary, *args], parameters=parameters)
        return self.index_path

    def _flags(self, kwargs: Mapping[str, Any]) -> list[str]:
        """Render extra keyword arguments as this aligner's long options.

        The one renderer both aligners use. A keyword becomes ``--`` plus its name with
        underscores replaced by :attr:`_flag_separator`, so ``{"min_frag_length": 30}``
        is ``["--min-frag-length", "30"]`` for chromap and ``["--min_frag_length",
        "30"]`` for STAR. A list or tuple value becomes several arguments after the
        one flag.

        Only the caller's own keywords are rendered. The knobs an :meth:`index` spells
        itself are written into ``args`` under the tool's flag name, and are recorded
        rather than rendered — see :meth:`_build`.
        """
        flags: list[str] = []
        for key, value in kwargs.items():
            flag = f"--{key.replace('_', self._flag_separator)}"
            if isinstance(value, (list, tuple)):
                flags += [flag, *(str(item) for item in value)]
            else:
                flags += [flag, str(value)]
        return flags

    def _build_command(self, *, overwrite: bool = False) -> str:
        """Return the :class:`~genome.assembly.genome.Genome` call that builds this index.

        Quoted verbatim in every error message this class raises: an error a
        caller cannot act on is a bug, not an error message. ``overwrite`` asks
        for the forcing form, which is the repair for a directory that cannot be
        trusted.
        """
        arguments = [argument for argument in (self._build_arguments,) if argument]
        if overwrite:
            arguments.append("overwrite=True")
        return f"Genome.build_{self.name}_index({', '.join(arguments)})"

    def _registration(self) -> CompletionRecord | None:
        """Return this index's completion record, or ``None`` when none was built.

        The three-way answer every build in the package asks for: a record when
        the index is finished and agrees with disk, ``None`` when the directory
        is absent or empty, and a raise when it holds something untrustworthy.

        An index that no longer matches the reference it was built from is one of
        the untrustworthy states, so the record is also held to the assembly's
        current digest (:meth:`_check_assembly_unchanged`) before it is handed
        back — a finished index over a rebuilt reference is not a finished index.

        Raises
        ------
        genome.store.completion.RegistrationError
            If the directory holds files without a record, a record that
            disagrees with what is on disk, or a record built from a different
            reference than the one registered now.
        """
        record = check_registration(self.index_dir, repair=self._build_command(overwrite=True))
        if record is not None:
            self._check_assembly_unchanged(record)
        return record

    def _assembly_digest(self) -> str | None:
        """Return the digest the assembly's own record pins for its FASTA, or ``None``.

        Read out of the completion record one directory up, never computed here:
        holding an index to its reference has to stay a string against a string,
        since hashing a FASTA would make opening an index cost a pass over a whole
        genome — the very cost a record exists to avoid.

        Returns
        -------
        str or None
            The assembly record's ``sha256``, or ``None`` when the assembly has no
            record, or has one that pins no digest. Both mean *unknown*.
        """
        record = self._genome.assembly_dir.read_record()
        return None if record is None else record.sha256

    def _check_assembly_unchanged(self, record: CompletionRecord) -> None:
        """Raise when ``record`` was built from a different reference than is here now.

        The index's own files can be intact while the assembly beneath them was
        re-registered and its sequences, or its chromosome names, changed —
        chromap is the plain case, its index storing no sequence names at all, so
        it stays byte-identical throughout. The two digests are what tells.

        Parameters
        ----------
        record : genome.store.completion.CompletionRecord
            This index's record, already held to the files it claims.

        Raises
        ------
        genome.store.completion.RegistrationMismatchError
            If both sides pin a digest and the two disagree. The message names
            both digests and the call that rebuilds the index.
        """
        built_from = record.details.get(_ASSEMBLY_DIGEST_KEY)
        current = self._assembly_digest()
        # An absent digest on either side means unknown, not wrong: an index built
        # before this was recorded pins none, and an assembly may pin none either.
        # The same reading tool_versions gives a version that could not be gathered.
        if built_from is None or current is None or built_from == current:
            return
        raise RegistrationMismatchError(
            f"the {self.name} index in {self.index_dir} was built from a different "
            f"{self.assembly} than the one registered now: the index pins {built_from}, "
            f"while {self._genome.assembly_dir.record_path} now pins {current}. "
            f"The reference was rebuilt after this index was, so the index no longer "
            f"matches the sequences — or the chromosome names — it would be mapped "
            f"against. Rebuild it with `{self._build_command(overwrite=True)}`."
        )

    def _begin_build(self) -> None:
        """Prepare :attr:`index_dir` for a fresh build, discarding any earlier record.

        Dropping the record before the tool runs is what keeps an interrupted
        rebuild honest: the directory then holds files nothing vouches for, which
        reads as unfinished, rather than an old record whose claims a partial
        rebuild might happen to still satisfy.
        """
        self.index_dir.mkdir(parents=True, exist_ok=True)
        record_path(self.index_dir).unlink(missing_ok=True)

    def _run(self, args: Sequence[str]) -> None:
        """Run the aligner binary with ``args`` inside :attr:`index_dir`.

        The index directory is used as the working directory so any log files
        the tool drops in the CWD stay co-located with the index. The tool's
        stdout/stderr are inherited, so its progress and any error messages
        stream live to the console rather than being captured — an index build
        runs for an hour and a caller watching it wants to see it move.

        Raises
        ------
        genome.external.ToolNotFoundError
            If the binary is not installed. The message names the command that
            installs it.
        RuntimeError
            If the tool exits non-zero. The tool's own output (already printed
            above) carries the diagnostic detail.
        """
        self._tool.run(args, cwd=self.index_dir, capture=False)

    def _claimed_files(self) -> list[Path]:
        """Return every file in :attr:`index_dir` the finished build is answerable for.

        An index is claimed whole. A STAR *genomeDir* is a dozen binary files that
        are worthless one without another, and the tool runs with the index
        directory as its CWD, so a log it drops there is part of what this build
        produced too — recording all of them makes reopening an index a real
        integrity check rather than a check of one representative file.

        Two things are left out: the record itself, which does not exist yet when
        this is called and would otherwise have to describe its own size, and the
        working area, which holds working state rather than claimed outputs.
        """
        return sorted(
            path
            for path in self.index_dir.rglob("*")
            if path.is_file()
            and path.name != RECORD_NAME
            and WORK_DIR_NAME not in path.relative_to(self.index_dir).parts
        )

    def _record_completion(self, *, command: Sequence[str], parameters: dict[str, Any]) -> None:
        """Write the completion record that makes this index finished.

        Called last, after the aligner has written every file it is going to:
        the record is the only thing that says a build finished, so writing it
        any earlier would bless a half-built index. It carries the exact command,
        the knobs the build was determined by (the contract is :meth:`_build`'s),
        the aligner version and the FASTA consumed, so the build can be explained
        and reproduced later.

        Two things say which reference was indexed, and they answer different
        questions: ``fasta`` is the path — *which file* — and is provenance, while
        the assembly's digest is *which bytes*, and is what a later reopening is
        held to. A digest that could not be read is left out rather than written
        as null, so an unknown fact is absent here exactly as it is in
        ``tool_versions``.
        """
        details: dict[str, Any] = {
            "aligner": self.name,
            "binary": self.binary,
            "assembly": self.assembly,
            "fasta": str(self._genome.files.fasta),
            "command": list(command),
            "parameters": parameters,
        }
        digest = self._assembly_digest()
        if digest is not None:
            details[_ASSEMBLY_DIGEST_KEY] = digest
        record = build_record(
            self.index_dir,
            kind="index",
            name=self.index_dir.name,
            files=self._claimed_files(),
            details=details,
        )
        # The tool remembers the version it reported, and the binary that just ran is
        # that same one — asking again would run it a second time for no new information.
        write_record(self.index_dir, replace(record, tool_versions={self.binary: self.version}))
