"""The one place this package shells out to a native binary.

An **External tool** is a binary the package drives rather than reimplements: ``samtools``,
``faToTwoBit`` and ``twoBitInfo`` prepare an assembly, STAR and chromap build an index.
:class:`ExternalTool` is the whole of what a caller needs from one — where it is, what
version it is, how to run it, and how to run it only when what it would build is out of
date — and two adapters implement it. :class:`InstalledTool` resolves on ``PATH`` and
shells out; :class:`RecordingTool` records the calls and runs nothing, which is what a
test binds an aligner to instead of a binary it does not have.

A tool that cannot be located raises :class:`ToolNotFoundError` whose message *is* its
install instructions, so no caller has to know which conda package carries which binary —
``faToTwoBit`` comes from ``ucsc-fatotwobit`` and nothing but this module needs to know it.

Side effects live here: nothing outside this module calls :mod:`subprocess`.

Examples
--------
>>> from genome.external import RecordingTool
>>> star = RecordingTool("STAR", version="2.7.11b")
>>> star.run(["--runMode", "genomeGenerate"], capture=False)
''
>>> star.version, star.calls[0].args
('2.7.11b', ('--runMode', 'genomeGenerate'))
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from abc import ABC, abstractmethod
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class _Installation:
    """How one binary is installed, and where its own documentation lives."""

    package: str
    homepage: str | None = None


#: Every binary this package drives, mapped to what installs it. A name absent from here
#: installs as ``pixi add <lowercased name>``, which is right for most; the entries are
#: the ones where it is not, plus the homepages worth quoting at someone who has to go
#: install something. Nothing outside this module spells a conda package name.
_INSTALLATIONS: dict[str, _Installation] = {
    "samtools": _Installation("samtools", "https://www.htslib.org"),
    "faToTwoBit": _Installation("ucsc-fatotwobit"),
    "twoBitInfo": _Installation("ucsc-twobitinfo"),
    "STAR": _Installation("star", "https://github.com/alexdobin/STAR"),
    "chromap": _Installation("chromap", "https://github.com/haowenz/chromap"),
    # Known but not required: nothing here shells out to bedtools today, and the rules
    # name it an **External tool** all the same. Knowing what installs it costs one line
    # and is what makes the first caller that does reach for it fail with a command
    # rather than a bare name.
    "bedtools": _Installation("bedtools", "https://bedtools.readthedocs.io"),
}

#: The **External tool**s every run of this package needs — the ones that prepare an
#: assembly, which is the path nothing can avoid. ``genome doctor`` checks exactly these.
#: STAR and chromap are deliberately absent: they are optional features, and an aligner
#: checks for its own binary when it is asked to build. So is ``bedtools``, which nothing
#: here shells out to yet — :data:`_INSTALLATIONS` still knows what installs it.
REQUIRED_TOOLS: tuple[str, ...] = ("samtools", "faToTwoBit", "twoBitInfo")

#: What :func:`doctor` reports for a tool that runs but will not identify itself —
#: several UCSC binaries reject ``--version`` outright. The tool is listed all the same,
#: because presence is the question ``doctor`` answers and presence is not in doubt.
NO_VERSION_REPORTED = "installed; reports no version"


class ToolNotFoundError(RuntimeError):
    """Raised when an **External tool** cannot be located on ``PATH``.

    The message is the tool's :meth:`ExternalTool.install_instructions`, so the next
    action is in the exception rather than somewhere the caller has to go and look.
    """


@dataclass(frozen=True)
class ToolCall:
    """One invocation of an **External tool** — what it was asked to do, and where.

    Attributes
    ----------
    args : tuple of str
        The arguments after the executable itself.
    cwd : pathlib.Path or None
        The directory the tool ran in, or ``None`` for the caller's own.
    capture : bool
        Whether the tool's output was captured rather than inherited.

    Examples
    --------
    >>> ToolCall(("faidx", "hg38.fa"), None, capture=True).args[0]
    'faidx'
    """

    args: tuple[str, ...]
    cwd: Path | None
    capture: bool


def is_fresh(output: Path, inputs: Sequence[Path]) -> bool:
    """Return whether ``output`` is fresh against ``inputs`` — the **Freshness** rule.

    Fresh means ``output`` exists, is non-empty, and is at least as new as every
    input — the same staleness rule ``make`` uses. Missing inputs are ignored;
    the caller validates that required inputs exist.
    """
    if not output.is_file() or output.stat().st_size == 0:
        return False
    out_mtime = output.stat().st_mtime
    return all(out_mtime >= inp.stat().st_mtime for inp in inputs if inp.is_file())


class ExternalTool(ABC):
    """One binary the package drives instead of reimplementing it.

    Four questions and nothing else: where the binary is (:attr:`path`), what version it
    is (:attr:`version`), how to run it (:meth:`run`), and how to run it only when what
    it would build is stale (:meth:`run_to`). Both are answered lazily and remembered, so
    holding a tool costs nothing and a caller that never runs it never needs it installed.

    Subclass it only to add an adapter: :class:`InstalledTool` shells out for real and
    :class:`RecordingTool` records instead, and between them they are the seam. Everything
    a caller relies on — the freshness rule, the failure message, the install
    instructions — lives here, once, so the two adapters cannot drift.

    Parameters
    ----------
    name : str
        The executable's name, spelled as it is on ``PATH`` — ``"faToTwoBit"``, not
        ``"fatotwobit"``.
    package : str, optional
        The conda package that installs it. Defaults to what :data:`_INSTALLATIONS`
        records, or the lowercased name.
    homepage : str, optional
        The tool's own documentation, quoted in :meth:`install_instructions`.

    Attributes
    ----------
    name : str
        The executable's name.
    package : str
        The conda package that installs it.
    homepage : str or None
        The tool's own documentation, when there is one worth quoting.

    Examples
    --------
    >>> tool = RecordingTool("samtools", version="samtools 1.21")
    >>> tool.run(["faidx", "hg38.fa"])
    ''
    >>> tool.calls[0].args
    ('faidx', 'hg38.fa')
    """

    def __init__(
        self, name: str, *, package: str | None = None, homepage: str | None = None
    ) -> None:
        known = _INSTALLATIONS.get(name)
        self.name = name
        self.package = package or (known.package if known else name.lower())
        self.homepage = homepage or (known.homepage if known else None)
        self._path: str | None = None
        self._version: str | None = None

    # -- what a caller asks of a tool ----------------------------------------

    @property
    def path(self) -> str:
        """Absolute path to the executable, located once and remembered.

        Returns
        -------
        str
            The path the tool will be run from.

        Raises
        ------
        ToolNotFoundError
            If the binary cannot be located. The message is
            :meth:`install_instructions`.

        Examples
        --------
        >>> RecordingTool("samtools").path
        '/fake/samtools'
        """
        if self._path is None:
            self._path = self._locate()
        return self._path

    @property
    def version(self) -> str:
        """The tool's version line, or ``""`` when it will not identify itself.

        Asked once, on first use, and remembered — so recording the version of a tool
        that just ran costs nothing, and constructing something that holds a tool runs
        no subprocess at all.

        An empty string means *the tool ran and declined*: several UCSC binaries reject
        ``--version`` outright. That is a different answer from a tool that is not there,
        which raises.

        Returns
        -------
        str
            The first non-empty line of the tool's ``--version`` output (stdout
            preferred, falling back to stderr — different tools choose differently), or
            ``""``.

        Raises
        ------
        ToolNotFoundError
            If the binary cannot be located.

        Examples
        --------
        >>> RecordingTool("STAR", version="2.7.11b").version
        '2.7.11b'
        """
        if self._version is None:
            self._version = self._detect_version()
        return self._version

    def install_instructions(self) -> str:
        """Return the text to put in front of someone whose binary is missing.

        Names the exact command, because an error a caller cannot act on is a bug rather
        than an error message. This is also the message :class:`ToolNotFoundError`
        carries, so the instructions travel with the failure.

        Returns
        -------
        str
            Several lines: what is missing, what installs it, and what to check if it
            is installed already.

        Examples
        --------
        >>> print(RecordingTool("STAR").install_instructions())
        STAR is not installed. Add it to the project environment with:
            pixi add star            # channels: conda-forge, bioconda
        Already installed? Activate the environment with `pixi shell`, or run via `pixi run`.
        See https://github.com/alexdobin/STAR for details.
        """
        lines = [
            f"{self.name} is not installed. Add it to the project environment with:",
            f"    pixi add {self.package}            # channels: conda-forge, bioconda",
            "Already installed? Activate the environment with `pixi shell`, or run via `pixi run`.",
        ]
        if self.homepage is not None:
            lines.append(f"See {self.homepage} for details.")
        return "\n".join(lines)

    def run(self, args: Sequence[str], *, cwd: Path | None = None, capture: bool = True) -> str:
        """Run the tool with ``args`` and return what it wrote to stdout.

        Parameters
        ----------
        args : sequence of str
            The arguments after the executable.
        cwd : pathlib.Path, optional
            The directory to run in. Defaults to the caller's own — pass one when the
            tool drops files in its working directory and they belong beside its output.
        capture : bool, default True
            Whether to capture the tool's output or let it inherit this process's stdout
            and stderr. The two are for different runs and both are wanted: capturing
            puts the tool's own diagnostics into the error raised on failure, which is
            what a command that finishes in a second wants; inheriting streams progress
            live, which is what an index build that runs for an hour wants, and leaves
            the failure message pointing at output the reader has already seen.

        Returns
        -------
        str
            The captured stdout, or ``""`` when ``capture`` is false.

        Raises
        ------
        ToolNotFoundError
            If the binary cannot be located.
        RuntimeError
            If the tool exits non-zero. The message names the tool, its exit code and
            the arguments, and carries the tool's own output when it was captured.

        Examples
        --------
        >>> RecordingTool("samtools").run(["faidx", "hg38.fa"])
        ''
        """
        completed = self._execute(list(args), cwd=cwd, capture=capture)
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout or "").strip() if capture else ""
            explain = f": {detail}" if detail else "; see its output above for the error"
            raise RuntimeError(
                f"{self.name} failed (exit {completed.returncode}) for args {list(args)!r}{explain}"
            )
        return completed.stdout or ""

    def run_to(
        self,
        args: Sequence[str],
        *,
        output: Path,
        inputs: Sequence[Path],
        overwrite: bool = False,
    ) -> Path:
        """Run the tool to build ``output``, skipping the call when ``output`` is fresh.

        The cached command-running primitive every preparation step is built from. When
        ``output`` is fresh relative to ``inputs`` the tool is not invoked at all and
        ``output`` is returned as it stands, so re-preparing an assembly costs a handful
        of ``stat`` calls rather than a second pass over a genome.

        The **Freshness** rule lives here rather than in the callers: it is one rule, and
        a caller that had to apply it would restate the branch at every step and could
        drift from what ``overwrite`` means. Running goes back out through :meth:`run`
        rather than around it to :meth:`_execute`, so a single
        ``monkeypatch.setattr(ExternalTool, "run", ...)`` catches every invocation this
        package makes, by either adapter — the property
        :func:`genome.io.download.fetch_url` is spelled for, and for the same reason.

        Parameters
        ----------
        args : sequence of str
            The arguments after the executable, written so the tool produces ``output``.
        output : pathlib.Path
            What this call builds, and what its freshness is judged by.
        inputs : sequence of pathlib.Path
            What ``output`` is built from. ``output`` is stale once any of them is newer.
            Inputs that do not exist are ignored; the caller validates the ones it needs.
        overwrite : bool, default False
            Run regardless of freshness.

        Returns
        -------
        pathlib.Path
            ``output``, whether it was rebuilt or reused.

        Raises
        ------
        ToolNotFoundError
            If the binary cannot be located — and only when the tool is actually run, so
            a fresh output is served without the tool being installed at all.
        RuntimeError
            If the tool exits non-zero.

        Examples
        --------
        >>> from pathlib import Path
        >>> tool = RecordingTool("samtools")
        >>> tool.run_to(                                  # doctest: +SKIP
        ...     ["faidx", "hg38.fa"], output=Path("hg38.fa.fai"), inputs=[Path("hg38.fa")]
        ... )
        PosixPath('hg38.fa.fai')
        """
        if overwrite or not is_fresh(output, inputs):
            self.run(args)
        return output

    # -- the seam an adapter fills -------------------------------------------

    @abstractmethod
    def _locate(self) -> str:
        """Return the absolute path to the executable, or raise :class:`ToolNotFoundError`."""

    @abstractmethod
    def _execute(
        self, args: Sequence[str], *, cwd: Path | None, capture: bool
    ) -> subprocess.CompletedProcess[str]:
        """Run the executable with ``args`` and return what it did, however that is done."""

    def _detect_version(self) -> str:
        """Ask the tool for its version, answering ``""`` when it will not say."""
        try:
            completed = self._execute(["--version"], cwd=None, capture=True)
        except (OSError, subprocess.SubprocessError):
            return ""
        if completed.returncode != 0:
            return ""
        text = (completed.stdout or completed.stderr or "").strip()
        return text.splitlines()[0] if text else ""


class InstalledTool(ExternalTool):
    """The **External tool** as it is installed on this machine.

    Resolution is ``shutil.which``, then the ``bin/`` directory of the running
    interpreter: in a conda/pixi environment the native tools are installed alongside
    ``python``, so the second lookup still finds them when a script is run with the
    environment's interpreter by absolute path and ``PATH`` therefore lacks its ``bin/``.

    Examples
    --------
    >>> InstalledTool("samtools").version                 # doctest: +SKIP
    'samtools 1.21 (using htslib 1.21)'
    """

    def _locate(self) -> str:
        """Return the path ``shutil.which`` finds, else the interpreter's own ``bin/``."""
        path = shutil.which(self.name)
        if path is not None:
            return path

        sibling = Path(sys.executable).parent / self.name
        if sibling.is_file() and os.access(sibling, os.X_OK):
            return str(sibling)

        raise ToolNotFoundError(self.install_instructions())

    def _execute(
        self, args: Sequence[str], *, cwd: Path | None, capture: bool
    ) -> subprocess.CompletedProcess[str]:
        """Shell out, letting a non-zero exit come back rather than raise."""
        return subprocess.run(
            [self.path, *args],
            cwd=cwd,
            capture_output=capture,
            text=True,
            check=False,
        )


class RecordingTool(ExternalTool):
    """An **External tool** that records what it was asked to do and runs nothing.

    The stand-in a test binds where a real binary would go, so an index build is
    exercised end to end on a machine that has no aligner installed. It is a full
    :class:`ExternalTool`, not a patched-out method: the freshness rule, the failure
    message and the version cache are the same code the real one runs, and only the
    execution is replaced.

    Parameters
    ----------
    name : str, default "tool"
        The executable's name, as :class:`ExternalTool` takes it.
    version : str, default "0.0-test"
        What :attr:`~ExternalTool.version` reports; ``""`` stands for a tool that runs
        but will not identify itself.
    path : str, optional
        What :attr:`~ExternalTool.path` reports. Defaults to ``/fake/<name>``.
    stdout : str, default ""
        What each captured run returns.
    on_run : callable, optional
        Called with each :class:`ToolCall` as it is made — how a test makes a stand-in
        leave behind the files a real build would have written.
    package, homepage : str, optional
        As :class:`ExternalTool`.

    Attributes
    ----------
    calls : list of ToolCall
        Every call made, in order.
    exit_code : int
        What the next run reports. Set it non-zero to make :meth:`~ExternalTool.run`
        fail exactly as a real tool failing would, which is how a test reaches the
        failure path without a binary that can fail.
    stdout : str
        What each captured run returns.
    on_run : callable or None
        Called with each :class:`ToolCall` as it is made.

    Examples
    --------
    >>> tool = RecordingTool("chromap")
    >>> tool.exit_code = 1
    >>> tool.run(["--build-index"])
    Traceback (most recent call last):
    RuntimeError: chromap failed (exit 1) for args ['--build-index']...
    """

    def __init__(
        self,
        name: str = "tool",
        *,
        version: str = "0.0-test",
        path: str | None = None,
        stdout: str = "",
        on_run: Callable[[ToolCall], None] | None = None,
        package: str | None = None,
        homepage: str | None = None,
    ) -> None:
        super().__init__(name, package=package, homepage=homepage)
        self.calls: list[ToolCall] = []
        self.exit_code = 0
        self.stdout = stdout
        self.on_run = on_run
        self._reported_path = path or f"/fake/{name}"
        self._reported_version = version

    def _locate(self) -> str:
        """Return the stand-in path; nothing is looked up on ``PATH``."""
        return self._reported_path

    def _detect_version(self) -> str:
        """Return the version this stand-in was told to report, running nothing."""
        return self._reported_version

    def _execute(
        self, args: Sequence[str], *, cwd: Path | None, capture: bool
    ) -> subprocess.CompletedProcess[str]:
        """Record the call and hand back the canned result."""
        call = ToolCall(tuple(args), cwd, capture)
        self.calls.append(call)
        if self.on_run is not None:
            self.on_run(call)
        return subprocess.CompletedProcess(
            [self.path, *args],
            self.exit_code,
            stdout=self.stdout if capture else None,
            stderr="" if capture else None,
        )


def doctor() -> dict[str, str]:
    """Verify every **External tool** the package needs and report what each one is.

    Returns
    -------
    dict of str to str
        Each name in :data:`REQUIRED_TOOLS` mapped to its version line, or to
        :data:`NO_VERSION_REPORTED` when the tool is there but will not identify itself.

    Raises
    ------
    ToolNotFoundError
        If any required tool is missing. The message names the missing tool and the
        command that installs it.

    Examples
    --------
    >>> doctor()                                          # doctest: +SKIP
    {'samtools': 'samtools 1.21', 'faToTwoBit': 'installed; reports no version', ...}
    """
    return {name: InstalledTool(name).version or NO_VERSION_REPORTED for name in REQUIRED_TOOLS}
