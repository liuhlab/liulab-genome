"""The record a finished build writes, and the working area it uses until it does.

I/O boundary module: it writes and reads one small JSON file per finished build, and
names the scratch directory a build downloads into. One shape serves every kind of
build — an assembly, an annotation, an aligner index — so a new kind inherits the
finished/unfinished question instead of answering it its own way.

A record is the **only** thing that says a build finished; no caller consults an output
file's mere existence. It is written last, after every file it claims exists, and it is
written atomically, so a record is never observed half-written. Reopening a build
compares the sizes it claims against what is on disk and reads no file contents, so
confirming a prepared genome costs one ``stat`` per claimed file rather than a pass over
a whole genome. Paths are stored relative to the record's own directory, so moving an
assembly does not invalidate it.

Everything a build downloads lands in the working area — :data:`WORK_DIR_NAME` inside
the same directory, hidden and disposable. It is on the same filesystem as the outputs,
so placing an unpacked file is a rename rather than a copy; it survives an interrupted
run, so a repeat costs no second download; and it is removed once the record is written.

Examples
--------
>>> from pathlib import Path
>>> from genome.io import completion
>>> completion.record_path(Path("/data/genome/hg38")).name
'.completion.json'
>>> completion.work_dir(Path("/data/genome/hg38")).name
'.work'
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from genome.external import ToolNotFoundError, tool_version

#: File name of the record a finished build writes — one spelling, everywhere.
RECORD_NAME = ".completion.json"

#: Directory name of the disposable working area a build downloads into, inside the
#: directory it is building. Hidden, holds working state rather than claimed outputs,
#: and is therefore never covered by a record.
WORK_DIR_NAME = ".work"


@dataclass(frozen=True)
class CompletionRecord:
    """What one finished build did, written down where it left its files.

    The provenance of everything in a directory: where the bytes came from, what they
    hashed to, which files were claimed and how big each is, which **External tool**
    versions and which package version produced them, and when it finished. It answers
    "is this finished?" and "how was this made?" with the same file.

    Attributes
    ----------
    kind : str
        What was built — ``"genome"``, ``"annotation"`` or ``"index"``.
    name : str
        The name that build is addressed by: an assembly name, a registered annotation
        name, or an index name.
    files : dict of str to int
        Every file the build claims, as a path **relative to the record's own
        directory** mapped to its size in bytes. Relative so the directory stays
        movable.
    source_url : str or None
        Where the bytes were fetched from, or ``None`` when nothing was fetched.
    sha256 : str or None
        Digest of the unpacked source (ADR-0006), or ``None`` when none was computed.
    tool_versions : dict of str to str
        Version line of each **External tool** the build ran, keyed by tool name. A
        tool that could not be run is simply absent: recording provenance never becomes
        a new hard dependency.
    package_version : str
        The version of this package that wrote the record.
    completed_at : str
        When the build finished, ISO-8601 in UTC.
    details : dict
        Kind-specific extras — the exact command and parameters for an index, say.

    Examples
    --------
    >>> record = CompletionRecord(
    ...     kind="genome",
    ...     name="hg38",
    ...     files={"hg38.fa": 3_099_922_541},
    ...     source_url="https://example.org/hg38.fa.gz",
    ...     sha256="1a2b3c",
    ...     tool_versions={"samtools": "samtools 1.21"},
    ...     package_version="2026.8.0",
    ...     completed_at="2026-08-12T09:00:00+00:00",
    ...     details={},
    ... )
    >>> record.files["hg38.fa"]
    3099922541
    """

    kind: str
    name: str
    files: dict[str, int]
    source_url: str | None
    sha256: str | None
    tool_versions: dict[str, str]
    package_version: str
    completed_at: str
    details: dict[str, Any]


@dataclass(frozen=True)
class FileDisagreement:
    """One file a :class:`CompletionRecord` claims that disk does not bear out.

    Attributes
    ----------
    path : str
        The claimed path, relative to the record's directory.
    expected : int
        The size the record claims, in bytes.
    actual : int or None
        The size on disk, or ``None`` when no file is there at all.

    Examples
    --------
    >>> print(FileDisagreement("hg38.fa", 3_099_922_541, 12))
    hg38.fa: recorded 3099922541 bytes, found 12
    >>> print(FileDisagreement("hg38.2bit", 841_756_144, None))
    hg38.2bit: recorded 841756144 bytes, missing
    """

    path: str
    expected: int
    actual: int | None

    def __str__(self) -> str:
        """Return a one-line description naming the file and how it disagrees."""
        found = "missing" if self.actual is None else f"found {self.actual}"
        return f"{self.path}: recorded {self.expected} bytes, {found}"


def record_path(directory: Path) -> Path:
    """Return the path of ``directory``'s completion record, whether or not it exists.

    Parameters
    ----------
    directory : pathlib.Path
        The directory a build filled.

    Returns
    -------
    pathlib.Path
        ``<directory>/.completion.json``.

    Examples
    --------
    >>> from pathlib import Path
    >>> record_path(Path("/data/genome/hg38"))
    PosixPath('/data/genome/hg38/.completion.json')
    """
    return directory / RECORD_NAME


def work_dir(directory: Path) -> Path:
    """Return ``directory``'s disposable working area, whether or not it exists.

    Parameters
    ----------
    directory : pathlib.Path
        The directory being built.

    Returns
    -------
    pathlib.Path
        ``<directory>/.work`` — see :data:`WORK_DIR_NAME`.

    Examples
    --------
    >>> from pathlib import Path
    >>> work_dir(Path("/data/genome/hg38"))
    PosixPath('/data/genome/hg38/.work')
    """
    return directory / WORK_DIR_NAME


def clear_work_dir(directory: Path) -> None:
    """Remove ``directory``'s working area and everything in it, if it is there.

    Called once the record is written, never before: while a run is still in progress
    the working area holds the archive it downloaded, so an interrupted run repairs
    without fetching it again.

    Parameters
    ----------
    directory : pathlib.Path
        The directory whose working area is discarded.

    Examples
    --------
    >>> from pathlib import Path
    >>> clear_work_dir(Path("/data/genome/hg38"))        # doctest: +SKIP
    """
    shutil.rmtree(work_dir(directory), ignore_errors=True)


def tool_versions(names: Iterable[str]) -> dict[str, str]:
    """Return each named **External tool**'s version line, skipping the ones that fail.

    A tool that is not installed is left out rather than raising: a record is
    provenance, and recording it must never turn a tool into a new hard dependency of
    the build that used it. A tool that is installed but declines to identify itself is
    left out for the same reason — several UCSC binaries reject ``--version`` outright —
    so this reports the versions that could be determined, not the tools that ran.

    Parameters
    ----------
    names : iterable of str
        Executable names, e.g. ``("samtools", "faToTwoBit")``.

    Returns
    -------
    dict of str to str
        Tool name to its reported version line, for every tool that answered.

    Examples
    --------
    >>> tool_versions(["samtools"])                      # doctest: +SKIP
    {'samtools': 'samtools 1.21'}
    >>> tool_versions(["definitelyNotInstalled"])
    {}
    """
    versions: dict[str, str] = {}
    for name in names:
        try:
            versions[name] = tool_version(name)
        except (ToolNotFoundError, OSError, subprocess.SubprocessError):
            continue
    return versions


def build_record(
    directory: Path,
    *,
    kind: str,
    name: str,
    files: Iterable[Path],
    source_url: str | None = None,
    sha256: str | None = None,
    tools: Iterable[str] = (),
    details: dict[str, Any] | None = None,
) -> CompletionRecord:
    """Describe a build that has just finished in ``directory``.

    Stats every file the build claims, asks each named **External tool** for its
    version, and stamps the package version and the current UTC time. Hand the result
    to :func:`write_record`; the two are separate so that what is recorded can be
    inspected before it lands.

    Parameters
    ----------
    directory : pathlib.Path
        The directory the build filled, and the one the record is written in. Every
        path in ``files`` must live inside it.
    kind : str
        What was built — ``"genome"``, ``"annotation"`` or ``"index"``.
    name : str
        The name the build is addressed by.
    files : iterable of pathlib.Path
        Every file the build claims. Each must already exist: the record is written
        last.
    source_url : str, optional
        Where the bytes were fetched from.
    sha256 : str, optional
        Digest of the unpacked source. Pass a digest already computed rather than
        hashing a large file a second time.
    tools : iterable of str, optional
        Names of the **External tools** the build ran, whose versions are recorded.
    details : dict, optional
        Kind-specific extras, such as an index build's command and parameters.

    Returns
    -------
    CompletionRecord
        The record describing this build. Nothing is written.

    Raises
    ------
    FileNotFoundError
        If a claimed file does not exist — the record is written after every file it
        claims, so this means the build is not finished.
    ValueError
        If a claimed file lies outside ``directory``; a record claims only files
        beneath itself, which is what keeps the directory movable.

    Examples
    --------
    >>> from pathlib import Path
    >>> build_record(                                    # doctest: +SKIP
    ...     Path("/data/genome/hg38"),
    ...     kind="genome",
    ...     name="hg38",
    ...     files=[Path("/data/genome/hg38/hg38.fa")],
    ...     tools=["samtools"],
    ... )
    CompletionRecord(kind='genome', name='hg38', ...)
    """
    from genome import __version__  # deferred: the package imports this module

    return CompletionRecord(
        kind=kind,
        name=name,
        files=_claimed_files(directory, files, kind=kind, name=name),
        source_url=source_url,
        sha256=sha256,
        tool_versions=tool_versions(tools),
        package_version=__version__,
        completed_at=datetime.now(UTC).isoformat(timespec="seconds"),
        details=dict(details or {}),
    )


def write_record(directory: Path, record: CompletionRecord) -> Path:
    """Write ``record`` into ``directory``, atomically, and return where it landed.

    The record goes to a temporary file in the same directory and is then renamed over
    its destination, so a reader sees either the previous record or the whole new one —
    never a half-written file, and never one from a run that died mid-write.

    Parameters
    ----------
    directory : pathlib.Path
        The directory the build filled. Created if it does not exist.
    record : CompletionRecord
        What to write, normally from :func:`build_record`.

    Returns
    -------
    pathlib.Path
        The path written — :func:`record_path` of ``directory``.

    Raises
    ------
    OSError
        If the directory cannot be written to.

    Examples
    --------
    >>> write_record(Path("/data/genome/hg38"), record)  # doctest: +SKIP
    PosixPath('/data/genome/hg38/.completion.json')
    """
    directory.mkdir(parents=True, exist_ok=True)
    destination = record_path(directory)
    handle, temporary_name = tempfile.mkstemp(dir=directory, prefix=RECORD_NAME, suffix=".tmp")
    temporary = Path(temporary_name)
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            json.dump(asdict(record), stream, indent=2)
            stream.write("\n")
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def read_record(directory: Path) -> CompletionRecord | None:
    """Return ``directory``'s completion record, or ``None`` when there is not one.

    ``None`` means *not finished*. A record that is absent, unreadable, or not the
    shape this package writes all read the same way here — telling those apart, and
    deciding which of them should raise, belongs to the caller that repairs them.

    Parameters
    ----------
    directory : pathlib.Path
        The directory to read the record from.

    Returns
    -------
    CompletionRecord or None
        The record, or ``None`` when the directory holds no usable one.

    Examples
    --------
    >>> from pathlib import Path
    >>> read_record(Path("/tmp/definitely-not-a-build")) is None
    True
    """
    try:
        payload = json.loads(record_path(directory).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(payload, dict):
        return None
    try:
        return CompletionRecord(**payload)
    except TypeError:
        return None


def disagreements(directory: Path, record: CompletionRecord) -> list[FileDisagreement]:
    """Return every file ``record`` claims that ``directory`` does not bear out.

    The cheap check that answers "is this still finished?": one ``stat`` per claimed
    file, comparing presence and size. **No file contents are read**, so the cost is
    set by how many files a build claims and not by how many bytes they hold — a
    prepared human genome is confirmed in milliseconds.

    An empty list means the directory agrees with its record. A caller that must repair
    a disagreement has, for each entry, which file and how (see
    :class:`FileDisagreement`).

    Parameters
    ----------
    directory : pathlib.Path
        The directory the record was written in; claimed paths are resolved against it.
    record : CompletionRecord
        The record to hold the directory to.

    Returns
    -------
    list of FileDisagreement
        One entry per claimed file that is missing or the wrong size, in the order the
        record claims them. Empty when everything agrees.

    Examples
    --------
    >>> disagreements(Path("/data/genome/hg38"), record)     # doctest: +SKIP
    [FileDisagreement(path='hg38.2bit', expected=841756144, actual=None)]
    """
    found: list[FileDisagreement] = []
    for relative, expected in record.files.items():
        path = directory / relative
        if not path.is_file():
            found.append(FileDisagreement(relative, expected, None))
            continue
        actual = path.stat().st_size
        if actual != expected:
            found.append(FileDisagreement(relative, expected, actual))
    return found


def _claimed_files(
    directory: Path, files: Iterable[Path], *, kind: str, name: str
) -> dict[str, int]:
    """Return each file's size keyed by its path relative to ``directory``."""
    claimed: dict[str, int] = {}
    for path in files:
        try:
            relative = path.relative_to(directory)
        except ValueError:
            raise ValueError(
                f"cannot record {path} for the {kind} {name!r}: a record claims only files "
                f"inside its own directory ({directory}), which is what keeps that directory "
                f"movable. Write the record in the directory the build filled."
            ) from None
        if not path.is_file():
            raise FileNotFoundError(
                f"cannot record a completed {kind} build of {name!r}: {path} does not exist. "
                f"The record is written last, after every file it claims — build the missing "
                f"file first."
            )
        claimed[relative.as_posix()] = path.stat().st_size
    return claimed
