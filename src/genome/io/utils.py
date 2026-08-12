"""Shared I/O helpers: running native tools, checksumming files, caching by freshness.

These back the ``io`` layer's shelling-out to pixi-managed binaries (``samtools``,
``faToTwoBit``, …). Format-specific logic lives in its own module (e.g.
:mod:`genome.io.fasta`); only format-agnostic plumbing belongs here — :func:`sha256_file`
and :class:`ChecksumMismatchError` know about files and hashes and nothing about
assemblies or annotations, so either kind of file is checked the same way.
"""

from __future__ import annotations

import gzip
import hashlib
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from genome.external import _resolve


class ChecksumMismatchError(ValueError):
    """A file's contents disagree with the checksum that was expected of them.

    Raised where a digest computed from bytes on disk is compared against a recorded
    one — the sha256 an assembly's metadata row pins, say. The message names all three
    things a caller needs to act: which file, what was expected, what was actually
    there.

    Parameters
    ----------
    path : pathlib.Path
        The file that was hashed.
    expected : str
        The digest that was expected of it, as a hex string.
    actual : str
        The digest actually computed from it.

    Attributes
    ----------
    path : pathlib.Path
        The file that was hashed.
    expected : str
        The expected digest.
    actual : str
        The computed digest.

    Examples
    --------
    >>> from pathlib import Path
    >>> raise ChecksumMismatchError(Path("hg38.fa"), "1a2b3c", "9f8e7d")
    Traceback (most recent call last):
    genome.io.utils.ChecksumMismatchError: sha256 mismatch for hg38.fa: expected 1a2b3c, got 9f8e7d. ...
    """

    def __init__(self, path: Path, expected: str, actual: str) -> None:
        super().__init__(
            f"sha256 mismatch for {path}: expected {expected}, got {actual}. "
            f"Delete the file and fetch it again; if the source legitimately changed, "
            f"update the recorded checksum instead."
        )
        self.path = path
        self.expected = expected
        self.actual = actual


def sha256_file(path: Path) -> str:
    """Return the sha256 hex digest of ``path``, reading it in chunks.

    Streams the file through :func:`hashlib.file_digest`, so a whole-genome FASTA is
    digested a buffer at a time and never lands in memory.

    Parameters
    ----------
    path : pathlib.Path
        The file to digest.

    Returns
    -------
    str
        The lower-case hex digest, without an algorithm prefix.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.

    Examples
    --------
    >>> from pathlib import Path
    >>> sha256_file(Path("tests/data/tiny.fa"))              # doctest: +SKIP
    '9316629bab14f9298a043f8b92e1e04a573b12d6a367ccc07c8f8040e5a13981'
    """
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _gunzip(src: Path, dest: Path) -> Path:
    """Stream-decompress gzip ``src`` into ``dest`` (chunked; never fully in memory)."""
    with gzip.open(src, "rb") as fin, dest.open("wb") as fout:
        shutil.copyfileobj(fin, fout)
    return dest


def _run(name: str, args: Sequence[str]) -> None:
    """Resolve ``name`` on ``PATH`` (via pixi) and run it with ``args``.

    Raises
    ------
    genome.external.ToolNotFoundError
        If ``name`` is not on ``PATH``.
    RuntimeError
        If the tool exits non-zero; the message includes its stderr.
    """
    executable = _resolve(name)
    try:
        subprocess.run([executable, *args], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as err:
        detail = (err.stderr or err.stdout or "").strip()
        raise RuntimeError(
            f"{name} failed (exit {err.returncode}) for args {list(args)!r}: {detail}"
        ) from err


def _is_fresh(output: Path, inputs: Sequence[Path]) -> bool:
    """Return whether ``output`` is an up-to-date cache built from ``inputs``.

    Fresh means ``output`` exists, is non-empty, and is at least as new as every
    input — the same staleness rule ``make`` uses. Missing inputs are ignored;
    the caller validates that required inputs exist.
    """
    if not output.is_file() or output.stat().st_size == 0:
        return False
    out_mtime = output.stat().st_mtime
    return all(out_mtime >= inp.stat().st_mtime for inp in inputs if inp.is_file())


def _run_to(
    name: str,
    args: Sequence[str],
    output: Path,
    inputs: Sequence[Path],
    *,
    overwrite: bool = False,
) -> Path:
    """Run ``name`` to build ``output``, skipping the call when ``output`` is fresh.

    The cached command-running primitive shared by every preparation step. When
    ``output`` is fresh relative to ``inputs`` (see :func:`_is_fresh`) the tool is
    not invoked and ``output`` is returned as is; pass ``overwrite=True`` to
    regenerate unconditionally. ``args`` must be written so the tool produces
    ``output``. Returns ``output``; raises as :func:`_run`.
    """
    if not overwrite and _is_fresh(output, inputs):
        return output
    _run(name, args)
    return output
