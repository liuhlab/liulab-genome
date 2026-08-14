"""Shared I/O helpers: checksumming files, and decompressing them.

Format-specific logic lives in its own module (e.g. :mod:`genome.io.fasta`); only
format-agnostic plumbing belongs here — :func:`sha256_file` and
:class:`ChecksumMismatchError` know about files and hashes and nothing about assemblies
or annotations, so either kind of file is checked the same way.

Running a native tool is :mod:`genome.external`'s job and not this module's, and there is
no name-addressed form of it here. A caller that needs one holds an
:class:`genome.external.ExternalTool` and asks it — which is what :mod:`genome.io.fasta`
does at each preparation step.
"""

from __future__ import annotations

import gzip
import hashlib
import shutil
from pathlib import Path


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
