"""FASTA processing — faidx indexing, 2bit conversion, and chrom.sizes.

This is an I/O boundary module: it shells out to native binaries managed by
pixi (``samtools``, ``faToTwoBit``, ``twoBitInfo``) rather than reimplementing
their work in Python. Functions take and return :class:`pathlib.Path` and leave
their outputs next to the input unless an explicit destination is given.

The headline entry point is :func:`prepare_fasta`, which runs all three steps:

1. ``samtools faidx``  →  ``<fasta>.fai``        (random-access index)
2. ``faToTwoBit``      →  ``<fasta>.2bit``       (compact 2bit encoding)
3. ``twoBitInfo``      →  ``<fasta>.chrom.sizes`` (``name<TAB>length`` per seq)

Examples
--------
>>> from genome.io.fasta import prepare_fasta
>>> files = prepare_fasta("hg38.fa")             # doctest: +SKIP
>>> files.chrom_sizes.name                       # doctest: +SKIP
'hg38.chrom.sizes'
"""

from __future__ import annotations

import subprocess
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from genome.external import _resolve

# FASTA suffixes we know how to strip when deriving sibling output names.
_FASTA_SUFFIXES: tuple[str, ...] = (
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
    ".fa",
    ".fasta",
    ".fna",
)


@dataclass(frozen=True)
class GenomeFiles:
    r"""A FASTA together with its derived index and companion files.

    Attributes
    ----------
    fasta : pathlib.Path
        The source FASTA file.
    fai : pathlib.Path
        The ``samtools faidx`` index, ``<fasta>.fai``.
    twobit : pathlib.Path
        The 2bit-encoded sequence.
    chrom_sizes : pathlib.Path
        Two-column ``<name>\t<length>`` chromosome sizes file.
    """

    fasta: Path
    fai: Path
    twobit: Path
    chrom_sizes: Path


def _strip_fasta_suffix(path: Path) -> Path:
    """Return ``path`` with a known FASTA suffix removed (no suffix → unchanged)."""
    name = path.name
    for ext in _FASTA_SUFFIXES:
        if name.endswith(ext):
            return path.with_name(name[: -len(ext)])
    return path.with_suffix("")


def _run(name: str, args: Sequence[str]) -> None:
    """Resolve and run a native tool, raising an actionable error on failure.

    Parameters
    ----------
    name
        Tool to run (e.g. ``"samtools"``); resolved on ``PATH`` via pixi.
    args
        Arguments passed after the executable.

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


def _require_file(path: Path) -> None:
    """Raise :class:`FileNotFoundError` if ``path`` is not an existing file."""
    if not path.is_file():
        raise FileNotFoundError(f"input FASTA not found: {path}")


def faidx(fasta_path: str | Path) -> Path:
    """Build a ``samtools faidx`` index for ``fasta_path``.

    Parameters
    ----------
    fasta_path : str or pathlib.Path
        FASTA to index. May be bgzipped (``.fa.gz``); plain gzip is not
        supported by ``samtools faidx`` — use bgzip for compressed input.

    Returns
    -------
    pathlib.Path
        Path to the written ``<fasta>.fai`` index.

    Raises
    ------
    FileNotFoundError
        If ``fasta_path`` does not exist.
    RuntimeError
        If ``samtools faidx`` fails.
    """
    fasta = Path(fasta_path)
    _require_file(fasta)
    _run("samtools", ["faidx", str(fasta)])
    return fasta.with_name(fasta.name + ".fai")


def fasta_to_2bit(fasta_path: str | Path, twobit_path: str | Path | None = None) -> Path:
    """Convert a FASTA to UCSC 2bit format via ``faToTwoBit``.

    Parameters
    ----------
    fasta_path : str or pathlib.Path
        Source FASTA (``faToTwoBit`` accepts plain or gzipped FASTA).
    twobit_path : str or pathlib.Path, optional
        Destination ``.2bit`` path. Defaults to the FASTA path with its
        FASTA suffix replaced by ``.2bit`` (e.g. ``hg38.fa`` → ``hg38.2bit``).

    Returns
    -------
    pathlib.Path
        Path to the written ``.2bit`` file.

    Raises
    ------
    FileNotFoundError
        If ``fasta_path`` does not exist.
    RuntimeError
        If ``faToTwoBit`` fails.
    """
    fasta = Path(fasta_path)
    _require_file(fasta)
    twobit = (
        Path(twobit_path)
        if twobit_path is not None
        else _strip_fasta_suffix(fasta).with_suffix(".2bit")
    )
    _run("faToTwoBit", [str(fasta), str(twobit)])
    return twobit


def twobit_to_chrom_sizes(twobit_path: str | Path, sizes_path: str | Path | None = None) -> Path:
    """Write a ``chrom.sizes`` file from a 2bit file via ``twoBitInfo``.

    Parameters
    ----------
    twobit_path : str or pathlib.Path
        Source ``.2bit`` file.
    sizes_path : str or pathlib.Path, optional
        Destination path. Defaults to the 2bit path with its suffix replaced
        by ``.chrom.sizes`` (e.g. ``hg38.2bit`` → ``hg38.chrom.sizes``).

    Returns
    -------
    pathlib.Path
        Path to the written ``chrom.sizes`` file — one ``name<TAB>length``
        line per sequence.

    Raises
    ------
    FileNotFoundError
        If ``twobit_path`` does not exist.
    RuntimeError
        If ``twoBitInfo`` fails.
    """
    twobit = Path(twobit_path)
    if not twobit.is_file():
        raise FileNotFoundError(f"2bit file not found: {twobit}")
    sizes = (
        Path(sizes_path)
        if sizes_path is not None
        else twobit.with_name(twobit.stem + ".chrom.sizes")
    )
    _run("twoBitInfo", [str(twobit), str(sizes)])
    return sizes


def prepare_fasta(
    fasta_path: str | Path,
    *,
    twobit_path: str | Path | None = None,
    sizes_path: str | Path | None = None,
) -> GenomeFiles:
    """Index a FASTA, convert it to 2bit, and write its ``chrom.sizes``.

    Convenience wrapper that runs :func:`faidx`, :func:`fasta_to_2bit`, and
    :func:`twobit_to_chrom_sizes` in sequence and collects their outputs.

    Parameters
    ----------
    fasta_path : str or pathlib.Path
        FASTA to process.
    twobit_path : str or pathlib.Path, optional
        Override for the ``.2bit`` destination; see :func:`fasta_to_2bit`.
    sizes_path : str or pathlib.Path, optional
        Override for the ``chrom.sizes`` destination; see
        :func:`twobit_to_chrom_sizes`.

    Returns
    -------
    GenomeFiles
        Paths to the source FASTA and the three generated files.

    Raises
    ------
    FileNotFoundError
        If ``fasta_path`` does not exist.
    RuntimeError
        If any of the underlying native tools fail.

    Examples
    --------
    >>> files = prepare_fasta("genome.fa")        # doctest: +SKIP
    >>> files.fai, files.twobit, files.chrom_sizes        # doctest: +SKIP
    (PosixPath('genome.fa.fai'), PosixPath('genome.2bit'), PosixPath('genome.chrom.sizes'))
    """
    fasta = Path(fasta_path)
    fai = faidx(fasta)
    twobit = fasta_to_2bit(fasta, twobit_path)
    chrom_sizes = twobit_to_chrom_sizes(twobit, sizes_path)
    return GenomeFiles(fasta=fasta, fai=fai, twobit=twobit, chrom_sizes=chrom_sizes)
