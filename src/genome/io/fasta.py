"""FASTA processing — faidx indexing, 2bit conversion, and chrom.sizes.

This is an I/O boundary module: it shells out to native binaries managed by
pixi (``samtools``, ``faToTwoBit``, ``twoBitInfo``) rather than reimplementing
their work in Python. Functions take and return :class:`pathlib.Path` and leave
their outputs next to the input unless an explicit destination is given.

The headline entry point is :func:`prepare_fasta`, which runs all three steps:

1. ``samtools faidx``  →  ``<fasta>.fai``        (random-access index)
2. ``faToTwoBit``      →  ``<fasta>.2bit``       (compact 2bit encoding)
3. ``twoBitInfo``      →  ``<fasta>.chrom.sizes`` (``name<TAB>length`` per seq)

Every step is **cached**: it runs through :func:`_run_to`, which skips the
native tool when its output already exists and is newer than its input (the same
freshness rule ``make`` uses). Re-running a preparation is therefore cheap; pass
``overwrite=True`` to force regeneration.

Examples
--------
>>> from genome.io.fasta import prepare_fasta
>>> files = prepare_fasta("hg38.fa")             # doctest: +SKIP
>>> files.chrom_sizes.name                       # doctest: +SKIP
'hg38.chrom.sizes'
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from genome.io.utils import _run_to

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


def _require_file(path: Path) -> None:
    """Raise :class:`FileNotFoundError` if ``path`` is not an existing file."""
    if not path.is_file():
        raise FileNotFoundError(f"input FASTA not found: {path}")


def faidx(fasta_path: str | Path, *, overwrite: bool = False) -> Path:
    """Build a ``samtools faidx`` index for ``fasta_path``.

    Parameters
    ----------
    fasta_path : str or pathlib.Path
        FASTA to index. May be bgzipped (``.fa.gz``); plain gzip is not
        supported by ``samtools faidx`` — use bgzip for compressed input.
    overwrite : bool, default False
        By default a fresh existing ``<fasta>.fai`` (newer than the FASTA) is
        reused without re-running ``samtools``. Set ``True`` to rebuild it.

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
    fai = fasta.with_name(fasta.name + ".fai")
    return _run_to("samtools", ["faidx", str(fasta)], fai, [fasta], overwrite=overwrite)


def fasta_to_2bit(
    fasta_path: str | Path,
    twobit_path: str | Path | None = None,
    *,
    overwrite: bool = False,
) -> Path:
    """Convert a FASTA to UCSC 2bit format via ``faToTwoBit``.

    Parameters
    ----------
    fasta_path : str or pathlib.Path
        Source FASTA (``faToTwoBit`` accepts plain or gzipped FASTA).
    twobit_path : str or pathlib.Path, optional
        Destination ``.2bit`` path. Defaults to the FASTA path with its
        FASTA suffix replaced by ``.2bit`` (e.g. ``hg38.fa`` → ``hg38.2bit``).
    overwrite : bool, default False
        By default a fresh existing ``.2bit`` (newer than the FASTA) is reused
        without re-running ``faToTwoBit``. Set ``True`` to rebuild it.

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
    return _run_to("faToTwoBit", [str(fasta), str(twobit)], twobit, [fasta], overwrite=overwrite)


def twobit_to_chrom_sizes(
    twobit_path: str | Path,
    sizes_path: str | Path | None = None,
    *,
    overwrite: bool = False,
) -> Path:
    """Write a ``chrom.sizes`` file from a 2bit file via ``twoBitInfo``.

    Parameters
    ----------
    twobit_path : str or pathlib.Path
        Source ``.2bit`` file.
    sizes_path : str or pathlib.Path, optional
        Destination path. Defaults to the 2bit path with its suffix replaced
        by ``.chrom.sizes`` (e.g. ``hg38.2bit`` → ``hg38.chrom.sizes``).
    overwrite : bool, default False
        By default a fresh existing ``chrom.sizes`` (newer than the 2bit) is
        reused without re-running ``twoBitInfo``. Set ``True`` to rebuild it.

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
    return _run_to("twoBitInfo", [str(twobit), str(sizes)], sizes, [twobit], overwrite=overwrite)


def prepare_fasta(
    fasta_path: str | Path,
    *,
    twobit_path: str | Path | None = None,
    sizes_path: str | Path | None = None,
    overwrite: bool = False,
) -> GenomeFiles:
    """Index a FASTA, convert it to 2bit, and write its ``chrom.sizes``.

    Convenience wrapper that runs :func:`faidx`, :func:`fasta_to_2bit`, and
    :func:`twobit_to_chrom_sizes` in sequence and collects their outputs. Each
    step is cached: an output already present and newer than its input is reused
    rather than regenerated, so re-running is cheap (see ``overwrite``).

    Parameters
    ----------
    fasta_path : str or pathlib.Path
        FASTA to process.
    twobit_path : str or pathlib.Path, optional
        Override for the ``.2bit`` destination; see :func:`fasta_to_2bit`.
    sizes_path : str or pathlib.Path, optional
        Override for the ``chrom.sizes`` destination; see
        :func:`twobit_to_chrom_sizes`.
    overwrite : bool, default False
        Force every step to rerun even when its output looks fresh.

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
    fai = faidx(fasta, overwrite=overwrite)
    twobit = fasta_to_2bit(fasta, twobit_path, overwrite=overwrite)
    chrom_sizes = twobit_to_chrom_sizes(twobit, sizes_path, overwrite=overwrite)
    return GenomeFiles(fasta=fasta, fai=fai, twobit=twobit, chrom_sizes=chrom_sizes)
