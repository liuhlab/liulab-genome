"""Shared I/O helpers: running native tools and caching by output freshness.

These back the ``io`` layer's shelling-out to pixi-managed binaries (``samtools``,
``faToTwoBit``, …). Format-specific logic lives in its own module (e.g.
:mod:`genome.io.fasta`); only format-agnostic plumbing belongs here.
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from genome.external import _resolve


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
