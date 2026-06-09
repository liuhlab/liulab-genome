"""Discovery and version queries for required native tools.

This module is the I/O boundary for shelling out to bundled binaries such as
``samtools`` and ``bedtools``. Per project policy these tools are managed by
pixi (conda-forge + bioconda) and are expected to be on ``PATH`` when the
package runs inside the project environment.

Examples
--------
>>> from genome.external import doctor
>>> versions = doctor()                      # doctest: +SKIP
>>> sorted(versions)                         # doctest: +SKIP
['bedtools', 'samtools']
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

REQUIRED_TOOLS: tuple[str, ...] = ("samtools", "bedtools")


class ToolNotFoundError(RuntimeError):
    """Raised when a required external tool cannot be located on ``PATH``."""


def _resolve(name: str) -> str:
    """Return the absolute path to ``name`` or raise :class:`ToolNotFoundError`.

    Resolution order:

    1. ``shutil.which(name)`` — the normal ``PATH`` lookup.
    2. The ``bin/`` directory of the running interpreter
       (``Path(sys.executable).parent``). In a conda/pixi environment the native
       tools are installed alongside ``python``, so this still finds them when the
       script is run with the environment's interpreter by absolute path without
       the environment being activated (so ``PATH`` lacks its ``bin/``).

    Raises
    ------
    ToolNotFoundError
        If ``name`` is found by neither lookup.
    """
    path = shutil.which(name)
    if path is not None:
        return path

    sibling = Path(sys.executable).parent / name
    if sibling.is_file() and os.access(sibling, os.X_OK):
        return str(sibling)

    raise ToolNotFoundError(
        f"{name!r} not found on PATH. Activate the project environment with "
        f"`pixi shell` (or run via `pixi run`), or add the tool with "
        f"`pixi add {name}` (channels: conda-forge, bioconda)."
    )


def tool_version(name: str) -> str:
    """Return the first line of ``<name> --version`` output.

    Parameters
    ----------
    name
        The executable to query (e.g. ``"samtools"``).

    Returns
    -------
    str
        The first non-empty line of the tool's ``--version`` output (stdout
        preferred, falling back to stderr — different tools choose differently).

    Raises
    ------
    ToolNotFoundError
        If ``name`` is not on ``PATH``.
    subprocess.CalledProcessError
        If the tool exits non-zero when asked for its version.
    """
    path = _resolve(name)
    result = subprocess.run(
        [path, "--version"],
        capture_output=True,
        text=True,
        check=True,
    )
    text = (result.stdout or result.stderr).strip()
    if not text:
        return ""
    return text.splitlines()[0]


def doctor() -> dict[str, str]:
    """Verify required native tools and return their versions.

    Returns
    -------
    dict[str, str]
        Mapping from tool name to its reported version line. The set of
        required tools is :data:`REQUIRED_TOOLS`.

    Raises
    ------
    ToolNotFoundError
        If any required tool is missing from ``PATH``; the message names the
        missing tool and explains how to fix it.

    Examples
    --------
    >>> doctor()                              # doctest: +SKIP
    {'samtools': 'samtools 1.21 ...', 'bedtools': 'bedtools v2.31.1'}
    """
    return {name: tool_version(name) for name in REQUIRED_TOOLS}
