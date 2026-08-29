"""The per-motif cutoffs a **Threshold** converts to, and the disk cache that keeps them.

A **Threshold** is one per-position p-value. What the engine wants is a *score* per matrix,
and converting the one into the other is the engine's one slow step — a few seconds for a
full vertebrate release, paid before a single base is scanned. It is also a pure function
of the matrices, the **Background** and the p-value, so paying it twice for the same triple
buys nothing.

So it is cached on disk, in the ``thresholds/`` subtree of the **Data dir**'s ``motif/``
directory — beside the JASPAR files, shared by every project on the machine exactly as they
are. The key is a digest of that triple. The background reaching here has already been
rounded onto the 0.001 grid by
:func:`~genome.tf.motif.background.quantise_background`, and *that* is what lets two peak
sets from one genome land on one entry instead of two: their compositions differ in the
fourth decimal, and the fourth decimal is gone.

**A cache is a speed-up and never a dependency.** An unreadable, truncated or stale entry
is a miss and not an error; a directory that cannot be written to makes the scan slow and
not broken. Nothing here raises.

Examples
--------
>>> import os, tempfile
>>> import numpy as np
>>> from genome.tf.motif.motif import Motif
>>> counts = np.zeros((4, 8))
>>> for column, base in enumerate("GATTACAG"):
...     counts["ACGT".index(base), column] = 100.0
>>> matrix = (Motif("MA9999.1", "Gattacag", counts).log_odds() * np.log(2)).tolist()
>>> with tempfile.TemporaryDirectory() as directory:
...     os.environ["LIULAB_DATA"] = directory
...     first = cutoffs_for([matrix], (0.25, 0.25, 0.25, 0.25), 1e-4)
...     second = cutoffs_for([matrix], (0.25, 0.25, 0.25, 0.25), 1e-4)   # off the disk
...     cached = len(list(threshold_cache_dir().glob("*.json")))
...     _ = os.environ.pop("LIULAB_DATA")
>>> first == second, cached
(True, 1)
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any

import MOODS.tools
import numpy as np

from genome.io.registration import motif_data_dir

#: Where the cache lives under the **Data dir**'s ``motif/`` tree. Its own directory rather
#: than a file, so an entry can be deleted one at a time and the whole cache with an ``rm``
#: of one path.
CACHE_SUBDIR = "thresholds"

#: Bumped whenever what an entry means changes, so yesterday's files miss rather than lie.
_CACHE_VERSION = 1


def threshold_cache_dir() -> Path:
    """Return the directory the per-motif cutoffs are cached in.

    ``<liulab_data>/motif/thresholds/``. Nothing is created by asking — the write creates
    what it needs. Delete the directory to force every scan on this machine to recompute.

    Returns
    -------
    pathlib.Path
        ``<liulab_data>/motif/thresholds``.

    Examples
    --------
    >>> import os
    >>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
    >>> threshold_cache_dir()
    PosixPath('/scratch/liulab/motif/thresholds')
    >>> del os.environ["LIULAB_DATA"]
    """
    return motif_data_dir() / CACHE_SUBDIR


def cutoffs_for(
    matrices: list[list[list[float]]],
    background: tuple[float, ...],
    threshold: float,
) -> list[float]:
    """Return the score each matrix must clear, from the cache or from the engine.

    One cutoff per entry of ``matrices``, in the same order — so a doubled list of forward
    matrices and their reverse complements gets a doubled list of cutoffs, and the index
    split that recovers **Strand** still works.

    Parameters
    ----------
    matrices : list of 4 x L float matrices
        Log-odds matrices **in nats**, exactly as the engine will be handed them.
    background : tuple of float
        The **Background** the cutoffs are computed against, already on the 0.001 grid.
    threshold : float
        The **Threshold**: one per-position p-value, in ``(0, 1)``.

    Returns
    -------
    list of float
        One cutoff per matrix, in nats.

    Examples
    --------
    >>> import numpy as np
    >>> from genome.tf.motif.motif import Motif
    >>> counts = np.zeros((4, 8))
    >>> for column, base in enumerate("GATTACAG"):
    ...     counts["ACGT".index(base), column] = 100.0
    >>> matrix = (Motif("MA9999.1", "Gattacag", counts).log_odds() * np.log(2)).tolist()
    >>> [round(cutoff, 3) for cutoff in cutoffs_for([matrix], (0.25,) * 4, 1e-4)]
    [0.492]
    """
    path = threshold_cache_dir() / f"{_digest(matrices, background, threshold)}.json"
    cached = _read(path, len(matrices))
    if cached is not None:
        return cached
    weights = list(background)
    cutoffs = [
        float(MOODS.tools.threshold_from_p(matrix, weights, threshold)) for matrix in matrices
    ]
    _write(
        path,
        {
            "version": _CACHE_VERSION,
            "background": list(background),
            "threshold": threshold,
            "cutoffs": cutoffs,
        },
    )
    return cutoffs


def _digest(
    matrices: list[list[list[float]]],
    background: tuple[float, ...],
    threshold: float,
) -> str:
    """Hash the triple a cutoff is a pure function of: the matrices, the background, the p."""
    running = hashlib.sha256()
    running.update(f"v{_CACHE_VERSION}|{list(background)}|{threshold!r}|{len(matrices)}".encode())
    for matrix in matrices:
        block = np.asarray(matrix, dtype=np.float64)
        running.update(f"|{block.shape}|".encode())
        running.update(block.tobytes())
    return running.hexdigest()


def _read(path: Path, count: int) -> list[float] | None:
    """Read ``count`` cutoffs back, or ``None`` for anything that is not exactly that."""
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
        if record["version"] != _CACHE_VERSION:
            return None
        cutoffs = [float(value) for value in record["cutoffs"]]
    except (OSError, ValueError, TypeError, KeyError):
        return None
    return cutoffs if len(cutoffs) == count else None


def _write(path: Path, record: dict[str, Any]) -> None:
    """Write one entry, atomically, and give up quietly if the cache cannot be written."""
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        # Written beside the entry and renamed onto it: two scans racing on one key then
        # leave a whole file rather than a half-written one for the next reader.
        partial = path.with_name(f"{path.name}.{os.getpid()}.partial")
        partial.write_text(json.dumps(record), encoding="utf-8")
        partial.replace(path)
    except OSError:
        return
