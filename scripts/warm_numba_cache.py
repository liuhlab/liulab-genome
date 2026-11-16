#!/usr/bin/env python
r"""Compile memelite's comparison engine before the tests run, and print how long it took.

Why this exists
---------------
tomtom, the engine behind ``MotifSet.compare``, is ``@njit(cache=True)``: the first call
in an environment compiles it, which takes 7 to 25 s, and every later call loads that
from numba's on-disk cache in well under a second. The test lane caps each test at a few
seconds, so the compile must not land inside a test. When it did, the timeout stopped it
halfway, nothing was cached, and the next test on that worker failed on numba's broken
state.

CI starts cold on every pull request. Its caches are scoped to the pull request, so one
pull request never reads another's, and a fresh ``pixi install`` gives memelite's files
new modification times, which numba treats as a changed source and discards the whole
cache for. So ``pixi run test`` runs this first, and the compile happens once, here, in
one process, before any worker starts.

One call is enough because ``genome.tf.motif.compare`` hands tomtom C-ordered float64
matrices whatever the motifs were built from, so numba compiles one specialisation and
every comparison reuses it. On a warm cache this takes about a second.
"""

from __future__ import annotations

import sys
import time


def main() -> int:
    """Compare two tiny motifs once, which compiles tomtom or loads it from the cache."""
    started = time.perf_counter()
    import numpy as np

    from genome.tf.motif import Motif, MotifSet

    def spelled(motif_id: str, bases: str) -> Motif:
        columns = [[19.0 if base == letter else 1.0 for letter in bases] for base in "ACGT"]
        return Motif(motif_id, "", np.array(columns))

    targets = MotifSet([spelled("first", "ACGTACGTA"), spelled("second", "TTTTTTTTT")])
    targets.compare(spelled("query", "ACGTAC"))
    print(f"memelite's tomtom compiled or loaded in {time.perf_counter() - started:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
