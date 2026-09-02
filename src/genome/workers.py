"""How many processes a run may use, and why the answer is never just the core count.

Lab-HPC infrastructure, not a domain module: it reads ``os.environ`` and asks the operating
system what this process may run on, so it belongs to no bounded context and imports nothing
else from this package. Scanning is its first caller, not its subject.

**One, unless asked.** The library default is a single worker, so importing this package
and calling into it never starts a process nobody asked for — under the spawn start method a
process pool re-imports the caller's script, and an unguarded script would then re-execute
itself. A console script is a proper entry point and can default to the whole machine; an
imported library cannot, which is why :data:`DEFAULT_WORKERS` is 1 and the command line
passes ``None`` instead.

**And "the whole machine" is not the core count.** On a cluster the process is inside an
allocation: a two-CPU Slurm job on a sixteen-core node would start fourteen workers too
many, thrash, and be slower than serial. So resolution reads the allocation first
(:data:`SLURM_CPU_VARS`), then the process's own CPU affinity, and only then the machine.

Examples
--------
>>> resolve_workers(1)              # what the library defaults to
1
>>> resolve_workers(4)              # an explicit request is honoured as given
4
>>> import os
>>> os.environ["SLURM_CPUS_PER_TASK"] = "3"
>>> resolve_workers(None)           # inside an allocation, the allocation wins
3
>>> del os.environ["SLURM_CPUS_PER_TASK"]
"""

from __future__ import annotations

import os

#: How many workers the **library** uses when the caller names none. One, on purpose: see
#: this module's own explanation. The command line resolves ``None`` instead.
DEFAULT_WORKERS = 1

#: The Slurm variables naming an allocation, most specific first. ``SLURM_CPUS_PER_TASK``
#: is what ``--cpus-per-task`` sets and is what one task may use; ``SLURM_CPUS_ON_NODE`` is
#: the node's share of the job and is the fallback for an allocation that named neither.
SLURM_CPU_VARS: tuple[str, str] = ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE")


def resolve_workers(workers: int | None = None) -> int:
    """Return how many processes to run with.

    A number is taken as given — a caller who says 2 gets 2, on a laptop or on a login
    node. ``None`` means *work it out*, and working it out is
    :data:`SLURM_CPU_VARS`, then process affinity, then the machine's cores. Never the
    machine's cores alone: that is the answer that ignores a cluster allocation.

    Nothing here starts a process, and nothing here is cached — an allocation is a property
    of the process, and reading it at the call site is what keeps it current.

    Parameters
    ----------
    workers : int, optional
        The count to use, or ``None`` to resolve it from the environment.

    Returns
    -------
    int
        At least 1.

    Raises
    ------
    ValueError
        If ``workers`` is below 1. Zero workers is not a serial scan, it is no scan.

    Examples
    --------
    >>> resolve_workers(2)
    2
    >>> resolve_workers(None) >= 1
    True
    """
    if workers is not None:
        if workers < 1:
            raise ValueError(
                f"workers must be at least 1, got {workers}. Pass 1 to scan serially, or "
                f"None to take the count from the Slurm allocation or the machine."
            )
        return int(workers)
    allocated = _slurm_cpus()
    return allocated if allocated is not None else _machine_cpus()


def _slurm_cpus() -> int | None:
    """Return the CPUs this Slurm allocation grants, or ``None`` outside an allocation."""
    for name in SLURM_CPU_VARS:
        raw = os.environ.get(name, "").strip()
        if not raw:
            continue
        try:
            found = int(raw)
        except ValueError:
            continue  # e.g. SLURM_CPUS_ON_NODE's "2(x3)" form on a heterogeneous job
        if found >= 1:
            return found
    return None


def _machine_cpus() -> int:
    """Return the CPUs this process may actually run on, affinity respected where it can be."""
    # The affinity-aware count, which the platform may decline to give — hence the two
    # fallbacks under it, the first of which exists only on Linux.
    counted = os.process_cpu_count()
    if counted:
        return counted
    affinity = getattr(os, "sched_getaffinity", None)
    if affinity is not None:
        found = len(affinity(0))
        if found:
            return found
    return os.cpu_count() or 1
