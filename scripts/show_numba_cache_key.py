#!/usr/bin/env python
r"""Print the two things numba keys its on-disk JIT cache on, and nothing else.

Why this exists
---------------
``actions/cache`` and numba disagree about what a cache hit means, and only one of
them says so in the log. ``actions/cache`` reports on the *files*: it restored them,
so it prints success. numba then decides per compiled overload whether those files
apply to *this* process, and when they do not it silently recompiles. The run reads
as a hit and the wall says miss, and the only surviving evidence is ``--durations``.

Two things in the cache index can differ between a run that saved and a run that
restored, and they fail in ways worth telling apart:

*The target description*, ``(triple, cpu_name, cpu_features)``, is the middle
element of every overload's key. The host CPU model is in there, GitHub's runners
are heterogeneous, and which one a job lands on is not ours to choose — so an
unpinned key varies per run. A mismatch loses only the overloads compiled for the
other CPU; the index keeps its old entries and grows new ones beside them.

*The source stamp*, ``(st_mtime, st_size)`` of the ``.py`` numba compiled, is checked
once for the whole file. numba stamps it at full float precision here — the locator
that floors the mtime is a different one, for in-tree caches — so sub-second drift
counts. A mismatch discards the **entire** index rather than one entry, so it costs
every signature at once.

Printing both means a run that pays JIT it should not have can be read directly:
the target line moving is one bug, the stamp line moving is the other, and the
shape of the loss in ``--durations`` says which without a second run.

Run it in the same environment and with the same variables as the test lane, or it
describes a process that is not the one compiling.
"""

from __future__ import annotations

import pkgutil
import sys
from importlib import import_module
from pathlib import Path

PACKAGE = "memelite"
"""The dependency whose JIT cache CI pays for. Its engines are ``@njit(cache=True)``."""


def _target_description() -> tuple[str, str, str]:
    """Return the target tuple numba puts in every overload's cache key.

    Read off the live codegen rather than rebuilt from ``config``, so it is what
    numba will actually key on and not a second guess at how it decides.

    Examples
    --------
    >>> triple, cpu, _features = _target_description()
    >>> isinstance(triple, str) and isinstance(cpu, str)
    True
    """
    from numba.core.registry import cpu_target

    return cpu_target.target_context.codegen().magic_tuple()


def _host_cpu_name() -> str:
    """Return the CPU model numba would have keyed on had nothing pinned it.

    Examples
    --------
    >>> _host_cpu_name() != ""
    True
    """
    from llvmlite import binding

    return binding.get_host_cpu_name()


def _compiled_sources() -> list[Path]:
    """Return the source files ``PACKAGE`` hands numba, deduplicated and sorted.

    Found by importing each submodule and asking the dispatchers it defines which
    file they came from, so a release that moves its engines is followed rather
    than missed.

    Examples
    --------
    >>> all(p.suffix == ".py" for p in _compiled_sources())
    True
    """
    from numba.core.dispatcher import Dispatcher

    package = import_module(PACKAGE)
    files: set[Path] = set()
    for info in pkgutil.iter_modules(package.__path__):
        try:
            module = import_module(f"{PACKAGE}.{info.name}")
        except Exception as exc:
            # Broad on purpose: a submodule this environment cannot import defines no
            # dispatcher either, so it contributes no stamp. A diagnostic that fails the
            # job it is diagnosing would be worse than one that reports a shorter list —
            # but it says which module it skipped, because a short list and a complete
            # one are otherwise the same list.
            print(f"  (skipped {PACKAGE}.{info.name}: {type(exc).__name__}: {exc})")
            continue
        for value in vars(module).values():
            if isinstance(value, Dispatcher):
                files.add(Path(value.py_func.__code__.co_filename))
    return sorted(files)


def main() -> int:
    """Print the target description and every source stamp.

    Exits 0 whatever happens, including when it cannot work out what to print. This
    reports on the cache; it does not gate anything, and a diagnostic that fails the
    lane it is diagnosing turns one confusing run into two.
    """
    try:
        triple, cpu_name, cpu_features = _target_description()
        print(f"target triple   : {triple}")
        print(f"cpu name        : {cpu_name}  (host reports {_host_cpu_name()})")
        print(f"cpu features    : {cpu_features!r}")
    except Exception as exc:
        print(f"target          : unreadable — {type(exc).__name__}: {exc}")

    try:
        sources = _compiled_sources()
    except Exception as exc:
        print(f"source stamps   : unreadable — {type(exc).__name__}: {exc}")
        return 0

    if not sources:
        print(f"source stamps   : none found — {PACKAGE} defines no cached dispatcher")
        return 0
    print("source stamps   : (st_mtime, st_size) per file numba compiles")
    for path in sources:
        try:
            stat = path.stat()
        except OSError as exc:
            print(f"  {path.name:<24} unreadable — {exc}")
            continue
        print(f"  {path.name:<24} ({stat.st_mtime!r}, {stat.st_size})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
