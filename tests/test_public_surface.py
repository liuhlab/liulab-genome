"""The public surface guard: an error a caller can be handed, they can also name.

Twelve exception classes were once raised at callers from modules that re-exported
nothing, so catching one meant importing from a path the API reference declares
internal and free to move. Two of the twelve were worse than inconvenient: they were
bases whose subclasses were already public, so ``except ShippedTableError`` — the whole
point of having a base — could not be written at all.

This module is the guard on that not coming back. It asserts the invariant rather than
the list: **every exception class this package defines is reachable from some public**
``__all__``. A thirteenth added tomorrow fails here, which is what a test naming the
twelve could never do.
"""

from __future__ import annotations

import importlib
import inspect
import pkgutil
from pathlib import Path
from types import ModuleType

import genome

PACKAGE = Path(__file__).resolve().parents[1] / "src" / "genome"


def _packages() -> list[str]:
    """Return every importable package under ``genome``, itself included.

    A package is a directory carrying ``__init__.py``, which is exactly what the
    repo's rule means by a public surface: only an ``__init__.py`` re-export and the
    CLI are public, so these directories are the whole of where a public name can live.
    """
    found = ["genome"]
    for init in sorted(PACKAGE.rglob("__init__.py")):
        parts = init.parent.relative_to(PACKAGE).parts
        if parts:
            found.append("genome." + ".".join(parts))
    return found


def _public_names() -> dict[str, list[str]]:
    """Map every name a package re-exports to the packages exporting it."""
    published: dict[str, list[str]] = {}
    for name in _packages():
        module = importlib.import_module(name)
        for exported in getattr(module, "__all__", []):
            published.setdefault(exported, []).append(name)
    return published


def _modules() -> list[ModuleType]:
    """Return every module in the package, so a class cannot hide in an unimported one."""
    loaded = [genome]
    for _, name, _ in pkgutil.walk_packages(genome.__path__, "genome."):
        loaded.append(importlib.import_module(name))
    return loaded


def _defined_exceptions() -> dict[str, str]:
    """Map each exception class this package *defines* to the module defining it.

    ``__module__`` is what makes this the definition site rather than a re-export, so a
    class imported into five modules is still counted once, where it was written.
    """
    defined: dict[str, str] = {}
    for module in _modules():
        for name, obj in vars(module).items():
            if (
                inspect.isclass(obj)
                and issubclass(obj, BaseException)
                and obj.__module__ == module.__name__
            ):
                defined[name] = module.__name__
    return defined


def test_every_error_this_package_raises_can_be_imported_from_a_public_name() -> None:
    """An error a caller can be handed is an error a caller can name in an ``except``.

    The failure names the offenders and the module each was defined in, because the fix
    is always the same shape: re-export it from the ``__init__.py`` of the package that
    owns it, and add it to that package's ``__all__``.
    """
    published = _public_names()
    unreachable = {
        name: where for name, where in _defined_exceptions().items() if name not in published
    }

    assert not unreachable, (
        "These exception classes are raised at callers but can be imported from no public "
        "name, so catching one means importing from a module the API reference declares "
        "internal:\n"
        + "\n".join(f"  {name} — defined in {where}" for name, where in sorted(unreachable.items()))
        + "\nRe-export each from the __init__.py of the package that owns it and add it to "
        "that package's __all__."
    )
