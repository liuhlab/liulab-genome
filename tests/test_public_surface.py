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

**Classes are compared by identity, never by name.** ``VersionedGeneIdError`` is defined
twice — once in :mod:`genome.tf.link` and once in :mod:`genome.homology.compara` — and
they are different classes. A guard keyed on the bare name would let either satisfy the
other, so making one private would still pass. Both are public today; the point is that
the guard would not notice if one stopped being.
"""

from __future__ import annotations

import importlib
import inspect
import pkgutil
from types import ModuleType

import genome
import genome.store


def _every_module() -> list[ModuleType]:
    """Import and return every module in the package, so none hides a class from the walk."""
    modules = [genome]
    for _, name, _ in pkgutil.walk_packages(genome.__path__, "genome."):
        modules.append(importlib.import_module(name))
    return modules


def _packages(modules: list[ModuleType]) -> list[ModuleType]:
    """Narrow ``modules`` to the packages among them — the whole of where a public name can live.

    A package is the thing carrying ``__path__``, which is exactly what the repo's rule
    means by a public surface: only an ``__init__.py`` re-export and the CLI are public.
    """
    return [module for module in modules if hasattr(module, "__path__")]


def _published_exceptions(modules: list[ModuleType]) -> set[type[BaseException]]:
    """Return the exception *classes* every package re-exports, resolved to objects."""
    published: set[type[BaseException]] = set()
    for package in _packages(modules):
        for name in getattr(package, "__all__", []):
            obj = getattr(package, name, None)
            if inspect.isclass(obj) and issubclass(obj, BaseException):
                published.add(obj)
    return published


def _defined_exceptions(modules: list[ModuleType]) -> dict[str, type[BaseException]]:
    """Map ``module.ClassName`` to each exception class this package *defines*.

    ``__module__`` is what makes this the definition site rather than a re-export, so a
    class imported into five modules is counted once, where it was written. The key is
    qualified because two modules may define the same name.
    """
    defined: dict[str, type[BaseException]] = {}
    for module in modules:
        for name, obj in vars(module).items():
            if (
                inspect.isclass(obj)
                and issubclass(obj, BaseException)
                and obj.__module__ == module.__name__
            ):
                defined[f"{module.__name__}.{name}"] = obj
    return defined


def test_every_error_this_package_raises_can_be_imported_from_a_public_name() -> None:
    """An error a caller can be handed is an error a caller can name in an ``except``.

    The failure names the offenders and where each was defined, because the fix is always
    the same shape: re-export it from the ``__init__.py`` of the package that owns it, and
    add it to that package's ``__all__``.
    """
    modules = _every_module()
    published = _published_exceptions(modules)
    unreachable = sorted(
        where for where, error in _defined_exceptions(modules).items() if error not in published
    )

    assert not unreachable, (
        "These exception classes are raised at callers but can be imported from no public "
        "name, so catching one means importing from a module the API reference declares "
        "internal:\n"
        + "\n".join(f"  {where}" for where in unreachable)
        + "\nRe-export each from the __init__.py of the package that owns it and add it to "
        "that package's __all__."
    )


def test_every_name_a_package_publishes_resolves_to_something() -> None:
    """A name in ``__all__`` that resolves to nothing is a broken import waiting to happen.

    It also silently weakens the guard above, which asks whether a class is published: a
    stale string would answer for a class that is no longer there.
    """
    dangling = [
        f"{package.__name__}.{name}"
        for package in _packages(_every_module())
        for name in getattr(package, "__all__", [])
        if not hasattr(package, name)
    ]

    assert not dangling, (
        "These names are listed in an __all__ but resolve to nothing, so `from <package> "
        "import <name>` raises:\n" + "\n".join(f"  {name}" for name in sorted(dangling))
    )


def test_the_store_publishes_exception_classes_and_nothing_else() -> None:
    """The narrowed rule that let :mod:`genome.store` re-export at all, held by a test.

    That package re-exports no callable on purpose: the suite's offline guard patches the
    fetch step on the module object, and a re-exported callable would be a second
    reference ``monkeypatch.setattr`` never reaches. Exception classes are exempt because
    nothing patches one — so the exemption has to stay exactly that narrow, and the
    docstring saying so is not enough on its own.
    """
    not_an_error = sorted(
        name
        for name in genome.store.__all__
        if not (
            inspect.isclass(getattr(genome.store, name, None))
            and issubclass(getattr(genome.store, name), BaseException)
        )
    )

    assert not not_an_error, (
        "genome.store may re-export exception classes and nothing else — these are not "
        "exceptions:\n" + "\n".join(f"  {name}" for name in not_an_error) + "\n"
        "A callable re-exported here is a second reference that monkeypatch.setattr on the "
        "module would never reach, which is the bug the offline guard exists to prevent. "
        "Callers should hold the module: `from genome.store import fetch`."
    )
