"""The shared source reader, tested — because four structural guards now believe it.

`tests/_sources.py` is not a guard itself; it is what three import bans and one
user-facing-text guard read the package through. A silent wrong answer here does not fail
a test, it *passes* four of them, so the answers are pinned rather than assumed.

Each case below is one a guard depends on and would have got wrong on its own: the two
evasions an import ban has to see through, the lookalike it must not fire on, and the
empty walk that would otherwise read as a clean package.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from ._sources import PACKAGE, imported_modules, imports_any, parse, sources


def _module(directory: Path, name: str, body: str) -> Path:
    """Write one throwaway module and return its path."""
    path = directory / name
    path.write_text(body, encoding="utf-8")
    return path


def test_the_package_root_is_the_installed_one() -> None:
    """A guard should ask about what ships, not about the checkout it happens to run in."""
    assert PACKAGE.name == "genome"
    assert (PACKAGE / "__init__.py").is_file()
    assert len(sources()) > 1


def test_an_import_deferred_into_a_function_body_still_counts(tmp_path: Path) -> None:
    """The first evasion. A ban only the laziest evasion defeats is not a guarantee."""
    path = _module(
        tmp_path,
        "deferred.py",
        "def later():\n    from genome.homology import compara\n    return compara\n",
    )

    assert imported_modules(path) == {"genome.homology"}
    assert imports_any(path, ["genome.homology"])


def test_an_import_hidden_under_type_checking_still_counts(tmp_path: Path) -> None:
    """The second evasion, and the one that looks legitimate while it happens."""
    path = _module(
        tmp_path,
        "typed.py",
        "from typing import TYPE_CHECKING\n"
        "if TYPE_CHECKING:\n"
        "    from genome.homology.compara import HomologySet\n",
    )

    assert imports_any(path, ["genome.homology"])


def test_a_module_that_merely_starts_with_the_banned_letters_is_not_banned(
    tmp_path: Path,
) -> None:
    """The dot is the boundary between a subpackage and a different package entirely.

    A bare ``startswith`` bans ``genome.homologyfoo`` too, and a guard that is wrong the
    first time it fires gets deleted rather than fixed.
    """
    path = _module(tmp_path, "lookalike.py", "import genome.homologyfoo\n")

    assert imported_modules(path) == {"genome.homologyfoo"}
    assert not imports_any(path, ["genome.homology"])


def test_a_relative_import_names_nothing_a_ban_can_match(tmp_path: Path) -> None:
    """Resolving one needs the package the file sits in; every ban here is absolute."""
    path = _module(tmp_path, "sibling.py", "from . import neighbour\n")

    assert imported_modules(path) == {""}
    assert not imports_any(path, ["genome"])


def test_a_walk_that_finds_nothing_is_refused_rather_than_passed(tmp_path: Path) -> None:
    """The one failure a structural guard cannot afford: passing by checking nothing."""
    with pytest.raises(AssertionError, match="nothing would be checked"):
        sources(tmp_path)


def test_sources_are_read_as_utf8_whatever_the_machine_prefers(tmp_path: Path) -> None:
    """What this package's own sources are written in is a fact about the repository."""
    path = _module(tmp_path, "accented.py", '"""Écrit en UTF-8, avec une élision."""\n')

    assert "Écrit" in (parse(path).body[0].value.value)  # type: ignore[attr-defined]
