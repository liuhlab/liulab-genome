"""Importing the package stays cheap: pandas, pyarrow, MOODS and gffutils load on first use.

pandas was the largest part of what ``import genome`` cost, and every entry point paid it,
because importing any submodule runs ``genome/__init__`` first. A caller that only
names an assembly or locates a tool never touches a frame. So the modules that use pandas
or pyarrow import them inside the functions that need them, and only for typing at module
scope.

MOODS and gffutils are not declared dependencies at all, so for them it is a matter of working
rather than of speed: ``import genome`` has to succeed in an environment that has neither.

The check runs in a fresh interpreter, because this test session has long since loaded
both libraries.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from genome.annotation.database import build_database
from genome.tf.motif.thresholds import cutoffs_for

#: The package, and the modules seqforge imports from it.
_ENTRY_POINTS = (
    "genome",
    "genome.external",
    "genome.assembly.metadata",
    "genome.assembly.chimera",
    "genome.cli",
)

#: Libraries an entry point must not load: heavy ones, and ones the package does not declare.
_LAZY = ("pandas", "pyarrow", "MOODS", "gffutils")


def test_importing_the_package_loads_none_of_the_lazy_libraries() -> None:
    probe = "\n".join(
        [
            *(f"import {name}" for name in _ENTRY_POINTS),
            "import sys",
            f"print(sorted({set(_LAZY)!r} & set(sys.modules)))",
        ]
    )

    loaded = subprocess.run(
        [sys.executable, "-c", probe], capture_output=True, text=True, check=True
    )

    assert loaded.stdout.strip() == "[]"


def test_a_call_that_needs_moods_names_how_to_get_it_when_it_is_missing(
    monkeypatch: pytest.MonkeyPatch, liulab_data: Path
) -> None:
    # None in sys.modules is how an import is told the library is not there.
    monkeypatch.setitem(sys.modules, "MOODS.tools", None)

    # A matrix no earlier run has cached, so the call reaches the library.
    with pytest.raises(ModuleNotFoundError, match="moods") as raised:
        cutoffs_for([[[0.1], [0.2], [0.3], [0.4]]], (0.25,) * 4, 1e-4)

    assert "bioconda" in str(raised.value)
    assert isinstance(raised.value.__cause__, ImportError)


def test_a_call_that_needs_gffutils_names_how_to_get_it_when_it_is_missing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setitem(sys.modules, "gffutils", None)

    with pytest.raises(ModuleNotFoundError, match="gffutils") as raised:
        build_database(tmp_path / "mine.gtf", tmp_path / "mine.db")

    assert isinstance(raised.value.__cause__, ImportError)
