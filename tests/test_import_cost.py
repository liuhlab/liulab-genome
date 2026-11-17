"""Importing the package stays cheap: pandas and pyarrow load on first use, not on import.

pandas was the largest part of what ``import genome`` cost, and every entry point paid it,
because importing any submodule runs ``genome/__init__`` first. A caller that only
names an assembly or locates a tool never touches a frame. So the modules that use pandas
or pyarrow import them inside the functions that need them, and only for typing at module
scope.

The check runs in a fresh interpreter, because this test session has long since loaded
both libraries.
"""

from __future__ import annotations

import subprocess
import sys

#: The package, and the modules seqforge imports from it.
_ENTRY_POINTS = ("genome", "genome.external", "genome.assembly.metadata", "genome.assembly.chimera")


def test_importing_the_package_loads_neither_pandas_nor_pyarrow() -> None:
    probe = "\n".join(
        [
            *(f"import {name}" for name in _ENTRY_POINTS),
            "import sys",
            "print(sorted({'pandas', 'pyarrow'} & set(sys.modules)))",
        ]
    )

    loaded = subprocess.run(
        [sys.executable, "-c", probe], capture_output=True, text=True, check=True
    )

    assert loaded.stdout.strip() == "[]"
