"""Smoke test — the package imports and exposes a version string.

Real feature tests arrive in Phase 3 alongside ``genome.seq``.
"""

import genome


def test_package_imports() -> None:
    """The package imports and surfaces a non-empty ``__version__``."""
    assert isinstance(genome.__version__, str)
    assert genome.__version__
