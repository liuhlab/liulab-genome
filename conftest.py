"""Wiring, and nothing else: the suite's two guards, declared where they reach ``src/``.

A conftest's fixtures reach the directory it sits in and everything under it, and nothing
above — so ``tests/conftest.py`` reaches tests and only tests. The package's docstring
examples are collected from ``src/`` (``--doctest-modules``, see ``pyproject.toml``), which
is outside that tree, and a doctest item is a test: it runs behind both autouse guards or
it is the hole in them. This file is the one conftest above both trees, so declaring the
guards here is what gives them the reach the promise already claimed. The guards themselves
stay in :mod:`tests._guards`, which is also where the rest of the suite reaches for
``install_network_guard``.

``pytester`` is enabled alongside them because a plugin list is declared in the root
conftest or nowhere: ``tests/test_doctest_guards.py`` runs a throwaway module in a pytest
of its own to prove the reach rather than assert the shape of this file.
"""

pytest_plugins = ["pytester", "tests._guards"]
