"""Shared fixtures for the test suite."""

from __future__ import annotations

import os
from collections.abc import Callable, Sequence
from pathlib import Path

import pytest

from genome.io import utils as utils_mod


@pytest.fixture
def touch_newer_than() -> Callable[..., None]:
    """Return a helper that sets ``path``'s mtime ``delta`` seconds after ``reference``'s."""

    def _touch(path: Path, reference: Path, *, delta: float = 10.0) -> None:
        mtime = reference.stat().st_mtime + delta
        os.utime(path, (mtime, mtime))

    return _touch


@pytest.fixture
def run_calls(monkeypatch: pytest.MonkeyPatch) -> list[tuple[str, list[str]]]:
    """Record (and suppress) every ``genome.io.utils._run`` call so caching is observable."""
    calls: list[tuple[str, list[str]]] = []

    def fake_run(name: str, args: Sequence[str]) -> None:
        calls.append((name, list(args)))

    monkeypatch.setattr(utils_mod, "_run", fake_run)
    return calls
