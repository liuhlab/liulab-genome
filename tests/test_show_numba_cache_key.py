"""Tests for scripts/show_numba_cache_key.py — the CI step that makes the cache key readable.

The script exists because numba's cache key is invisible: a run can restore every file
and recompile anyway, and nothing says so. Two things are therefore worth pinning here.

*It reports the key numba will actually use*, not a reconstruction of how numba decides.
The pinned-target test drives a real subprocess with ``NUMBA_CPU_NAME`` set, because
numba resolves that at import and no in-process fixture can move it afterwards — and the
whole point of the CI step is that the pin is observable in the log.

*It cannot fail the lane it diagnoses.* A diagnostic that raises turns one confusing run
into two, so the failure paths are driven directly with the sources unreadable.
"""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "show_numba_cache_key.py"


@pytest.fixture(scope="module")
def script() -> ModuleType:
    """Import the script by path — ``scripts/`` is deliberately not a package."""
    spec = importlib.util.spec_from_file_location("show_numba_cache_key", _SCRIPT)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_sources_it_stamps_are_the_files_numba_compiles(script: ModuleType) -> None:
    sources = script._compiled_sources()
    assert sources, "memelite defines cached dispatchers; finding none means the probe broke"
    assert all(p.suffix == ".py" for p in sources)
    assert all(p.is_file() for p in sources)
    assert sources == sorted(set(sources)), "deduplicated and ordered, so runs are comparable"


def test_the_target_it_reports_is_the_one_numba_keys_on(script: ModuleType) -> None:
    from numba.core.registry import cpu_target

    # Untyped: `target_context` is a custom cached-property descriptor pyright cannot
    # follow. The point of the assertion is that the script reads the live codegen
    # rather than rebuilding numba's decision from `config`.
    context: Any = cpu_target.target_context
    assert script._target_description() == context.codegen().magic_tuple()


def test_it_prints_the_three_parts_of_the_key_and_one_stamp_per_source(
    script: ModuleType, capsys: pytest.CaptureFixture[str]
) -> None:
    assert script.main() == 0
    out = capsys.readouterr().out
    assert "target triple" in out
    assert "cpu name" in out
    assert "cpu features" in out
    for path in script._compiled_sources():
        assert path.name in out


def test_a_source_it_cannot_stat_is_named_rather_than_raised(
    script: ModuleType, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(script, "_compiled_sources", lambda: [Path("/nonexistent/gone.py")])
    assert script.main() == 0
    assert "gone.py" in capsys.readouterr().out


def test_a_package_that_compiles_nothing_is_reported_not_crashed(
    script: ModuleType, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(script, "_compiled_sources", lambda: [])
    assert script.main() == 0
    assert "none found" in capsys.readouterr().out


def test_a_probe_that_cannot_answer_still_leaves_the_lane_green(
    script: ModuleType, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    def _boom() -> None:
        raise RuntimeError("llvm went away")

    monkeypatch.setattr(script, "_target_description", _boom)
    monkeypatch.setattr(script, "_compiled_sources", _boom)
    assert script.main() == 0
    out = capsys.readouterr().out
    assert "unreadable" in out
    assert "llvm went away" in out


def test_the_pin_ci_sets_is_visible_in_what_the_script_prints() -> None:
    """The contract the CI step exists for: pin the target, and the log says so."""
    env = {**os.environ, "NUMBA_CPU_NAME": "generic", "NUMBA_CPU_FEATURES": ""}
    done = subprocess.run(
        [sys.executable, str(_SCRIPT)], env=env, capture_output=True, text=True, check=False
    )
    assert done.returncode == 0, done.stderr
    assert "cpu name        : generic" in done.stdout
    assert "cpu features    : ''" in done.stdout
    # The contrast is the point: an unpinned run would key on this instead.
    assert "host reports" in done.stdout
