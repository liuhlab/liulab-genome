"""Tests for genome.io.utils — running native tools and output-freshness caching.

None of these need the native binaries: ``_run``'s success path is exercised by
the end-to-end tests in test_fasta, and here ``_run`` is either stubbed (for the
caching tests, via the ``run_calls`` fixture) or pointed at the interpreter to
drive its error handling.
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from pathlib import Path

import pytest

from genome.external import ToolNotFoundError
from genome.io import utils as utils_mod
from genome.io.utils import _is_fresh, _run, _run_to


def test_is_fresh_rules(tmp_path: Path, touch_newer_than: Callable[..., None]) -> None:
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"

    assert _is_fresh(out, [src]) is False  # missing output
    out.write_text("")
    assert _is_fresh(out, [src]) is False  # empty output
    out.write_text("y")
    touch_newer_than(out, src)
    assert _is_fresh(out, [src]) is True  # non-empty and newer
    touch_newer_than(src, out)
    assert _is_fresh(out, [src]) is False  # input regenerated -> stale


def test_run_to_runs_when_output_missing(
    tmp_path: Path, run_calls: list[tuple[str, list[str]]]
) -> None:
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"

    result = _run_to("tool", ["build", str(out)], out, [src])

    assert result == out
    assert run_calls == [("tool", ["build", str(out)])]


def test_run_to_skips_when_output_fresh(
    tmp_path: Path,
    run_calls: list[tuple[str, list[str]]],
    touch_newer_than: Callable[..., None],
) -> None:
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"
    out.write_text("cached")
    touch_newer_than(out, src)

    result = _run_to("tool", ["build"], out, [src])

    assert result == out
    assert run_calls == []  # served from the fresh cache


def test_run_to_overwrite_forces_run(
    tmp_path: Path,
    run_calls: list[tuple[str, list[str]]],
    touch_newer_than: Callable[..., None],
) -> None:
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"
    out.write_text("cached")
    touch_newer_than(out, src)

    _run_to("tool", ["build"], out, [src], overwrite=True)

    assert run_calls == [("tool", ["build"])]


def test_run_to_reruns_when_input_is_newer(
    tmp_path: Path,
    run_calls: list[tuple[str, list[str]]],
    touch_newer_than: Callable[..., None],
) -> None:
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"
    out.write_text("cached")
    touch_newer_than(src, out)  # input regenerated after the output

    _run_to("tool", ["build"], out, [src])

    assert run_calls == [("tool", ["build"])]


def test_run_raises_when_tool_missing() -> None:
    with pytest.raises(ToolNotFoundError):
        _run("definitely-not-a-real-tool-xyz", [])


def test_run_wraps_nonzero_exit_in_runtime_error(monkeypatch: pytest.MonkeyPatch) -> None:
    # Resolve to the interpreter and run a snippet that fails with known stderr;
    # _run should surface that stderr in an actionable RuntimeError.
    monkeypatch.setattr(utils_mod, "_resolve", lambda _name: sys.executable)

    with pytest.raises(RuntimeError, match="boom"):
        _run("python", ["-c", "import sys; sys.stderr.write('boom'); sys.exit(1)"])
