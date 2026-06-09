"""Tests for genome.external — gated on the presence of native binaries.

These tests run unconditionally inside the pixi env (samtools/bedtools are
runtime conda deps). Outside that env they skip cleanly.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

import genome.external as external
from genome.external import REQUIRED_TOOLS, ToolNotFoundError, _resolve, doctor, tool_version

_BINARIES_PRESENT = all(shutil.which(t) is not None for t in REQUIRED_TOOLS)


def test_missing_tool_raises_actionable() -> None:
    with pytest.raises(ToolNotFoundError, match="pixi"):
        tool_version("definitely-not-a-real-tool-xyz")


def test_resolve_falls_back_to_interpreter_bin(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # Simulate running the env's interpreter without PATH activated: with an empty
    # PATH the normal lookup misses, but the tool sits beside sys.executable (the
    # conda/pixi bin/), so resolution falls back to that directory.
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    tool = bin_dir / "faToTwoBit"
    tool.write_text("#!/bin/sh\n")
    tool.chmod(0o755)

    monkeypatch.setenv("PATH", "")
    monkeypatch.setattr(external.sys, "executable", str(bin_dir / "python"))

    assert _resolve("faToTwoBit") == str(tool)


def test_resolve_raises_when_neither_lookup_finds_tool(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # Empty PATH and an interpreter bin/ that holds no such tool -> hard failure.
    monkeypatch.setenv("PATH", "")
    monkeypatch.setattr(external.sys, "executable", str(tmp_path / "bin" / "python"))

    with pytest.raises(ToolNotFoundError, match="pixi"):
        _resolve("definitely-not-a-real-tool-xyz")


@pytest.mark.skipif(not _BINARIES_PRESENT, reason="samtools/bedtools not on PATH")
def test_samtools_version_returns_nonempty_string() -> None:
    out = tool_version("samtools")
    assert isinstance(out, str)
    assert out.strip() != ""
    assert "samtools" in out.lower()


@pytest.mark.skipif(not _BINARIES_PRESENT, reason="samtools/bedtools not on PATH")
def test_bedtools_version_returns_nonempty_string() -> None:
    out = tool_version("bedtools")
    assert isinstance(out, str)
    assert out.strip() != ""


@pytest.mark.skipif(not _BINARIES_PRESENT, reason="samtools/bedtools not on PATH")
def test_doctor_returns_all_required_tools() -> None:
    result = doctor()
    assert set(result.keys()) == set(REQUIRED_TOOLS)
    for name, ver in result.items():
        assert ver.strip() != "", f"{name} returned empty version"
