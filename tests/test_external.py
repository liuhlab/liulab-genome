"""Tests for genome.external — gated on the presence of native binaries.

These tests run unconditionally inside the pixi env (samtools/bedtools are
runtime conda deps). Outside that env they skip cleanly.
"""

from __future__ import annotations

import shutil

import pytest

from genome.external import REQUIRED_TOOLS, ToolNotFoundError, doctor, tool_version

_BINARIES_PRESENT = all(shutil.which(t) is not None for t in REQUIRED_TOOLS)


def test_missing_tool_raises_actionable() -> None:
    with pytest.raises(ToolNotFoundError, match="pixi"):
        tool_version("definitely-not-a-real-tool-xyz")


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
