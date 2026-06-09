"""Tests for the Typer CLI (``genome``)."""

from __future__ import annotations

import json as _json
import shutil

import pytest
from typer.testing import CliRunner

from genome.cli import app
from genome.external import REQUIRED_TOOLS

runner = CliRunner()
_BINARIES_PRESENT = all(shutil.which(t) is not None for t in REQUIRED_TOOLS)


class TestVersion:
    def test_version_prints_string(self) -> None:
        result = runner.invoke(app, ["version"])
        assert result.exit_code == 0
        assert result.stdout.strip() != ""


class TestRevcomp:
    def test_basic(self) -> None:
        result = runner.invoke(app, ["revcomp", "ATCG"])
        assert result.exit_code == 0
        assert result.stdout.strip() == "CGAT"

    def test_preserves_case(self) -> None:
        result = runner.invoke(app, ["revcomp", "aTcG"])
        assert result.exit_code == 0
        assert result.stdout.strip() == "CgAt"

    def test_json(self) -> None:
        result = runner.invoke(app, ["revcomp", "ATCG", "--json"])
        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload == {"input": "ATCG", "reverse_complement": "CGAT"}

    def test_invalid_input_exits_2(self) -> None:
        result = runner.invoke(app, ["revcomp", "ATCX"])
        assert result.exit_code == 2
        # Typer's CliRunner merges stderr into output by default in newer versions;
        # check either source for the error message.
        combined = (result.stdout or "") + (result.stderr or "")
        assert "error" in combined.lower()
        assert "X" in combined


@pytest.mark.skipif(not _BINARIES_PRESENT, reason="samtools/bedtools not on PATH")
class TestDoctor:
    def test_doctor_text(self) -> None:
        result = runner.invoke(app, ["doctor"])
        assert result.exit_code == 0
        for tool in REQUIRED_TOOLS:
            assert tool in result.stdout

    def test_doctor_json(self) -> None:
        result = runner.invoke(app, ["doctor", "--json"])
        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert set(payload.keys()) == set(REQUIRED_TOOLS)
