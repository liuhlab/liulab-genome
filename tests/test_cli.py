"""Tests for the Typer CLI (``genome``)."""

from __future__ import annotations

import json as _json
import shutil
from pathlib import Path

import pytest
from typer.testing import CliRunner

from genome.cli import app
from genome.external import REQUIRED_TOOLS

from .conftest import FakeFetch

runner = CliRunner()
_BINARIES_PRESENT = all(shutil.which(t) is not None for t in REQUIRED_TOOLS)

#: sha256 of the committed ``tiny.fa``, which the fake fetch serves as any assembly.
_TINY_FA_SHA256 = "9316629bab14f9298a043f8b92e1e04a573b12d6a367ccc07c8f8040e5a13981"


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


class TestTableRow:
    """``genome table-row`` — download an assembly and print its finished table row.

    Offline throughout: ``fake_fetch`` serves the committed ``tiny.fa.gz`` in place of
    any download, and ``LIULAB_DATA`` points the assembly directory at a temp dir. hg38
    and sacCer3 are used because the shipped table pins a source URL for both, which
    also skips the network name check.
    """

    @pytest.fixture(autouse=True)
    def _offline(
        self, fake_fetch: FakeFetch, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        fake_fetch.serve("tiny.fa.gz")
        monkeypatch.setenv("LIULAB_DATA", str(tmp_path))

    def test_prints_the_row_to_paste(self) -> None:
        result = runner.invoke(app, ["table-row", "hg38"])

        assert result.exit_code == 0
        assert result.stdout.strip().split("\t") == [
            "hg38",
            "Homo sapiens",
            "hg38",
            "GRCh38",
            "GCF_000001405.40",
            "9606",
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
            _TINY_FA_SHA256,
        ]

    def test_json(self) -> None:
        result = runner.invoke(app, ["table-row", "hg38", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["assembly_name"] == "hg38"
        assert payload["ncbi_taxid"] == 9606
        assert payload["sha256"] == _TINY_FA_SHA256

    def test_a_checksum_mismatch_exits_nonzero_naming_both_values(self) -> None:
        # sacCer3's row pins the real genome's digest; the fixture is a subsample of it.
        result = runner.invoke(app, ["table-row", "sacCer3"])

        assert result.exit_code == 1
        combined = (result.stdout or "") + (result.stderr or "")
        assert "6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3" in combined
        assert _TINY_FA_SHA256 in combined
