"""Tests for the Typer CLI (``genome``)."""

from __future__ import annotations

import json as _json
import shutil
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

from genome import metadata
from genome.cli import app
from genome.external import REQUIRED_TOOLS
from genome.io import download as download_mod
from genome.io.fasta import GenomeFiles

from .conftest import FakeFetch

runner = CliRunner()
_BINARIES_PRESENT = all(shutil.which(t) is not None for t in REQUIRED_TOOLS)

#: sha256 of the committed ``tiny.fa``, which the fake fetch serves as any assembly.
_TINY_FA_SHA256 = "9316629bab14f9298a043f8b92e1e04a573b12d6a367ccc07c8f8040e5a13981"

#: sha256 of the committed ``tiny.gtf`` — the unpacked bytes ``tiny.gtf.gz`` yields.
_TINY_GTF_SHA256 = "255f43bd9abef76424d1c2d89a40cccc1a36215409bbc8f32dcead49ca3baf5e"

#: The URL the stood-in annotation row pins, served from ``tests/data``.
_ANNOTATION_URL = "https://mirror.example.invalid/annotations/tiny.gtf.gz"


def _output(result: object) -> str:
    """Return a result's stdout and stderr together, wherever the runner put them."""
    return (getattr(result, "stdout", "") or "") + (getattr(result, "stderr", "") or "")


@dataclass
class _OkResponse:
    """Stand-in for the ``HEAD`` response the golden-path name check reads."""

    status_code: int = 200

    def raise_for_status(self) -> None:
        """Succeed, as a 200 does."""


@pytest.fixture
def offline_prepare(monkeypatch: pytest.MonkeyPatch) -> None:
    """Stub the native preparation and the UCSC name check, so a CLI run needs neither.

    The three derived files are written as ``prepare_fasta`` names them, which is all
    the completion record claims of them.
    """

    def fake_prepare_fasta(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        fasta = Path(fasta_path)
        files = GenomeFiles(
            fasta=fasta,
            fai=fasta.with_name(fasta.name + ".fai"),
            twobit=fasta.with_name(fasta.stem + ".2bit"),
            chrom_sizes=fasta.with_name(fasta.stem + ".chrom.sizes"),
        )
        for derived in (files.fai, files.twobit, files.chrom_sizes):
            derived.write_text(f"derived from {fasta.name}\n")
        return files

    monkeypatch.setattr(download_mod, "prepare_fasta", fake_prepare_fasta)
    monkeypatch.setattr(download_mod.requests, "head", lambda url, **kwargs: _OkResponse())


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


class TestRegister:
    """``genome register`` — prepare an assembly and say what landed.

    Offline throughout: ``fake_fetch`` serves the committed ``tiny.fa.gz`` in place of
    any download, the ``HEAD`` name check is stubbed, and ``LIULAB_DATA`` points the
    assembly directory at a temp dir. The assembly is ``tiny``, which no shipped row
    lists, so nothing is pinned for the fixture to disagree with.
    """

    @pytest.fixture(autouse=True)
    def _offline(
        self,
        fake_fetch: FakeFetch,
        offline_prepare: None,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        fake_fetch.serve("tiny.fa.gz")
        monkeypatch.setenv("LIULAB_DATA", str(tmp_path))

    def test_registers_and_reports_where_it_landed(self, tmp_path: Path) -> None:
        result = runner.invoke(app, ["register", "tiny"])

        assert result.exit_code == 0
        assert str(tmp_path / "genome" / "tiny") in result.stdout
        assert _TINY_FA_SHA256 in result.stdout
        assert (tmp_path / "genome" / "tiny" / "tiny.fa").is_file()

    def test_json(self, tmp_path: Path) -> None:
        result = runner.invoke(app, ["register", "tiny", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["assembly"] == "tiny"
        assert payload["directory"] == str(tmp_path / "genome" / "tiny")
        assert payload["sha256"] == _TINY_FA_SHA256
        assert sorted(payload["files"]) == [
            "tiny.2bit",
            "tiny.chrom.sizes",
            "tiny.fa",
            "tiny.fa.fai",
        ]

    def test_a_broken_directory_exits_non_zero_naming_the_repair(self, tmp_path: Path) -> None:
        directory = tmp_path / "genome" / "tiny"
        directory.mkdir(parents=True)
        (directory / "tiny.fa").write_text("half a genome\n")

        result = runner.invoke(app, ["register", "tiny"])

        assert result.exit_code == 1
        assert "genome register tiny --force" in _output(result)

    def test_force_repairs_what_the_error_named(self, tmp_path: Path) -> None:
        directory = tmp_path / "genome" / "tiny"
        directory.mkdir(parents=True)
        (directory / "tiny.fa").write_text("half a genome\n")

        result = runner.invoke(app, ["register", "tiny", "--force", "--json"])

        assert result.exit_code == 0
        assert _json.loads(result.stdout)["sha256"] == _TINY_FA_SHA256

    def test_registering_from_a_source_never_asks_ucsc(self, data_dir: Path) -> None:
        result = runner.invoke(
            app, ["register", "tiny", "--source", str(data_dir / "tiny.fa.gz"), "--json"]
        )

        assert result.exit_code == 0
        assert _json.loads(result.stdout)["source_url"] == str(data_dir / "tiny.fa.gz")


class TestVerify:
    """``genome verify`` — re-read a FASTA and check it against the official row."""

    @pytest.fixture(autouse=True)
    def _offline(
        self,
        fake_fetch: FakeFetch,
        offline_prepare: None,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        fake_fetch.serve("tiny.fa.gz")
        monkeypatch.setenv("LIULAB_DATA", str(tmp_path))

    def test_reports_the_digest_of_a_registered_assembly(self, tmp_path: Path) -> None:
        assert runner.invoke(app, ["register", "tiny"]).exit_code == 0

        result = runner.invoke(app, ["verify", "tiny"])

        assert result.exit_code == 0
        assert _TINY_FA_SHA256 in result.stdout

    def test_json(self) -> None:
        assert runner.invoke(app, ["register", "tiny"]).exit_code == 0

        result = runner.invoke(app, ["verify", "tiny", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["sha256"] == _TINY_FA_SHA256
        assert payload["expected"] is None  # no row lists "tiny"
        assert payload["verified"] is False

    def test_a_hand_copied_fasta_is_checkable_against_the_official_row(
        self, data_dir: Path
    ) -> None:
        # sacCer3's row pins the real genome's digest; the fixture is a subsample of it,
        # so this is the mismatch a copy from a bad mirror would produce.
        result = runner.invoke(app, ["verify", "sacCer3", "--fasta", str(data_dir / "tiny.fa")])

        assert result.exit_code == 1
        assert "sha256 mismatch" in _output(result)

    def test_nothing_registered_exits_non_zero_naming_the_command(self) -> None:
        result = runner.invoke(app, ["verify", "tiny"])

        assert result.exit_code == 1
        assert "genome register tiny" in _output(result)


class TestRegisterAnnotation:
    """``genome register-annotation`` — fetch, verify and build an annotation by name.

    The CLI is a thin client over the shipped table and takes no metadata argument by
    design, so the table itself is what is stood in for here: one row pointing at the
    committed ``tiny.gtf.gz`` and pinning its digest. The registration it drives — the
    fetch, the checksum, the real gffutils build, the record — is the shipped code.
    """

    @pytest.fixture(autouse=True)
    def _offline(
        self, fake_fetch: FakeFetch, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        fake_fetch.serve("tiny.gtf.gz")
        monkeypatch.setenv("LIULAB_DATA", str(tmp_path))
        table = pd.DataFrame(
            [
                {
                    "assembly": "tiny",
                    "name": "ensgene_v101",
                    "provider": "UCSC",
                    "version": "ensGene.v101",
                    "url": _ANNOTATION_URL,
                    "sha256": _TINY_GTF_SHA256,
                    "default": "yes",
                }
            ],
            dtype=str,
        )
        monkeypatch.setattr(metadata, "_annotation_table", lambda: table)

    def test_registers_and_reports_where_it_landed(self, tmp_path: Path) -> None:
        result = runner.invoke(app, ["register-annotation", "tiny", "ensgene_v101"])

        directory = tmp_path / "genome" / "tiny" / "gtf" / "ensgene_v101"
        assert result.exit_code == 0
        assert str(directory) in result.stdout
        assert _TINY_GTF_SHA256 in result.stdout
        assert (directory / "ensgene_v101.gtf").is_file()
        assert (directory / "ensgene_v101.db").is_file()

    def test_json(self, tmp_path: Path) -> None:
        result = runner.invoke(app, ["register-annotation", "tiny", "ensgene_v101", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["assembly"] == "tiny"
        assert payload["name"] == "ensgene_v101"
        assert payload["directory"] == str(tmp_path / "genome" / "tiny" / "gtf" / "ensgene_v101")
        assert payload["source_url"] == _ANNOTATION_URL
        assert payload["sha256"] == _TINY_GTF_SHA256
        assert sorted(payload["files"]) == ["ensgene_v101.db", "ensgene_v101.gtf"]

    def test_a_name_no_row_lists_exits_non_zero_saying_what_is_offered(self) -> None:
        result = runner.invoke(app, ["register-annotation", "tiny", "nope"])

        assert result.exit_code == 1
        assert "ensgene_v101" in _output(result)

    def test_a_broken_directory_exits_non_zero_naming_the_repair(self, tmp_path: Path) -> None:
        directory = tmp_path / "genome" / "tiny" / "gtf" / "ensgene_v101"
        directory.mkdir(parents=True)
        (directory / "ensgene_v101.db").write_bytes(b"half a database")

        result = runner.invoke(app, ["register-annotation", "tiny", "ensgene_v101"])

        assert result.exit_code == 1
        assert "genome register-annotation tiny ensgene_v101 --force" in _output(result)

    def test_force_repairs_what_the_error_named(self, tmp_path: Path) -> None:
        directory = tmp_path / "genome" / "tiny" / "gtf" / "ensgene_v101"
        directory.mkdir(parents=True)
        (directory / "ensgene_v101.db").write_bytes(b"half a database")

        result = runner.invoke(
            app, ["register-annotation", "tiny", "ensgene_v101", "--force", "--json"]
        )

        assert result.exit_code == 0
        assert _json.loads(result.stdout)["sha256"] == _TINY_GTF_SHA256


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

    def test_an_existing_pin_is_reported_rather_than_enforced(self) -> None:
        # sacCer3's row already pins the real genome's digest, and the fixture is a
        # subsample of it, so the two disagree. This is the command a maintainer runs
        # precisely when an upstream file has changed and the pin must be regenerated,
        # so it prints what actually arrived instead of refusing. Checking a FASTA
        # against the official row is what verifying an assembly is for.
        result = runner.invoke(app, ["table-row", "sacCer3"])

        assert result.exit_code == 0
        row = result.stdout.strip().split("\t")
        assert row[0] == "sacCer3"
        assert row[-1] == _TINY_FA_SHA256
