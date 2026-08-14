"""Tests for the Typer CLI (``genome``)."""

from __future__ import annotations

import json as _json
import shutil
from dataclasses import dataclass
from pathlib import Path

import gffutils
import pandas as pd
import pytest
from typer.testing import CliRunner

from genome import metadata
from genome.cli import app
from genome.external import REQUIRED_TOOLS
from genome.io import download as download_mod
from genome.io.completion import record_path
from genome.io.fasta import GenomeFiles
from genome.io.gtf import annotation_dir, register_gtf

from .conftest import FakeFetch

runner = CliRunner()
_BINARIES_PRESENT = all(shutil.which(t) is not None for t in REQUIRED_TOOLS)

#: sha256 of the committed ``tiny.fa``, which the fake fetch serves as any assembly.
_TINY_FA_SHA256 = "9316629bab14f9298a043f8b92e1e04a573b12d6a367ccc07c8f8040e5a13981"

#: sha256 of the committed ``tiny.gtf`` — the unpacked bytes ``tiny.gtf.gz`` yields.
_TINY_GTF_SHA256 = "255f43bd9abef76424d1c2d89a40cccc1a36215409bbc8f32dcead49ca3baf5e"

#: The URL the stood-in annotation row pins, served from ``tests/data``.
_ANNOTATION_URL = "https://mirror.example.invalid/annotations/tiny.gtf.gz"

# A bare exon-level GTF — exon lines and nothing else, which is what gene/transcript
# inference exists for. Built with inference off, its database holds exons alone.
_BARE_GTF = (
    "\n".join(
        [
            'chrI\ttest\texon\t1\t50\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
            'chrI\ttest\texon\t60\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
        ]
    )
    + "\n"
)


def _output(result: object) -> str:
    """Return a result's stdout and stderr together, wherever the runner put them."""
    return (getattr(result, "stdout", "") or "") + (getattr(result, "stderr", "") or "")


def _feature_types(database_path: Path) -> list[str]:
    """The kinds of feature a built database holds, with the connection closed behind us."""
    database = gffutils.FeatureDB(str(database_path))
    try:
        return sorted(database.featuretypes())
    finally:
        database.conn.close()


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
        # …and the command that registers a GTF the table does not list, since that is
        # what a caller who named an unlisted annotation is most likely reaching for.
        assert "genome register-gtf tiny" in _output(result)

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

    def test_the_chromosome_check_is_stood_down_from_the_command_line(
        self, fake_fetch: FakeFetch, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # The committed Ensembl-spelled GTF (I, II, III) against a UCSC-spelled
        # assembly (chrI, chrII, chrIII): refused by default, and registered anyway
        # once the caller says they have looked at the mismatch and accept it.
        fake_fetch.serve("ensembl_style.gtf")
        assembly_dir = tmp_path / "genome" / "tiny"
        assembly_dir.mkdir(parents=True)
        (assembly_dir / "tiny.chrom.sizes").write_text("chrI\t10000\nchrII\t10000\nchrIII\t10000\n")
        table = pd.DataFrame(
            [
                {
                    "assembly": "tiny",
                    "name": "ensgene_v101",
                    "provider": "UCSC",
                    "version": "ensGene.v101",
                    "url": "https://mirror.example.invalid/annotations/ensembl_style.gtf",
                    "sha256": None,
                    "default": "yes",
                }
            ],
            dtype=str,
        )
        monkeypatch.setattr(metadata, "_annotation_table", lambda: table)

        refused = runner.invoke(app, ["register-annotation", "tiny", "ensgene_v101"])

        assert refused.exit_code == 1
        assert "chromosome" in _output(refused)

        result = runner.invoke(
            app,
            ["register-annotation", "tiny", "ensgene_v101", "--no-check-chromosomes", "--json"],
        )

        assert result.exit_code == 0
        details = _json.loads(result.stdout)["details"]
        assert details["chromosomes_checked"] is False
        # The reason rides in `details` as the record holds it — no second spelling of it
        # for the JSON surface to drift from.
        assert details["chromosomes_unchecked_because"] == "caller-override"

    def test_feature_inference_is_reachable_for_a_listed_annotation_too(
        self, fake_fetch: FakeFetch, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # Nothing says a listed annotation declares genes and transcripts, so the
        # inference the API exposes has to be reachable on this command as well.
        bare = tmp_path / "bare.gtf"
        bare.write_text(_BARE_GTF)
        fake_fetch.serve(bare)
        table = pd.DataFrame(
            [
                {
                    "assembly": "tiny",
                    "name": "bare",
                    "provider": "somebody",
                    "version": "1",
                    "url": "https://mirror.example.invalid/annotations/bare.gtf",
                    "sha256": None,
                    "default": "yes",
                }
            ],
            dtype=str,
        )
        monkeypatch.setattr(metadata, "_annotation_table", lambda: table)

        result = runner.invoke(
            app, ["register-annotation", "tiny", "bare", "--infer-genes", "--infer-transcripts"]
        )

        assert result.exit_code == 0
        database = tmp_path / "genome" / "tiny" / "gtf" / "bare" / "bare.db"
        assert _feature_types(database) == ["exon", "gene", "transcript"]


class TestRegisterGtf:
    """``genome register-gtf`` — register a GTF the annotation table does not list.

    The by-path way in, from a shell: no table row, no download, no checksum to compare
    against — the caller says where the file is. Offline by construction, since the GTF
    is a local one; ``LIULAB_DATA`` points the assembly directory at a temp dir.
    """

    @pytest.fixture(autouse=True)
    def _offline(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("LIULAB_DATA", str(tmp_path))

    def test_registers_a_gtf_no_row_lists_and_reports_where_it_landed(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        source = data_dir / "tiny.gtf"

        result = runner.invoke(app, ["register-gtf", "tiny", str(source), "mine"])

        directory = tmp_path / "genome" / "tiny" / "gtf" / "mine"
        assert result.exit_code == 0
        assert str(directory) in result.stdout
        assert str(source) in result.stdout
        assert (directory / "mine.gtf").is_file()
        assert (directory / "mine.db").is_file()

    def test_json(self, tmp_path: Path, data_dir: Path) -> None:
        source = data_dir / "tiny.gtf"

        result = runner.invoke(app, ["register-gtf", "tiny", str(source), "mine", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["assembly"] == "tiny"
        assert payload["name"] == "mine"
        assert payload["directory"] == str(tmp_path / "genome" / "tiny" / "gtf" / "mine")
        assert payload["source_url"] == str(source)
        assert payload["sha256"] == _TINY_GTF_SHA256
        assert sorted(payload["files"]) == ["mine.db", "mine.gtf"]

    def test_it_is_then_listed_as_registered_but_not_offered(self, data_dir: Path) -> None:
        assert (
            runner.invoke(
                app, ["register-gtf", "tiny", str(data_dir / "tiny.gtf"), "mine"]
            ).exit_code
            == 0
        )

        result = runner.invoke(app, ["annotations", "tiny", "--json"])

        assert result.exit_code == 0
        rows = _json.loads(result.stdout)["annotations"]
        assert [(row["name"], row["offered"], row["registered"]) for row in rows] == [
            ("mine", False, True)
        ]

    def test_a_gtf_that_is_not_there_exits_non_zero_saying_what_to_pass(
        self, tmp_path: Path
    ) -> None:
        result = runner.invoke(app, ["register-gtf", "tiny", str(tmp_path / "nope.gtf"), "mine"])

        assert result.exit_code == 1
        assert "GTF file not found" in _output(result)

    def test_a_broken_directory_exits_non_zero_naming_the_repair(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        source = data_dir / "tiny.gtf"
        directory = tmp_path / "genome" / "tiny" / "gtf" / "mine"
        directory.mkdir(parents=True)
        (directory / "mine.db").write_bytes(b"half a database")

        result = runner.invoke(app, ["register-gtf", "tiny", str(source), "mine"])

        assert result.exit_code == 1
        assert f"genome register-gtf tiny {source} mine --force" in _output(result)

    def test_force_repairs_what_the_error_named(self, tmp_path: Path, data_dir: Path) -> None:
        directory = tmp_path / "genome" / "tiny" / "gtf" / "mine"
        directory.mkdir(parents=True)
        (directory / "mine.db").write_bytes(b"half a database")

        result = runner.invoke(
            app, ["register-gtf", "tiny", str(data_dir / "tiny.gtf"), "mine", "--force", "--json"]
        )

        assert result.exit_code == 0
        assert _json.loads(result.stdout)["sha256"] == _TINY_GTF_SHA256

    def test_the_chromosome_check_is_stood_down_from_the_command_line(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        # The committed Ensembl-spelled GTF (I, II, III) against a UCSC-spelled assembly
        # (chrI, chrII, chrIII): the assembly's chrom.sizes is found from its name, so
        # this way in checks the names too — and stands the check down when asked.
        source = data_dir / "ensembl_style.gtf"
        assembly_dir = tmp_path / "genome" / "tiny"
        assembly_dir.mkdir(parents=True)
        (assembly_dir / "tiny.chrom.sizes").write_text("chrI\t10000\nchrII\t10000\nchrIII\t10000\n")

        refused = runner.invoke(app, ["register-gtf", "tiny", str(source), "mine"])

        assert refused.exit_code == 1
        assert "chromosome" in _output(refused)

        result = runner.invoke(
            app, ["register-gtf", "tiny", str(source), "mine", "--no-check-chromosomes", "--json"]
        )

        assert result.exit_code == 0
        details = _json.loads(result.stdout)["details"]
        assert details["chromosomes_checked"] is False
        assert details["chromosomes_unchecked_because"] == "caller-override"

    def test_a_bare_exon_level_gtf_is_registrable_with_feature_inference(
        self, tmp_path: Path
    ) -> None:
        # Without the flags the database holds exons and nothing else — genes and
        # transcripts are what a caller registers an annotation for.
        source = tmp_path / "bare.gtf"
        source.write_text(_BARE_GTF)
        gtf_root = tmp_path / "genome" / "tiny" / "gtf"

        assert runner.invoke(app, ["register-gtf", "tiny", str(source), "exons"]).exit_code == 0
        assert _feature_types(gtf_root / "exons" / "exons.db") == ["exon"]

        result = runner.invoke(
            app,
            ["register-gtf", "tiny", str(source), "genes", "--infer-genes", "--infer-transcripts"],
        )

        assert result.exit_code == 0
        assert _feature_types(gtf_root / "genes" / "genes.db") == ["exon", "gene", "transcript"]


class TestWhatARegistrationSaysAboutTheChromosomes:
    """Both registration commands say which of four things happened to the name check.

    ``--no-check-chromosomes`` used to be answered with "register the assembly first",
    which the caller may well have done already: the record could not tell *nothing to
    check against* from *the caller stood the check down*, so the surface picked one and
    was wrong half the time. Each state now has its own sentence, and the one that had
    advice to give is the only one that gives any.
    """

    #: The advice that belongs to exactly one of the four states.
    _ADVICE = "register the assembly first"

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

    @staticmethod
    def _prepare_assembly(tmp_path: Path) -> None:
        """Put the assembly's ``chrom.sizes`` where the check looks for it."""
        assembly_dir = tmp_path / "genome" / "tiny"
        assembly_dir.mkdir(parents=True, exist_ok=True)
        (assembly_dir / "tiny.chrom.sizes").write_text("chrI\t10000\nchrII\t10000\nchrIII\t10000\n")

    def test_a_check_that_ran_is_reported_by_both_commands(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        # Reported rather than left to silence, which would read exactly the same as a
        # surface that had nothing good to say.
        self._prepare_assembly(tmp_path)

        by_name = runner.invoke(app, ["register-annotation", "tiny", "ensgene_v101"])
        by_path = runner.invoke(app, ["register-gtf", "tiny", str(data_dir / "tiny.gtf"), "mine"])

        for result in (by_name, by_path):
            assert result.exit_code == 0
            assert "chromosomes checked" in result.stdout
            assert self._ADVICE not in result.stdout

    def test_nothing_to_check_against_is_the_state_that_advises(self) -> None:
        # No chrom.sizes: the assembly is not registered here yet, and registering it is
        # exactly what would let the names be verified.
        result = runner.invoke(app, ["register-annotation", "tiny", "ensgene_v101"])

        assert result.exit_code == 0
        assert "chromosomes not checked" in result.stdout
        assert self._ADVICE in result.stdout

    def test_standing_the_check_down_is_advised_nothing_by_either_command(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        # The bug this fixes: the assembly is registered, the caller turned the check off
        # deliberately, and being told to register the assembly first is wrong.
        self._prepare_assembly(tmp_path)

        by_name = runner.invoke(
            app, ["register-annotation", "tiny", "ensgene_v101", "--no-check-chromosomes"]
        )
        by_path = runner.invoke(
            app,
            [
                "register-gtf",
                "tiny",
                str(data_dir / "ensembl_style.gtf"),
                "mine",
                "--no-check-chromosomes",
            ],
        )

        for result in (by_name, by_path):
            assert result.exit_code == 0
            assert "stood down" in result.stdout
            assert self._ADVICE not in result.stdout

    def test_a_record_written_before_the_reason_existed_reports_it_as_unknown(
        self, tmp_path: Path
    ) -> None:
        # An annotation registered by an older version, reported by re-running the
        # command over it: the record returned is the one already on disk, whose bare
        # `false` stands for either reason. Neither may be claimed, and neither raises.
        self._prepare_assembly(tmp_path)
        assert runner.invoke(app, ["register-annotation", "tiny", "ensgene_v101"]).exit_code == 0
        path = record_path(annotation_dir(tmp_path / "genome" / "tiny", "ensgene_v101"))
        written = _json.loads(path.read_text())
        written["details"] = {"chromosomes_checked": False}
        path.write_text(_json.dumps(written))

        result = runner.invoke(app, ["register-annotation", "tiny", "ensgene_v101"])

        assert result.exit_code == 0
        assert "does not say why" in result.stdout
        assert self._ADVICE not in result.stdout
        assert "stood down" not in result.stdout


class TestAnnotations:
    """``genome annotations`` — what the tables offer, set against what is registered here.

    The shipped table answers here rather than one stood up for the test: reporting it
    is this command's whole job. hg38 is the assembly, whose row offers one annotation
    and flags it as the default.
    """

    @pytest.fixture(autouse=True)
    def _offline(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("LIULAB_DATA", str(tmp_path))

    def test_an_assembly_with_nothing_registered_is_the_case_it_serves(
        self, tmp_path: Path
    ) -> None:
        result = runner.invoke(app, ["annotations", "hg38"])

        assert result.exit_code == 0
        assert "gencode_v50" in result.stdout
        assert "offered, not registered" in result.stdout
        assert "genome register-annotation hg38 gencode_v50" in result.stdout
        # Nothing was prepared to answer the question — the assembly is not even there.
        assert not (tmp_path / "genome" / "hg38").exists()

    def test_json(self, tmp_path: Path) -> None:
        result = runner.invoke(app, ["annotations", "hg38", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["assembly"] == "hg38"
        assert payload["directory"] == str(tmp_path / "genome" / "hg38")
        assert payload["default_annotation"] == "gencode_v50"
        assert [
            (row["name"], row["offered"], row["registered"]) for row in payload["annotations"]
        ] == [("gencode_v50", True, False)]

    def test_it_sets_what_is_registered_here_against_what_is_offered(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        register_gtf(tmp_path / "genome" / "hg38", data_dir / "tiny.gtf", "mine")

        result = runner.invoke(app, ["annotations", "hg38", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert [
            (row["name"], row["offered"], row["registered"]) for row in payload["annotations"]
        ] == [
            ("gencode_v50", True, False),
            ("mine", False, True),
        ]
        # The table's flag decides the default, whatever this machine happens to hold.
        assert payload["default_annotation"] == "gencode_v50"

    def test_an_assembly_no_row_lists_reports_an_empty_answer_rather_than_failing(self) -> None:
        result = runner.invoke(app, ["annotations", "tiny"])

        assert result.exit_code == 0
        assert "tiny" in result.stdout

    def test_a_broken_offered_annotation_reads_as_broken_and_names_its_repair(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        # It used to read `offered, not registered` — indistinguishable from one nobody
        # had ever fetched — and the closing line sent the reader to a command that
        # would itself raise and demand --force.
        assembly_dir = tmp_path / "genome" / "hg38"
        register_gtf(assembly_dir, data_dir / "tiny.gtf", "gencode_v50")
        record_path(annotation_dir(assembly_dir, "gencode_v50")).unlink()

        result = runner.invoke(app, ["annotations", "hg38"])

        assert result.exit_code == 0
        assert "offered, not registered" not in result.stdout
        assert "broken" in result.stdout
        assert "genome register-annotation hg38 gencode_v50 --force" in result.stdout
        default_line = next(
            line for line in result.stdout.splitlines() if line.startswith("default:")
        )
        assert "--force" in default_line

    def test_a_broken_unlisted_annotation_is_listed_at_all(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        assembly_dir = tmp_path / "genome" / "hg38"
        annotation = register_gtf(assembly_dir, data_dir / "tiny.gtf", "mine")
        annotation.db.write_bytes(b"truncated")

        result = runner.invoke(app, ["annotations", "hg38"])

        assert result.exit_code == 0
        assert "mine" in result.stdout
        assert "broken" in result.stdout
        assert f"genome register-gtf hg38 {data_dir / 'tiny.gtf'} mine --force" in result.stdout

    def test_json_carries_the_broken_state_and_the_repair(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        assembly_dir = tmp_path / "genome" / "hg38"
        register_gtf(assembly_dir, data_dir / "tiny.gtf", "healthy")
        annotation = register_gtf(assembly_dir, data_dir / "tiny.gtf", "mine")
        annotation.db.write_bytes(b"truncated")

        result = runner.invoke(app, ["annotations", "hg38", "--json"])

        assert result.exit_code == 0
        rows = {row["name"]: row for row in _json.loads(result.stdout)["annotations"]}
        assert rows["mine"]["broken"] is True
        assert rows["mine"]["registered"] is False
        assert rows["mine"]["repair"].endswith("mine --force")
        assert "mine.db" in rows["mine"]["problem"]
        # One broken annotation costs neither the exit code nor the ones beside it.
        assert rows["healthy"]["broken"] is False
        assert rows["healthy"]["registered"] is True
        assert rows["gencode_v50"]["broken"] is False


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
