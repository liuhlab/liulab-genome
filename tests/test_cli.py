"""Tests for the Typer CLI (``genome``)."""

from __future__ import annotations

import json as _json
import shutil
from dataclasses import dataclass, replace
from pathlib import Path

import gffutils
import pytest
from typer.testing import CliRunner

from genome import __version__ as genome_version
from genome import metadata
from genome.cli import app
from genome.external import REQUIRED_TOOLS
from genome.external import doctor as doctor_api
from genome.io import download as download_mod
from genome.io.completion import read_record, record_path, write_record
from genome.io.fasta import PREPARATION_TOOLS, GenomeFiles
from genome.io.gtf import annotation_dir, register_gtf
from genome.metadata import AnnotationMetadata, AssemblyMetadata
from genome.seq import DNA

from .conftest import CHIMERA_COMPONENTS, COMPONENT_ANNOTATION, FakeFetch

runner = CliRunner()
_BINARIES_PRESENT = all(shutil.which(t) is not None for t in REQUIRED_TOOLS)
_PREPARATION_PRESENT = all(shutil.which(t) is not None for t in PREPARATION_TOOLS)

#: sha256 of the committed ``tiny.fa``, which the fake fetch serves as any assembly.
_TINY_FA_SHA256 = "9316629bab14f9298a043f8b92e1e04a573b12d6a367ccc07c8f8040e5a13981"

#: sha256 of the committed ``tiny.gtf`` — the unpacked bytes ``tiny.gtf.gz`` yields.
_TINY_GTF_SHA256 = "255f43bd9abef76424d1c2d89a40cccc1a36215409bbc8f32dcead49ca3baf5e"

#: The URL the stood-in annotation row pins, served from ``tests/data``.
_ANNOTATION_URL = "https://mirror.example.invalid/annotations/tiny.gtf.gz"

#: The row the annotation table is stood in with: the committed ``tiny.gtf.gz``, pinned.
#: The CLI takes no metadata argument by design, so the table it reads is what moves.
_TINY_ANNOTATION = AnnotationMetadata(
    assembly="tiny",
    name="ensgene_v101",
    provider="UCSC",
    version="ensGene.v101",
    url=_ANNOTATION_URL,
    sha256=_TINY_GTF_SHA256,
    default=True,
)

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

    def test_version_json(self) -> None:
        result = runner.invoke(app, ["version", "--json"])
        assert result.exit_code == 0
        assert _json.loads(result.stdout) == {"version": genome_version}

    def test_json_and_text_report_the_same_version(self) -> None:
        text = runner.invoke(app, ["version"]).stdout.strip()
        payload = _json.loads(runner.invoke(app, ["version", "--json"]).stdout)
        assert payload["version"] == text


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
        assert "error" in _output(result).lower()
        assert "X" in _output(result)  # the base it could not complement

    def test_error_names_the_alphabet_it_was_held_to(self) -> None:
        result = runner.invoke(app, ["revcomp", "ATCX"])
        assert "{ACGT}" in _output(result)

    def test_alphabet_comes_from_the_type(self, monkeypatch: pytest.MonkeyPatch) -> None:
        # The defect this pins: the edge used to spell A/C/G/T itself, so changing what a
        # DNA may contain reached the CLI not at all. Widen the type's alphabet and the
        # command follows — both the check and the message naming it.
        monkeypatch.setattr(DNA, "ALPHABET", frozenset("ACGTN"))
        accepted = runner.invoke(app, ["revcomp", "ATCN"])
        assert accepted.exit_code == 0
        assert accepted.stdout.strip() == "NGAT"
        refused = runner.invoke(app, ["revcomp", "ATCX"])
        assert refused.exit_code == 2
        assert "{ACGNT}" in _output(refused)  # rendered sorted, a set having no order


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


class TestRegisterResolvesTheName:
    """What a name means, settled by four checks in order — and the two refusals.

    A record already here, then a source the caller named, then a name whose every part
    is prepared here or listed in the shipped table, then the download that was always
    the answer. Offline throughout, and the fetch step is recorded rather than merely
    stubbed: what a refusal is asserted on is that nothing was fetched at all, since
    turning one mistyped string into a whole-genome download per part is the failure the
    gate exists to prevent.
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

    def test_a_source_the_caller_named_settles_a_name_that_would_read_as_a_chimera(
        self, data_dir: Path, fake_fetch: FakeFetch
    ) -> None:
        # hg38 and mm10 are both listed, so the name alone reads as two assemblies and
        # would refuse on a machine holding neither. Saying where the bytes come from is
        # the caller answering the question, and it is believed.
        source = data_dir / "tiny.fa.gz"

        result = runner.invoke(app, ["register", "hg38_mm10", "--source", str(source), "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["source_url"] == str(source)
        assert "components" not in payload["details"]
        assert fake_fetch.calls == []

    def test_an_existing_record_says_what_to_rebuild_rather_than_the_name(
        self, data_dir: Path, fake_fetch: FakeFetch
    ) -> None:
        # The clause that stops a plain hg38_mm10 seeded years ago from silently becoming
        # a chimera: it was registered as an ordinary assembly, so that is what --force
        # registers again.
        assert (
            runner.invoke(
                app, ["register", "hg38_mm10", "--source", str(data_dir / "tiny.fa.gz")]
            ).exit_code
            == 0
        )

        result = runner.invoke(app, ["register", "hg38_mm10", "--force", "--json"])

        assert result.exit_code == 0
        assert "components" not in _json.loads(result.stdout)["details"]
        assert fake_fetch.last.url.endswith("hg38_mm10.fa.gz")

    def test_only_a_lost_record_falls_back_to_the_name(
        self, data_dir: Path, tmp_path: Path
    ) -> None:
        # …and with the record gone, the name is all that is left: the same directory now
        # reads as a chimera of hg38 and mm10, neither of which this machine has.
        assert (
            runner.invoke(
                app, ["register", "hg38_mm10", "--source", str(data_dir / "tiny.fa.gz")]
            ).exit_code
            == 0
        )
        record_path(tmp_path / "genome" / "hg38_mm10").unlink()

        result = runner.invoke(app, ["register", "hg38_mm10", "--force"])

        assert result.exit_code == 1
        assert "`genome register hg38`" in _output(result)
        assert "`genome register mm10`" in _output(result)

    def test_a_mis_ordered_name_is_refused_by_naming_the_canonical_spelling(
        self, fake_fetch: FakeFetch
    ) -> None:
        # What a --component flag would have bought, bought for less: the components are
        # in the name, so typing them in the wrong order is detectable.
        result = runner.invoke(app, ["register", "ecHT115_ce11"])

        assert result.exit_code == 1
        assert "`genome register ce11_ecHT115`" in _output(result)
        assert fake_fetch.calls == []

    def test_a_cold_machine_names_the_missing_component_rather_than_downloading(
        self, fake_fetch: FakeFetch
    ) -> None:
        result = runner.invoke(app, ["register", "ce11_ecHT115"])

        assert result.exit_code == 1
        assert "`genome register ce11`" in _output(result)
        assert "`genome register ecHT115`" in _output(result)
        assert fake_fetch.calls == []

    def test_force_is_not_a_bypass_of_the_gate(self, fake_fetch: FakeFetch) -> None:
        # It repairs a directory; it does not answer the question of what belongs in one.
        result = runner.invoke(app, ["register", "ce11_ecHT115", "--force"])

        assert result.exit_code == 1
        assert "`genome register ce11`" in _output(result)
        assert fake_fetch.calls == []

    def test_a_name_neither_prepared_nor_listed_is_downloaded_as_it_always_was(
        self, fake_fetch: FakeFetch
    ) -> None:
        # The whole separation between ce11_ecHT115 and a free-form local key: neither
        # `my` nor `ref` is an assembly here or in the table, so `my_ref` is one name
        # somebody chose and the download is the answer it always was (ADR-0003).
        result = runner.invoke(app, ["register", "my_ref", "--json"])

        assert result.exit_code == 0
        assert fake_fetch.last.url.endswith("my_ref.fa.gz")
        assert "components" not in _json.loads(result.stdout)["details"]


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
        # No row lists "tiny", so what it is held to is the digest its own registration
        # recorded — the fallback, and the payload says which answered.
        assert payload["expected"] == _TINY_FA_SHA256
        assert payload["expected_from"] == "record"
        assert payload["verified"] is True
        # An assembly that is not a chimera has no components to be asked about — null,
        # rather than a status that would read as a check somebody made.
        assert payload["components"] is None

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


class TestWhatAVerifiedDigestWasHeldTo:
    """Three answers, three sentences — and never one wording covering two of them.

    Being held to the digest the lab pinned, being held to the one this machine last
    produced, and being held to nothing at all are different results, and a caller who
    cannot tell them apart reads the weakest as the strongest.
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

    def test_the_row_pinned_it(self, monkeypatch: pytest.MonkeyPatch) -> None:
        row = AssemblyMetadata(
            assembly_name="tiny",
            species="Testus minimus",
            ucsc_name="tiny",
            ncbi_name="TINY.1",
            ncbi_assembly_id="GCF_0.0",
            ncbi_taxid=1,
            source_url="https://mirror.example.invalid/tiny.fa.gz",
            sha256=_TINY_FA_SHA256,
        )
        monkeypatch.setattr(metadata, "assembly_table", lambda: (row,))

        assert runner.invoke(app, ["register", "tiny"]).exit_code == 0
        result = runner.invoke(app, ["verify", "tiny"])

        assert result.exit_code == 0
        assert "matches the digest the metadata table pins for it" in result.stdout

    def test_the_record_pinned_it(self) -> None:
        assert runner.invoke(app, ["register", "tiny"]).exit_code == 0

        result = runner.invoke(app, ["verify", "tiny"])

        assert result.exit_code == 0
        assert "own registration recorded" in result.stdout
        assert "not an independent pin" in result.stdout

    def test_nothing_pinned_it(self, data_dir: Path) -> None:
        # A FASTA handed over by hand, checked against an assembly whose row pins nothing
        # and which is not registered here: there is no digest to be held to at all.
        result = runner.invoke(
            app, ["verify", "ce11_ecHT115", "--fasta", str(data_dir / "tiny.fa"), "--json"]
        )

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert (payload["expected"], payload["expected_from"], payload["verified"]) == (
            None,
            None,
            False,
        )
        # Nothing is asked about components either: the assembly's own registration is
        # not what is being verified.
        assert payload["components"] is None

        human = runner.invoke(app, ["verify", "ce11_ecHT115", "--fasta", str(data_dir / "tiny.fa")])

        assert "nothing to check it against" in human.stdout
        assert "components" not in human.stdout


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
        monkeypatch.setattr(metadata, "annotation_table", lambda: (_TINY_ANNOTATION,))

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
        row = replace(
            _TINY_ANNOTATION,
            url="https://mirror.example.invalid/annotations/ensembl_style.gtf",
            sha256=None,
        )
        monkeypatch.setattr(metadata, "annotation_table", lambda: (row,))

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
        row = replace(
            _TINY_ANNOTATION,
            name="bare",
            provider="somebody",
            version="1",
            url="https://mirror.example.invalid/annotations/bare.gtf",
            sha256=None,
        )
        monkeypatch.setattr(metadata, "annotation_table", lambda: (row,))

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
        monkeypatch.setattr(metadata, "annotation_table", lambda: (_TINY_ANNOTATION,))

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

    def test_a_chimera_is_refused_before_anything_is_downloaded(
        self, fake_fetch: FakeFetch
    ) -> None:
        # A chimera pins nothing, so this command has no job to do for one. The refusal
        # describes the row — the name, and every other column blank — rather than
        # printing a line that would look like one it computed something for.
        result = runner.invoke(app, ["table-row", "ce11_ecHT115"])

        assert result.exit_code == 1
        assert "ce11, ecHT115" in _output(result)
        assert "no sha256" in _output(result)
        assert "genome verify ce11_ecHT115" in _output(result)
        assert fake_fetch.calls == []


class TestTheSurfacesThatDidNotChange:
    """``annotations`` and ``doctor`` change nothing at all, asserted line for line.

    Both are pinned to their whole output rather than to a phrase inside it, because
    "nothing at all" is the claim: a line added to either of them fails here. The
    ``--json`` half is pinned the same way and for a stronger reason: a script parses it
    positionally as often as by key, so a reordered key is a break nobody would see.
    """

    @pytest.fixture(autouse=True)
    def _offline(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("LIULAB_DATA", str(tmp_path))

    def test_annotations_json_is_the_same_keys_in_the_same_order(self) -> None:
        result = runner.invoke(app, ["annotations", "hg38", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert list(payload) == ["assembly", "directory", "default_annotation", "annotations"]
        assert [list(row) for row in payload["annotations"]] == [
            [
                "name",
                "offered",
                "registered",
                "broken",
                "default",
                "provider",
                "version",
                "url",
                "sha256",
                "path",
                "problem",
                "repair",
            ]
        ]

    def test_register_gtf_json_is_the_record_then_the_two_facts_it_lacks(
        self, data_dir: Path
    ) -> None:
        result = runner.invoke(
            app, ["register-gtf", "tiny", str(data_dir / "tiny.gtf"), "mine", "--json"]
        )

        assert result.exit_code == 0
        assert list(_json.loads(result.stdout)) == [
            "kind",
            "name",
            "files",
            "source_url",
            "sha256",
            "tool_versions",
            "package_version",
            "completed_at",
            "details",
            "assembly",
            "directory",
        ]

    def test_annotations_prints_exactly_what_it_printed_before(self, tmp_path: Path) -> None:
        result = runner.invoke(app, ["annotations", "hg38"])

        assert result.exit_code == 0
        assert result.stdout == (
            f"annotations for hg38 in {tmp_path / 'genome' / 'hg38'}\n"
            f"  gencode_v50  offered, not registered  GENCODE v50\n"
            f"default: gencode_v50 — not registered here; register it with "
            f"`genome register-annotation hg38 gencode_v50`\n"
        )

    @pytest.mark.skipif(not _BINARIES_PRESENT, reason="samtools/bedtools not on PATH")
    def test_doctor_prints_one_line_per_tool_and_nothing_else(self) -> None:
        result = runner.invoke(app, ["doctor"])

        assert result.exit_code == 0
        assert result.stdout == "".join(f"{name}: {ver}\n" for name, ver in doctor_api().items())


#: The repair every chimera error names. Quoted here so the tests below can assert that a
#: message carries it and then run exactly it — a message naming a command nobody can
#: follow is worse than no message.
_CHIMERA_REPAIR = "genome register tinyCe_tinySc --force"


def _corrected_component(destination: Path) -> Path:
    """Write a valid FASTA carrying tinySc's chromosome names and different bases.

    What re-registering a component looks like: the same sequences by name, not by
    content, so a chimera built from the old bytes holds a copy of bytes that are no
    longer anywhere.
    """
    destination.write_text(">I\nACGTACGTAC\n>II\nACGTACGTAC\n>III\nACGTACGTAC\n")
    return destination


@pytest.mark.skipif(not _PREPARATION_PRESENT, reason="samtools/faToTwoBit/twoBitInfo not on PATH")
class TestChimeraFromTheCommandLine:
    """``genome register <name>`` is the only build spelling, end to end.

    The tiny components are registered through this same command and land under the
    shared root, which is where a component is looked for by name. Offline by
    construction — every source is a committed file, and a chimera fetches nothing — but
    the native tools are real, so what is asserted is what actually got written.
    """

    @pytest.fixture(autouse=True)
    def _root(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("LIULAB_DATA", str(tmp_path))

    @staticmethod
    def _register_components(*names: str, annotate: bool = False) -> None:
        """Register each named tiny component — and its own annotation when asked."""
        for name in names:
            component = CHIMERA_COMPONENTS[name]
            seeded = runner.invoke(app, ["register", name, "--source", str(component.fasta)])
            assert seeded.exit_code == 0, _output(seeded)
            if annotate:
                gtf = component.gtf
                assert gtf is not None
                built = runner.invoke(app, ["register-gtf", name, str(gtf), COMPONENT_ANNOTATION])
                assert built.exit_code == 0, _output(built)

    def test_naming_a_chimera_builds_it_and_the_report_says_what_it_is_made_of(
        self, tmp_path: Path
    ) -> None:
        # The line that used to read `source None` is the components, and the closing one
        # names the annotation this same command registered.
        self._register_components("tinyCe", "tinySc", annotate=True)

        result = runner.invoke(app, ["register", "tinyCe_tinySc"])

        assert result.exit_code == 0, _output(result)
        assert "  components  tinyCe, tinySc" in result.stdout
        assert "source" not in result.stdout
        assert f"  annotation  {COMPONENT_ANNOTATION}+{COMPONENT_ANNOTATION}" in result.stdout
        assert (tmp_path / "genome" / "tinyCe_tinySc" / "tinyCe_tinySc.fa").is_file()

    def test_a_build_with_nothing_to_merge_says_so_rather_than_saying_nothing(self) -> None:
        self._register_components("tinyCe", "tinySc")

        result = runner.invoke(app, ["register", "tinyCe_tinySc"])

        assert result.exit_code == 0, _output(result)
        assert "  annotation  none" in result.stdout

    def test_the_json_payload_is_the_record_untouched(self) -> None:
        # It already carried the components, which is why nothing was added to it.
        self._register_components("tinyCe", "tinySc")

        result = runner.invoke(app, ["register", "tinyCe_tinySc", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["source_url"] is None
        assert [entry["name"] for entry in payload["details"]["components"]] == [
            "tinyCe",
            "tinySc",
        ]

    def test_verify_says_the_components_are_unchanged(self) -> None:
        self._register_components("tinyCe", "tinySc")
        assert runner.invoke(app, ["register", "tinyCe_tinySc"]).exit_code == 0

        payload = _json.loads(runner.invoke(app, ["verify", "tinyCe_tinySc", "--json"]).stdout)
        human = runner.invoke(app, ["verify", "tinyCe_tinySc"])

        assert payload["components"] == "unchanged"
        assert "components  unchanged" in human.stdout

    def test_a_component_that_pinned_nothing_reads_as_unknown_rather_than_as_a_pass(
        self, tmp_path: Path
    ) -> None:
        # The line prints either way: a chimera whose components could not be compared is
        # unproven, and silence would be exactly what a pass looks like.
        self._register_components("tinyCe", "tinySc")
        assert runner.invoke(app, ["register", "tinyCe_tinySc"]).exit_code == 0
        directory = tmp_path / "genome" / "tinyCe_tinySc"
        record = read_record(directory)
        assert record is not None
        for entry in record.details["components"]:
            entry["sha256"] = None
        write_record(directory, record)

        payload = _json.loads(runner.invoke(app, ["verify", "tinyCe_tinySc", "--json"]).stdout)
        human = runner.invoke(app, ["verify", "tinyCe_tinySc"])

        assert payload["components"] == "unknown"
        assert human.exit_code == 0
        assert "components  unknown" in human.stdout

    def test_opening_by_name_catches_a_component_registered_again_underneath(
        self, tmp_path: Path
    ) -> None:
        # The hole this closes: opening by name returned from the chimera's own record,
        # which vouches for its files and can say nothing about the ones they were
        # copied from. Only building and verifying used to ask.
        self._register_components("tinyCe", "tinySc")
        assert runner.invoke(app, ["register", "tinyCe_tinySc"]).exit_code == 0
        corrected = str(_corrected_component(tmp_path / "corrected.fa"))
        assert (
            runner.invoke(app, ["register", "tinySc", "--force", "--source", corrected]).exit_code
            == 0
        )

        result = runner.invoke(app, ["register", "tinyCe_tinySc"])

        assert result.exit_code == 1
        assert _CHIMERA_REPAIR in _output(result)

    def test_the_repair_a_chimera_error_names_is_the_command_that_repairs_it(
        self, tmp_path: Path
    ) -> None:
        # Run verbatim, not paraphrased: this command used to route to the downloader and
        # fail with "Unknown UCSC assembly", so every chimera error quoted a repair nobody
        # could follow.
        self._register_components("tinyCe", "tinySc")
        assert runner.invoke(app, ["register", "tinyCe_tinySc"]).exit_code == 0
        corrected = str(_corrected_component(tmp_path / "corrected.fa"))
        assert (
            runner.invoke(app, ["register", "tinySc", "--force", "--source", corrected]).exit_code
            == 0
        )
        refused = runner.invoke(app, ["register", "tinyCe_tinySc"])
        assert _CHIMERA_REPAIR in _output(refused)

        repaired = runner.invoke(app, _CHIMERA_REPAIR.split()[1:])

        assert repaired.exit_code == 0, _output(repaired)
        # Rebuilt, not merely re-recorded: the corrected component's bases are in it.
        fasta = (tmp_path / "genome" / "tinyCe_tinySc" / "tinyCe_tinySc.fa").read_text()
        assert "ACGTACGTAC" in fasta
        verified = runner.invoke(app, ["verify", "tinyCe_tinySc", "--json"])
        assert _json.loads(verified.stdout)["components"] == "unchanged"

    def test_a_lost_record_is_rebuilt_from_the_name_by_the_command_it_names(
        self, tmp_path: Path
    ) -> None:
        # The residual a lost record leaves: the name is the only surviving information
        # about what this directory was, and it is enough.
        self._register_components("tinyCe", "tinySc")
        assert runner.invoke(app, ["register", "tinyCe_tinySc"]).exit_code == 0
        record_path(tmp_path / "genome" / "tinyCe_tinySc").unlink()
        refused = runner.invoke(app, ["register", "tinyCe_tinySc"])
        assert refused.exit_code == 1
        assert _CHIMERA_REPAIR in _output(refused)

        repaired = runner.invoke(app, _CHIMERA_REPAIR.split()[1:])

        assert repaired.exit_code == 0, _output(repaired)
        assert "  components  tinyCe, tinySc" in repaired.stdout
