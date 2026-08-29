"""Tests for the Typer CLI (``genome``)."""

from __future__ import annotations

import json as _json
import shutil
from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType

import gffutils
import pytest
from typer.testing import CliRunner

from genome import __version__ as genome_version
from genome import metadata
from genome.cli import app
from genome.external import REQUIRED_TOOLS
from genome.external import doctor as doctor_api
from genome.gene_list import curated_gene_list
from genome.io import download as download_mod
from genome.io.completion import read_record, record_path, write_record
from genome.io.fasta import PREPARATION_TOOLS, GenomeFiles
from genome.io.gtf import (
    AnnotationRegistry,
    GtfAnnotation,
    MergeSource,
    annotation_dir,
    register_merged_gtf,
)
from genome.metadata import AnnotationMetadata, AssemblyMetadata
from genome.seq import DNA
from genome.tf.gene import TFGeneTable, tf_gene_table
from genome.tf.motif import MIN_MOTIF_LENGTH, hit_count, provenance_of, read_hits
from genome.tf.motif import jaspar as jaspar_mod
from genome.tf.motif.jaspar import MOTIF_COUNTS

from .conftest import CHIMERA_COMPONENTS, COMPONENT_ANNOTATION, FakeFetch
from .test_jaspar import FIXTURE as _MOTIF_FIXTURE
from .test_jaspar import FIXTURE_COUNT as _MOTIF_COUNT
from .test_jaspar import FIXTURE_MOTIFS as _MOTIF_RECORDS
from .test_scan import FIXTURE as _PLANTED_FASTA

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

#: The two assemblies a census ships for, with the annotation each one defaults to. Named
#: rather than derived: which species has a census is exactly what the TF command is about,
#: so the pairing is written down where the test controls it. Two of them because the
#: censuses are two publishers' and their columns differ — mouse's has none beyond the
#: uniform four — and a command that only ever met Lambert's would not know that.
_TF_ASSEMBLY, _TF_ANNOTATION = "hg38", "gencode_v50"
_MOUSE_ASSEMBLY, _MOUSE_ANNOTATION = "mm39", "gencode_vM39"

#: Which of the committed motifs a scan leaves out, and how many are left to scan with.
#: Read off the fixture table rather than written down again, so a changed fixture moves the
#: expected summary with it: a motif under the minimum length cannot reach the default
#: **Threshold** at all, and is named among the skipped rather than called at something looser.
_MOTIFS_TOO_SHORT = [
    motif_id for motif_id, _name, length, _tax in _MOTIF_RECORDS if length < MIN_MOTIF_LENGTH
]
_MOTIFS_LONG_ENOUGH = _MOTIF_COUNT - len(_MOTIFS_TOO_SHORT)


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


def _register(assembly: str, assembly_dir: Path, gtf: Path, name: str) -> GtfAnnotation:
    """Register ``gtf`` under ``assembly_dir``, so a command has something to report on."""
    return AnnotationRegistry.locate(assembly, assembly_dir).register_path(gtf, name)


def _gtf_declaring(*gene_ids: str) -> str:
    """A GTF declaring one gene, transcript and exon per id, in the shape GENCODE has.

    The ids are the only thing that varies: nothing reading it looks at a coordinate, so
    every gene sits on the same interval rather than inviting anyone to.
    """
    return "".join(
        f'chrI\ttest\t{feature}\t1\t100\t.\t+\t.\tgene_id "{gene_id}"; transcript_id "{gene_id}_t";\n'
        for gene_id in gene_ids
        for feature in ("gene", "transcript", "exon")
    )


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


@pytest.fixture
def motif_release(fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch) -> FakeFetch:
    """Serve the committed transfac records as whichever **Release** is asked for.

    The arrangement ``tests/test_jaspar.py`` uses: the count check that stands where a
    **Completion marker** stands elsewhere is never switched off, only pointed at what the
    fake fetch actually serves.
    """
    monkeypatch.setattr(
        jaspar_mod, "MOTIF_COUNTS", MappingProxyType(dict.fromkeys(MOTIF_COUNTS, _MOTIF_COUNT))
    )
    fake_fetch.serve(_MOTIF_FIXTURE)
    return fake_fetch


@pytest.fixture
def planted_fasta(data_dir: Path) -> Path:
    """The committed FASTA with motifs planted at positions ``tests/data/README.md`` lists."""
    return data_dir / _PLANTED_FASTA


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
    any download, the ``HEAD`` name check is stubbed, and the shared ``liulab_data``
    fixture puts the assembly directory under this test's own root. The assembly is
    ``tiny``, which no shipped row lists, so nothing is pinned for the fixture to
    disagree with.
    """

    @pytest.fixture(autouse=True)
    def _offline(self, fake_fetch: FakeFetch, offline_prepare: None) -> None:
        fake_fetch.serve("tiny.fa.gz")

    def test_registers_and_reports_where_it_landed(self, liulab_data: Path) -> None:
        result = runner.invoke(app, ["register", "tiny"])

        assert result.exit_code == 0
        assert str(liulab_data / "genome" / "tiny") in result.stdout
        assert _TINY_FA_SHA256 in result.stdout
        assert (liulab_data / "genome" / "tiny" / "tiny.fa").is_file()

    def test_json(self, liulab_data: Path) -> None:
        result = runner.invoke(app, ["register", "tiny", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["assembly"] == "tiny"
        assert payload["directory"] == str(liulab_data / "genome" / "tiny")
        assert payload["sha256"] == _TINY_FA_SHA256
        assert sorted(payload["files"]) == [
            "tiny.2bit",
            "tiny.chrom.sizes",
            "tiny.fa",
            "tiny.fa.fai",
        ]

    def test_a_broken_directory_exits_non_zero_naming_the_repair(self, liulab_data: Path) -> None:
        directory = liulab_data / "genome" / "tiny"
        directory.mkdir(parents=True)
        (directory / "tiny.fa").write_text("half a genome\n")

        result = runner.invoke(app, ["register", "tiny"])

        assert result.exit_code == 1
        assert "genome register tiny --force" in _output(result)

    def test_force_repairs_what_the_error_named(self, liulab_data: Path) -> None:
        directory = liulab_data / "genome" / "tiny"
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
    def _offline(self, fake_fetch: FakeFetch, offline_prepare: None) -> None:
        fake_fetch.serve("tiny.fa.gz")

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
        self, data_dir: Path, liulab_data: Path
    ) -> None:
        # …and with the record gone, the name is all that is left: the same directory now
        # reads as a chimera of hg38 and mm10, neither of which this machine has.
        assert (
            runner.invoke(
                app, ["register", "hg38_mm10", "--source", str(data_dir / "tiny.fa.gz")]
            ).exit_code
            == 0
        )
        record_path(liulab_data / "genome" / "hg38_mm10").unlink()

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
    def _offline(self, fake_fetch: FakeFetch, offline_prepare: None) -> None:
        fake_fetch.serve("tiny.fa.gz")

    def test_reports_the_digest_of_a_registered_assembly(self) -> None:
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
    def _offline(self, fake_fetch: FakeFetch, offline_prepare: None) -> None:
        fake_fetch.serve("tiny.fa.gz")

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
    def _offline(self, fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch) -> None:
        fake_fetch.serve("tiny.gtf.gz")
        monkeypatch.setattr(metadata, "annotation_table", lambda: (_TINY_ANNOTATION,))

    def test_registers_and_reports_where_it_landed(self, liulab_data: Path) -> None:
        result = runner.invoke(app, ["register-annotation", "tiny", "ensgene_v101"])

        directory = liulab_data / "genome" / "tiny" / "gtf" / "ensgene_v101"
        assert result.exit_code == 0
        assert str(directory) in result.stdout
        assert _TINY_GTF_SHA256 in result.stdout
        assert (directory / "ensgene_v101.gtf").is_file()
        assert (directory / "ensgene_v101.db").is_file()

    def test_json(self, liulab_data: Path) -> None:
        result = runner.invoke(app, ["register-annotation", "tiny", "ensgene_v101", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["assembly"] == "tiny"
        assert payload["name"] == "ensgene_v101"
        assert payload["directory"] == str(liulab_data / "genome" / "tiny" / "gtf" / "ensgene_v101")
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

    def test_a_broken_directory_exits_non_zero_naming_the_repair(self, liulab_data: Path) -> None:
        directory = liulab_data / "genome" / "tiny" / "gtf" / "ensgene_v101"
        directory.mkdir(parents=True)
        (directory / "ensgene_v101.db").write_bytes(b"half a database")

        result = runner.invoke(app, ["register-annotation", "tiny", "ensgene_v101"])

        assert result.exit_code == 1
        assert "genome register-annotation tiny ensgene_v101 --force" in _output(result)

    def test_force_repairs_what_the_error_named(self, liulab_data: Path) -> None:
        directory = liulab_data / "genome" / "tiny" / "gtf" / "ensgene_v101"
        directory.mkdir(parents=True)
        (directory / "ensgene_v101.db").write_bytes(b"half a database")

        result = runner.invoke(
            app, ["register-annotation", "tiny", "ensgene_v101", "--force", "--json"]
        )

        assert result.exit_code == 0
        assert _json.loads(result.stdout)["sha256"] == _TINY_GTF_SHA256

    def test_the_chromosome_check_is_stood_down_from_the_command_line(
        self, fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch, liulab_data: Path
    ) -> None:
        # The committed Ensembl-spelled GTF (I, II, III) against a UCSC-spelled
        # assembly (chrI, chrII, chrIII): refused by default, and registered anyway
        # once the caller says they have looked at the mismatch and accept it.
        fake_fetch.serve("ensembl_style.gtf")
        assembly_dir = liulab_data / "genome" / "tiny"
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
        self,
        fake_fetch: FakeFetch,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        liulab_data: Path,
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
        database = liulab_data / "genome" / "tiny" / "gtf" / "bare" / "bare.db"
        assert _feature_types(database) == ["exon", "gene", "transcript"]


class TestRegisterGtf:
    """``genome register-gtf`` — register a GTF the annotation table does not list.

    The by-path way in, from a shell: no table row, no download, no checksum to compare
    against — the caller says where the file is. Offline by construction, since the GTF
    is a local one, and the shared ``liulab_data`` fixture puts the assembly directory
    under this test's own root.
    """

    def test_registers_a_gtf_no_row_lists_and_reports_where_it_landed(
        self, data_dir: Path, liulab_data: Path
    ) -> None:
        source = data_dir / "tiny.gtf"

        result = runner.invoke(app, ["register-gtf", "tiny", str(source), "mine"])

        directory = liulab_data / "genome" / "tiny" / "gtf" / "mine"
        assert result.exit_code == 0
        assert str(directory) in result.stdout
        assert str(source) in result.stdout
        assert (directory / "mine.gtf").is_file()
        assert (directory / "mine.db").is_file()

    def test_json(self, data_dir: Path, liulab_data: Path) -> None:
        source = data_dir / "tiny.gtf"

        result = runner.invoke(app, ["register-gtf", "tiny", str(source), "mine", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["assembly"] == "tiny"
        assert payload["name"] == "mine"
        assert payload["directory"] == str(liulab_data / "genome" / "tiny" / "gtf" / "mine")
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
        self, data_dir: Path, liulab_data: Path
    ) -> None:
        source = data_dir / "tiny.gtf"
        directory = liulab_data / "genome" / "tiny" / "gtf" / "mine"
        directory.mkdir(parents=True)
        (directory / "mine.db").write_bytes(b"half a database")

        result = runner.invoke(app, ["register-gtf", "tiny", str(source), "mine"])

        assert result.exit_code == 1
        assert f"genome register-gtf tiny {source} mine --force" in _output(result)

    def test_force_repairs_what_the_error_named(self, data_dir: Path, liulab_data: Path) -> None:
        directory = liulab_data / "genome" / "tiny" / "gtf" / "mine"
        directory.mkdir(parents=True)
        (directory / "mine.db").write_bytes(b"half a database")

        result = runner.invoke(
            app, ["register-gtf", "tiny", str(data_dir / "tiny.gtf"), "mine", "--force", "--json"]
        )

        assert result.exit_code == 0
        assert _json.loads(result.stdout)["sha256"] == _TINY_GTF_SHA256

    def test_the_chromosome_check_is_stood_down_from_the_command_line(
        self, data_dir: Path, liulab_data: Path
    ) -> None:
        # The committed Ensembl-spelled GTF (I, II, III) against a UCSC-spelled assembly
        # (chrI, chrII, chrIII): the assembly's chrom.sizes is found from its name, so
        # this way in checks the names too — and stands the check down when asked.
        source = data_dir / "ensembl_style.gtf"
        assembly_dir = liulab_data / "genome" / "tiny"
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
        self, tmp_path: Path, liulab_data: Path
    ) -> None:
        # Without the flags the database holds exons and nothing else — genes and
        # transcripts are what a caller registers an annotation for.
        source = tmp_path / "bare.gtf"
        source.write_text(_BARE_GTF)
        gtf_root = liulab_data / "genome" / "tiny" / "gtf"

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
    def _offline(self, fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch) -> None:
        fake_fetch.serve("tiny.gtf.gz")
        monkeypatch.setattr(metadata, "annotation_table", lambda: (_TINY_ANNOTATION,))

    @staticmethod
    def _prepare_assembly(liulab_data: Path) -> None:
        """Put the assembly's ``chrom.sizes`` where the check looks for it."""
        assembly_dir = liulab_data / "genome" / "tiny"
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
        self, tmp_path: Path, liulab_data: Path
    ) -> None:
        # An annotation registered by an older version, reported by re-running the
        # command over it: the record returned is the one already on disk, whose bare
        # `false` stands for either reason. Neither may be claimed, and neither raises.
        self._prepare_assembly(tmp_path)
        assert runner.invoke(app, ["register-annotation", "tiny", "ensgene_v101"]).exit_code == 0
        path = record_path(annotation_dir(liulab_data / "genome" / "tiny", "ensgene_v101"))
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

    def test_an_assembly_with_nothing_registered_is_the_case_it_serves(
        self, liulab_data: Path
    ) -> None:
        result = runner.invoke(app, ["annotations", "hg38"])

        assert result.exit_code == 0
        assert "gencode_v50" in result.stdout
        assert "offered, not registered" in result.stdout
        assert "genome register-annotation hg38 gencode_v50" in result.stdout
        # Nothing was prepared to answer the question — the assembly is not even there.
        assert not (liulab_data / "genome" / "hg38").exists()

    def test_json(self, liulab_data: Path) -> None:
        result = runner.invoke(app, ["annotations", "hg38", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["assembly"] == "hg38"
        assert payload["directory"] == str(liulab_data / "genome" / "hg38")
        assert payload["default_annotation"] == "gencode_v50"
        assert [
            (row["name"], row["offered"], row["registered"]) for row in payload["annotations"]
        ] == [("gencode_v50", True, False)]

    def test_it_sets_what_is_registered_here_against_what_is_offered(
        self, data_dir: Path, liulab_data: Path
    ) -> None:
        _register("hg38", liulab_data / "genome" / "hg38", data_dir / "tiny.gtf", "mine")

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
        self, data_dir: Path, liulab_data: Path
    ) -> None:
        # It used to read `offered, not registered` — indistinguishable from one nobody
        # had ever fetched — and the closing line sent the reader to a command that
        # would itself raise and demand --force.
        assembly_dir = liulab_data / "genome" / "hg38"
        _register("hg38", assembly_dir, data_dir / "tiny.gtf", "gencode_v50")
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
        self, data_dir: Path, liulab_data: Path
    ) -> None:
        assembly_dir = liulab_data / "genome" / "hg38"
        annotation = _register("hg38", assembly_dir, data_dir / "tiny.gtf", "mine")
        annotation.db.write_bytes(b"truncated")

        result = runner.invoke(app, ["annotations", "hg38"])

        assert result.exit_code == 0
        assert "mine" in result.stdout
        assert "broken" in result.stdout
        assert f"genome register-gtf hg38 {data_dir / 'tiny.gtf'} mine --force" in result.stdout

    def test_json_carries_the_broken_state_and_the_repair(
        self, data_dir: Path, liulab_data: Path
    ) -> None:
        assembly_dir = liulab_data / "genome" / "hg38"
        _register("hg38", assembly_dir, data_dir / "tiny.gtf", "healthy")
        annotation = _register("hg38", assembly_dir, data_dir / "tiny.gtf", "mine")
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
    any download, and the shared ``liulab_data`` fixture puts the assembly directory
    under this test's own root. hg38 and sacCer3 are used because the shipped table pins a
    source URL for both, which also skips the network name check.
    """

    @pytest.fixture(autouse=True)
    def _offline(self, fake_fetch: FakeFetch) -> None:
        fake_fetch.serve("tiny.fa.gz")

    def test_prints_the_row_to_paste(self) -> None:
        result = runner.invoke(app, ["table-row", "hg38"])

        assert result.exit_code == 0
        row = result.stdout.strip().split("\t")
        assert row[:8] == [
            "hg38",
            "Homo sapiens",
            "hg38",
            "GRCh38",
            "GCF_000001405.40",
            "9606",
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
            _TINY_FA_SHA256,
        ]
        # Every curated column rides along, not only the two this command computes: the
        # line is pasted over the row it replaces, so a cell dropped here is a cell
        # deleted from the table. hg38's intron bound is one nothing could recompute.
        assert row[metadata.METADATA_FIELDS.index("intron_length_cap")] == "1000000"
        assert row[metadata.METADATA_FIELDS.index("intron_length_cap_rationale")]

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
        assert row[metadata.METADATA_FIELDS.index("sha256")] == _TINY_FA_SHA256

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


class TestGeneCategoryCommands:
    """``genome gene-list`` and ``genome gene-categories`` — a category's genes, and which exist.

    The shipped curated gene lists answer, since which categories exist is data and no
    fixture may pretend otherwise: the tests read them off the shipped file rather than
    naming any. ``sacCer3``/``ensgene_v101`` is the pair used throughout, because its
    curated list is what makes an annotation with no biotype attribute answerable at all.
    """

    #: The category names a shipped curated list declares, in the order it spells them.
    def _declared(self, annotation: str) -> tuple[str, ...]:
        listed = curated_gene_list(annotation)
        assert listed is not None, f"no curated gene list ships for {annotation}"
        return tuple(listed.categories)

    def _ids(self, annotation: str, category: str) -> tuple[str, ...]:
        listed = curated_gene_list(annotation)
        assert listed is not None
        return listed.categories[category].gene_ids

    def _registered(self, liulab_data: Path, data_dir: Path) -> None:
        """Register the fixture GTF as sacCer3's ``ensgene_v101``, where the layout puts it."""
        _register(
            "sacCer3", liulab_data / "genome" / "sacCer3", data_dir / "tiny.gtf", "ensgene_v101"
        )

    def _merged(self, liulab_data: Path, tmp_path: Path) -> None:
        """Register a merged annotation of worm and its food, where the layout puts it."""
        source = tmp_path / "one.gtf"
        source.write_text('chrI\ttest\texon\t1\t50\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n')
        assembly_dir = liulab_data / "genome" / "ce11_ecHT115"
        assembly_dir.mkdir(parents=True, exist_ok=True)
        chrom_sizes = assembly_dir / "ce11_ecHT115.chrom.sizes"
        chrom_sizes.write_text("chrI__ce11\t10000\nchrI__ecHT115\t10000\n")
        register_merged_gtf(
            assembly_dir,
            "wormbase_ws298+refseq_rs_2025_06_26",
            [
                MergeSource("ce11", "wormbase_ws298", source),
                MergeSource("ecHT115", "refseq_rs_2025_06_26", source),
            ],
            separator="__",
            chrom_sizes=chrom_sizes,
        )

    def test_only_the_gene_ids_reach_stdout_so_the_output_pipes(
        self, liulab_data: Path, data_dir: Path
    ) -> None:
        self._registered(liulab_data, data_dir)
        category = self._declared("ensgene_v101")[0]

        result = runner.invoke(app, ["gene-list", "sacCer3", category])

        assert result.exit_code == 0
        assert result.stdout == "".join(
            f"{gene_id}\n" for gene_id in self._ids("ensgene_v101", category)
        )
        # The attribution is worth printing and must not cost the pipe, so it goes beside it.
        assert category in result.stderr
        assert "ensgene_v101" in result.stderr

    def test_the_annotation_may_be_named_instead_of_defaulted(
        self, liulab_data: Path, data_dir: Path
    ) -> None:
        self._registered(liulab_data, data_dir)
        category = self._declared("ensgene_v101")[0]

        result = runner.invoke(
            app, ["gene-list", "sacCer3", category, "--annotation", "ensgene_v101"]
        )

        assert result.exit_code == 0
        assert result.stdout.splitlines() == list(self._ids("ensgene_v101", category))

    def test_gene_list_json_carries_the_ids_and_keeps_the_sources_apart(
        self, liulab_data: Path, data_dir: Path
    ) -> None:
        self._registered(liulab_data, data_dir)
        category = self._declared("ensgene_v101")[0]

        result = runner.invoke(app, ["gene-list", "sacCer3", category, "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert list(payload) == ["assembly", "annotation", "category", "gene_ids", "sources"]
        assert payload["assembly"] == "sacCer3"
        assert payload["gene_ids"] == list(self._ids("ensgene_v101", category))
        assert [list(source) for source in payload["sources"]] == [
            ["component", "annotation", "description", "source", "gene_ids"]
        ]
        assert payload["sources"][0]["component"] is None

    def test_gene_categories_prints_one_row_per_category_with_its_count(
        self, liulab_data: Path, data_dir: Path
    ) -> None:
        self._registered(liulab_data, data_dir)
        declared = self._declared("ensgene_v101")

        result = runner.invoke(app, ["gene-categories", "sacCer3"])

        assert result.exit_code == 0
        lines = result.stdout.splitlines()
        assert lines[0] == "categories for sacCer3 / ensgene_v101"
        assert [line.split()[0] for line in lines[1:]] == list(declared)
        assert [line.split()[1] for line in lines[1:]] == [
            str(len(self._ids("ensgene_v101", category))) for category in declared
        ]

    def test_gene_categories_json_is_every_category_as_gene_list_answers_one(
        self, liulab_data: Path, data_dir: Path
    ) -> None:
        self._registered(liulab_data, data_dir)

        result = runner.invoke(app, ["gene-categories", "sacCer3", "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert [entry["category"] for entry in payload] == list(self._declared("ensgene_v101"))
        assert all(entry["gene_ids"] for entry in payload)

    def test_a_merged_annotation_shows_the_per_component_split(
        self, liulab_data: Path, tmp_path: Path
    ) -> None:
        # The case #111 was opened over: worm rRNA and its food's arrive as one category
        # and must stay distinguishable inside it.
        self._merged(liulab_data, tmp_path)

        result = runner.invoke(app, ["gene-categories", "ce11_ecHT115"])

        assert result.exit_code == 0
        lines = result.stdout.splitlines()
        assert lines[0] == "categories for ce11_ecHT115 / wormbase_ws298+refseq_rs_2025_06_26"
        assert any("(ce11: " in line and "ecHT115: " in line for line in lines[1:])

    def test_a_merged_annotations_gene_list_attributes_every_source(
        self, liulab_data: Path, tmp_path: Path
    ) -> None:
        self._merged(liulab_data, tmp_path)
        shared = next(
            category
            for category in self._declared("wormbase_ws298")
            if category in self._declared("refseq_rs_2025_06_26")
        )

        result = runner.invoke(app, ["gene-list", "ce11_ecHT115", shared, "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert [source["component"] for source in payload["sources"]] == ["ce11", "ecHT115"]
        assert payload["gene_ids"] == [
            *self._ids("wormbase_ws298", shared),
            *self._ids("refseq_rs_2025_06_26", shared),
        ]

    def test_an_annotation_no_curated_list_ships_for_exits_one_saying_so(
        self, liulab_data: Path, data_dir: Path
    ) -> None:
        # Not an empty list of genes and not exit 0: the caller has to be able to tell
        # *nothing is known here* from *there are none of these genes*.
        _register("tiny", liulab_data / "genome" / "tiny", data_dir / "tiny.gtf", "mine")

        result = runner.invoke(app, ["gene-list", "tiny", "rRNA"])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert "no curated gene list ships" in _output(result)
        assert "ensgene_v101" in _output(result)  # …and which annotations do have one

    def test_a_category_the_annotation_does_not_declare_exits_one_listing_the_ones_it_does(
        self, liulab_data: Path, data_dir: Path
    ) -> None:
        self._registered(liulab_data, data_dir)

        result = runner.invoke(app, ["gene-list", "sacCer3", "no_such_category"])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert "no_such_category" in _output(result)
        for category in self._declared("ensgene_v101"):
            assert category in _output(result)

    def test_an_unregistered_annotation_exits_one_naming_the_command_that_registers_it(
        self, liulab_data: Path
    ) -> None:
        result = runner.invoke(app, ["gene-categories", "sacCer3"])

        assert result.exit_code == 1
        assert "genome register-annotation sacCer3 ensgene_v101" in _output(result)


class TestTFGeneListCommand:
    """``genome tf-gene-list`` — an assembly's TF genes, shaped like ``genome gene-list``.

    The shipped censuses answer, as the shipped curated lists answer for gene categories:
    which genes are transcription factors is the census's judgement and no fixture stands
    in for it, so every expectation below is read off the shipped file. What is asserted
    here is the command and not the crossing — ``tests/test_gtf.py`` owns that — so: the
    stdout/stderr split that makes the output pipe, the record ``--json`` emits, and a
    non-zero exit naming the next action for each of the three ways it can fail.
    """

    def _census(self, assembly: str) -> TFGeneTable:
        """The census shipped for ``assembly``'s species, whichever publisher wrote it."""
        species = metadata.assembly_metadata(assembly).species
        assert species is not None, f"{assembly} has no species in the assembly table"
        census = tf_gene_table(species)
        assert census is not None, f"no census ships for {species}"
        return census

    def _registered(
        self,
        liulab_data: Path,
        *gene_ids: str,
        assembly: str = _TF_ASSEMBLY,
        name: str = _TF_ANNOTATION,
    ) -> None:
        """Register a GTF declaring ``gene_ids`` as ``assembly``'s ``name``, where it lives."""
        source = liulab_data / f"{assembly}.{name}.gtf"
        source.write_text(_gtf_declaring(*gene_ids))
        _register(assembly, liulab_data / "genome" / assembly, source, name)

    def _positive(self, assembly: str, count: int) -> list[str]:
        """``count`` gene ids, one per assessed-positive stem, versioned as GENCODE spells them."""
        return [f"{stem}.1" for stem in self._census(assembly).assessed_positive[:count]]

    def test_only_the_gene_ids_reach_stdout_so_the_output_pipes(self, liulab_data: Path) -> None:
        gene_ids = self._positive(_TF_ASSEMBLY, 2)
        self._registered(liulab_data, *gene_ids)

        result = runner.invoke(app, ["tf-gene-list", _TF_ASSEMBLY])

        assert result.exit_code == 0
        assert result.stdout == "".join(f"{gene_id}\n" for gene_id in gene_ids)
        # Whose judgement it is must be printed and must not cost the pipe, so it goes
        # beside the ids: the heading, the census's own attribution, and what the crossing
        # cost — the stems the census holds that this annotation carries no gene for.
        assert f"{_TF_ASSEMBLY} / {_TF_ANNOTATION}" in result.stderr
        assert "Homo sapiens" in result.stderr
        assert self._census(_TF_ASSEMBLY).provenance.attribution() in result.stderr
        unresolved = len(self._census(_TF_ASSEMBLY).assessed_positive) - len(gene_ids)
        assert f"2 genes, 2 gene ids, {unresolved} stems" in result.stderr

    def test_the_annotation_may_be_named_instead_of_defaulted(self, liulab_data: Path) -> None:
        gene_ids = self._positive(_TF_ASSEMBLY, 1)
        self._registered(liulab_data, *gene_ids, name="mine")

        result = runner.invoke(app, ["tf-gene-list", _TF_ASSEMBLY, "--annotation", "mine"])

        assert result.exit_code == 0
        assert result.stdout.splitlines() == gene_ids
        assert f"{_TF_ASSEMBLY} / mine" in result.stderr

    def test_json_carries_the_genes_the_provenance_and_the_unresolved_stems(
        self, liulab_data: Path
    ) -> None:
        census = self._census(_TF_ASSEMBLY)
        stem = census.assessed_positive[0]
        self._registered(liulab_data, f"{stem}.1")

        result = runner.invoke(app, ["tf-gene-list", _TF_ASSEMBLY, "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert list(payload) == [
            "assembly",
            "annotation",
            "species",
            "provenance",
            "genes",
            "gene_ids",
            "unresolved",
        ]
        assert (payload["assembly"], payload["annotation"]) == (_TF_ASSEMBLY, _TF_ANNOTATION)
        assert payload["gene_ids"] == [f"{stem}.1"]
        assert payload["provenance"]["pubmed_id"] == census.provenance.pubmed_id
        assert payload["unresolved"]  # what this annotation carries no gene for, not dropped
        cells = dict(
            zip(census.columns, census.rows[census.gene_id_stems.index(stem)], strict=True)
        )
        gene = payload["genes"][0]
        assert gene["dbd_family"] == cells["dbd_family"]
        assert gene["judgements"]["tf_assessment"] == cells["tf_assessment"]

    def test_a_mouse_assembly_is_answered_by_the_census_published_for_mouse(
        self, liulab_data: Path
    ) -> None:
        # The verdict travels with the census that reached it, so which one spoke is a
        # fact about the assembly's species and never about which one came first.
        gene_ids = self._positive(_MOUSE_ASSEMBLY, 1)
        self._registered(liulab_data, *gene_ids, assembly=_MOUSE_ASSEMBLY, name=_MOUSE_ANNOTATION)

        result = runner.invoke(app, ["tf-gene-list", _MOUSE_ASSEMBLY])

        assert result.exit_code == 0
        assert result.stdout.splitlines() == gene_ids
        assert "Mus musculus" in result.stderr
        assert self._census(_MOUSE_ASSEMBLY).provenance.attribution() in result.stderr
        assert self._census(_TF_ASSEMBLY).provenance.publisher not in result.stderr

    def test_a_census_recording_nothing_beyond_the_uniform_columns_still_emits_json(
        self, liulab_data: Path
    ) -> None:
        # AnimalTFDB ships the four uniform columns and no more, so every mouse gene's
        # judgements are empty. A surface reaching for a **TF assessment** that is not
        # there would raise here rather than print, which is why nothing does.
        gene_ids = self._positive(_MOUSE_ASSEMBLY, 1)
        self._registered(liulab_data, *gene_ids, assembly=_MOUSE_ASSEMBLY, name=_MOUSE_ANNOTATION)

        result = runner.invoke(app, ["tf-gene-list", _MOUSE_ASSEMBLY, "--json"])

        assert result.exit_code == 0
        gene = _json.loads(result.stdout)["genes"][0]
        assert gene["judgements"] == {}
        assert gene["dbd_family"]  # the uniform four are there all the same

    def test_an_unregistered_annotation_exits_one_naming_the_command_that_registers_it(
        self, liulab_data: Path
    ) -> None:
        result = runner.invoke(app, ["tf-gene-list", _TF_ASSEMBLY])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert f"genome register-annotation {_TF_ASSEMBLY} {_TF_ANNOTATION}" in _output(result)

    def test_a_species_no_census_ships_for_exits_one_naming_the_species_that_have_one(
        self, liulab_data: Path
    ) -> None:
        # Human gene ids registered for a worm assembly: the species is the assembly's own
        # and never what the GTF happens to hold, so this is refused rather than answered.
        self._registered(
            liulab_data, *self._positive(_TF_ASSEMBLY, 1), assembly="ce11", name="wormbase_ws298"
        )

        result = runner.invoke(app, ["tf-gene-list", "ce11"])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert "no TF census ships" in _output(result)
        assert "Caenorhabditis elegans" in _output(result)
        assert "Homo sapiens" in _output(result)  # …and what may be asked about instead

    def test_an_assembly_nothing_names_a_species_for_exits_one_saying_so(
        self, liulab_data: Path
    ) -> None:
        # Not the same fact as no census ships: the question was which species this is,
        # and no row answered it. An unlisted local key is the ordinary way in.
        self._registered(
            liulab_data, *self._positive(_TF_ASSEMBLY, 1), assembly="tiny", name="mine"
        )

        result = runner.invoke(app, ["tf-gene-list", "tiny", "--annotation", "mine"])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert "nothing says what species 'tiny' is" in _output(result)
        assert "Homo sapiens" in _output(result)


class TestTheSurfacesThatDidNotChange:
    """``annotations`` and ``doctor`` change nothing at all, asserted line for line.

    Both are pinned to their whole output rather than to a phrase inside it, because
    "nothing at all" is the claim: a line added to either of them fails here. The
    ``--json`` half is pinned the same way and for a stronger reason: a script parses it
    positionally as often as by key, so a reordered key is a break nobody would see.
    """

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
        self, liulab_data: Path
    ) -> None:
        # The line that used to read `source None` is the components, and the closing one
        # names the annotation this same command registered.
        self._register_components("tinyCe", "tinySc", annotate=True)

        result = runner.invoke(app, ["register", "tinyCe_tinySc"])

        assert result.exit_code == 0, _output(result)
        assert "  components  tinyCe, tinySc" in result.stdout
        assert "source" not in result.stdout
        assert f"  annotation  {COMPONENT_ANNOTATION}+{COMPONENT_ANNOTATION}" in result.stdout
        assert (liulab_data / "genome" / "tinyCe_tinySc" / "tinyCe_tinySc.fa").is_file()

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
        self, liulab_data: Path
    ) -> None:
        # The line prints either way: a chimera whose components could not be compared is
        # unproven, and silence would be exactly what a pass looks like.
        self._register_components("tinyCe", "tinySc")
        assert runner.invoke(app, ["register", "tinyCe_tinySc"]).exit_code == 0
        directory = liulab_data / "genome" / "tinyCe_tinySc"
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
        self, tmp_path: Path, liulab_data: Path
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
        fasta = (liulab_data / "genome" / "tinyCe_tinySc" / "tinyCe_tinySc.fa").read_text()
        assert "ACGTACGTAC" in fasta
        verified = runner.invoke(app, ["verify", "tinyCe_tinySc", "--json"])
        assert _json.loads(verified.stdout)["components"] == "unchanged"

    def test_a_lost_record_is_rebuilt_from_the_name_by_the_command_it_names(
        self, liulab_data: Path
    ) -> None:
        # The residual a lost record leaves: the name is the only surviving information
        # about what this directory was, and it is enough.
        self._register_components("tinyCe", "tinySc")
        assert runner.invoke(app, ["register", "tinyCe_tinySc"]).exit_code == 0
        record_path(liulab_data / "genome" / "tinyCe_tinySc").unlink()
        refused = runner.invoke(app, ["register", "tinyCe_tinySc"])
        assert refused.exit_code == 1
        assert _CHIMERA_REPAIR in _output(refused)

        repaired = runner.invoke(app, _CHIMERA_REPAIR.split()[1:])

        assert repaired.exit_code == 0, _output(repaired)
        assert "  components  tinyCe, tinySc" in repaired.stdout


# ---------------------------------------------------------------------------
# `genome motif-scan` — a FASTA in, a Parquet file out, a summary on stdout
# ---------------------------------------------------------------------------


class TestMotifScan:
    """The batch scan: what the summary says, where each half of the answer goes.

    Offline like every other test here — the ``fake_fetch`` fixture serves the committed
    transfac records as whatever release is asked for, with the count check pointed at
    them — and unmarked, since a process pool is not a binary this package ships.
    """

    def _summary(self, result: object) -> dict[str, object]:
        """Parse the JSON summary, which is the whole of what the command puts on stdout."""
        payload = _json.loads(getattr(result, "stdout", ""))
        assert isinstance(payload, dict)
        return payload

    def test_the_json_summary_carries_the_run(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        output = tmp_path / "hits.parquet"

        result = runner.invoke(
            app,
            [
                "motif-scan",
                str(planted_fasta),
                str(output),
                "--release",
                "2024",
                "--tax-group",
                "all",
                "--workers",
                "1",
                "--json",
            ],
        )

        assert result.exit_code == 0, _output(result)
        summary = self._summary(result)
        assert summary == {
            "release": "2024",
            "tax_group": "all",
            "motifs_scanned": _MOTIFS_LONG_ENOUGH,
            "motifs_skipped": _MOTIFS_TOO_SHORT,
            # Two 600-base records are far under the derivation floor, so 'auto' is uniform
            # — and says so rather than leaving it to be assumed.
            "background": [0.25, 0.25, 0.25, 0.25],
            "threshold": 1e-4,
            "sequences_scanned": 2,
            "hits_written": hit_count(output),
            "workers": 1,
            "output": str(output),
        }

    def test_the_summary_goes_to_stdout_and_the_hits_to_the_named_file(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        # The whole reason the hits are written rather than printed: stdout is one JSON
        # document and nothing else, whatever the scan found.
        output = tmp_path / "hits.parquet"

        result = runner.invoke(
            app, ["motif-scan", str(planted_fasta), str(output), "--workers", "1", "--json"]
        )

        summary = self._summary(result)
        written = read_hits(output)
        assert len(written) > 0
        assert summary["hits_written"] == len(written)
        assert "plantedI" in set(written["sequence_name"])
        # The table's own provenance is what the summary was read off, so the two agree.
        assert written.attrs["release"] == summary["release"]

    def test_the_human_summary_says_what_was_scanned_and_where_it_went(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        output = tmp_path / "hits.parquet"

        result = runner.invoke(
            app, ["motif-scan", str(planted_fasta), str(output), "--workers", "1"]
        )

        assert result.exit_code == 0, _output(result)
        assert f"scanned 2 sequences with {_MOTIFS_LONG_ENOUGH} motifs" in result.stdout
        assert str(output) in result.stdout
        # Which motifs were left out is printed, not silently dropped: an absent factor is
        # explainable only if the scan says it never scanned for it.
        for motif_id in _MOTIFS_TOO_SHORT:
            assert motif_id in result.stdout

    def test_a_scan_that_found_nothing_still_writes_a_file_and_says_so(
        self, motif_release: FakeFetch, tmp_path: Path
    ) -> None:
        empty = tmp_path / "unreadable.fa"
        empty.write_text(">nothing\n" + "N" * 400 + "\n")
        output = tmp_path / "hits.parquet"

        result = runner.invoke(
            app, ["motif-scan", str(empty), str(output), "--workers", "1", "--json"]
        )

        assert result.exit_code == 0, _output(result)
        assert self._summary(result)["hits_written"] == 0
        assert output.is_file()

    def test_a_missing_fasta_exits_non_zero_naming_the_file(
        self, motif_release: FakeFetch, tmp_path: Path
    ) -> None:
        missing = tmp_path / "nowhere.fa"
        output = tmp_path / "hits.parquet"

        result = runner.invoke(
            app, ["motif-scan", str(missing), str(output), "--workers", "1", "--json"]
        )

        assert result.exit_code == 1
        assert "not found" in _output(result)
        assert str(missing) in _output(result)
        # Nothing half-written to be mistaken for an answer.
        assert not output.exists()

    def test_a_release_this_package_does_not_prepare_is_refused_by_name(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        result = runner.invoke(
            app,
            ["motif-scan", str(planted_fasta), str(tmp_path / "hits.parquet"), "--release", "2019"],
        )

        assert result.exit_code == 1
        assert "2024, 2026" in _output(result)

    def test_a_threshold_that_is_not_a_p_value_is_refused(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        result = runner.invoke(
            app,
            ["motif-scan", str(planted_fasta), str(tmp_path / "hits.parquet"), "--threshold", "5"],
        )

        assert result.exit_code == 1
        assert "p-value" in _output(result)

    def test_zero_workers_is_refused_before_anything_is_scanned(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        output = tmp_path / "hits.parquet"

        result = runner.invoke(
            app, ["motif-scan", str(planted_fasta), str(output), "--workers", "0"]
        )

        assert result.exit_code == 1
        assert "at least 1" in _output(result)
        assert not output.exists()

    def test_a_background_mode_reaches_the_scan(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        # The parameter that decides the answer more than any other, and the summary
        # reports the one actually used rather than the one asked for.
        output = tmp_path / "hits.parquet"

        result = runner.invoke(
            app,
            [
                "motif-scan",
                str(planted_fasta),
                str(output),
                "--background",
                "derive",
                "--workers",
                "1",
                "--json",
            ],
        )

        background = self._summary(result)["background"]
        assert background != [0.25, 0.25, 0.25, 0.25]
        assert list(provenance_of(output)["background"]) == background

    def test_a_background_that_is_not_a_mode_is_refused_before_anything_runs(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        output = tmp_path / "hits.parquet"

        result = runner.invoke(
            app, ["motif-scan", str(planted_fasta), str(output), "--background", "gc"]
        )

        assert result.exit_code == 2
        assert not output.exists()
        assert not motif_release.calls  # not even the release was fetched

    def test_the_worker_count_defaults_to_the_allocation(
        self,
        motif_release: FakeFetch,
        planted_fasta: Path,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        # The one place the command deliberately differs from the library, which defaults
        # to serial: a console script is a proper entry point, so it takes what it was
        # given — the allocation first, and never the machine's cores over it.
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
        shared, serial = tmp_path / "shared.parquet", tmp_path / "serial.parquet"

        allocated = runner.invoke(app, ["motif-scan", str(planted_fasta), str(shared), "--json"])
        alone = runner.invoke(
            app, ["motif-scan", str(planted_fasta), str(serial), "--workers", "1", "--json"]
        )

        assert self._summary(allocated)["workers"] == 2
        assert self._summary(alone)["workers"] == 1
        # And the choice is about wall time and nothing else.
        assert read_hits(shared).equals(read_hits(serial))

    def test_the_progress_display_is_suppressed_under_the_json_flag(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        runner.invoke(
            app,
            [
                "motif-scan",
                str(planted_fasta),
                str(tmp_path / "hits.parquet"),
                "--workers",
                "1",
                "--json",
            ],
        )

        assert motif_release.last.progressbar is False

    def test_the_progress_display_is_drawn_without_it(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        runner.invoke(
            app,
            ["motif-scan", str(planted_fasta), str(tmp_path / "hits.parquet"), "--workers", "1"],
        )

        assert motif_release.last.progressbar is True
