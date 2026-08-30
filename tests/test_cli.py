"""Tests for the Typer CLI (``genome``)."""

from __future__ import annotations

import gzip
import hashlib
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
from genome.homology import (
    DEFAULT_RELEASE as HOMOLOGY_RELEASE,
)
from genome.homology import (
    QUALITY_SCORE_COLUMNS,
    HomologySet,
    homology_metadata,
    homology_prepare_command,
    homology_species,
)
from genome.io import download as download_mod
from genome.io import fetch as fetch_mod
from genome.io.completion import RECORD_NAME, read_record, record_path, write_record
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
from genome.tf.cofactor import CofactorTable, cofactor_table
from genome.tf.gene import TFGeneTable, tf_gene_table
from genome.tf.motif import MIN_MOTIF_LENGTH, hit_count, provenance_of, read_hits
from genome.tf.motif import jaspar as jaspar_mod
from genome.tf.motif.jaspar import MOTIF_COUNTS
from genome.xref import (
    ALLIANCE,
    ENSEMBL,
    ENTREZ,
    HGNC,
    MGI,
    UNIPROT,
    XrefSet,
    lookup_xref,
    xref_prepare_command,
    xref_set_dir,
    xref_slice_name,
    xref_species,
    xref_table,
)
from genome.xref import metadata as xref_metadata_mod

from .conftest import CHIMERA_COMPONENTS, COMPONENT_ANNOTATION, FakeFetch
from .test_homology import ABSENT as _NO_HOMOLOG
from .test_homology import FIXTURES as _COMPARA_FIXTURES
from .test_homology import ONE2MANY_HUMAN as _THREE_WORM_HOMOLOGS
from .test_homology import ONE2MANY_WORMS as _THE_THREE_WORMS
from .test_homology import ONE2ONE_HUMAN as _ONE_WORM_HOMOLOG
from .test_homology import ONE2ONE_WORM as _THE_ONE_WORM
from .test_homology import PAIRS as _HOMOLOGY_PAIRS
from .test_homology import _stems as _homology_stems
from .test_jaspar import FIXTURE as _MOTIF_FIXTURE
from .test_jaspar import FIXTURE_COUNT as _MOTIF_COUNT
from .test_jaspar import FIXTURE_MOTIFS as _MOTIF_RECORDS
from .test_scan import FIXTURE as _PLANTED_FASTA
from .test_xref import FIXTURE as _XREF_FIXTURE
from .test_xref import HUMAN_GENE_WITHOUT_A_HUB as _XREF_NO_HUB
from .test_xref import RELEASE as _XREF_RELEASE

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

#: The species ``genome xref`` is exercised against. Human, because the committed Alliance
#: fixture carries the two things the command's shape exists for: a foreign id naming two
#: **Gene id stem**s, and a real gene with no Ensembl cross-reference to reach at all.
_XREF_SPECIES = "Homo sapiens"

#: The three species ``genome homologs`` is exercised against, spelled as the shipped
#: provenance table spells them. All three pairings among them must answer, which is an
#: acceptance criterion of its own rather than a sample.
_HUMAN, _MOUSE, _WORM = "Homo sapiens", "Mus musculus", "Caenorhabditis elegans"

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


def _serve_compara(fake_fetch: FakeFetch, species: str, other: str) -> FakeFetch:
    """Serve the subsample of whichever published dump the shipped table names for a pair.

    Which of the two per-species files holds a pair is exactly what that table records, so
    a test serves what the package asked for rather than deciding for itself — the test
    that serves the *wrong* file does so on purpose and says why.
    """
    row = homology_metadata(species, other, HOMOLOGY_RELEASE)
    assert row is not None
    fake_fetch.serve(_COMPARA_FIXTURES[row.holding_species])
    return fake_fetch


@pytest.fixture
def xref_pinned(monkeypatch: pytest.MonkeyPatch, data_dir: Path) -> None:
    """Pin the curated **Xref source** rows to the committed Alliance fixture's digest.

    The arrangement ``tests/test_xref.py`` uses: the checksum check that holds a truncated
    download to be an error rather than a quietly short answer is never switched off, only
    pointed at what the fake fetch actually serves. Every other cell survives, the real URL
    among them, so what the command prints as provenance is the shipped row's own.
    """
    with gzip.open(data_dir / _XREF_FIXTURE, "rb") as handle:
        # md5 because that is the algorithm Alliance publishes, not a choice made here.
        digest = f"md5:{hashlib.md5(handle.read()).hexdigest()}"
    rows = tuple(replace(row, source_checksum=digest) for row in xref_table())
    monkeypatch.setattr(xref_metadata_mod, "xref_table", lambda: rows)


@pytest.fixture
def xref_release(fake_fetch: FakeFetch, xref_pinned: None) -> FakeFetch:
    """Serve the committed Alliance slice as the publisher's file, pinned to it."""
    fake_fetch.serve(_XREF_FIXTURE)
    return fake_fetch


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


class TestTFCofactorListCommand:
    """``genome tf-cofactor-list`` — an assembly's cofactors, shaped like ``tf-gene-list``.

    The shipped tables answer, as the shipped censuses answer for TF genes: which genes a
    publisher lists as cofactors is that publisher's judgement and no fixture stands in for
    it, so every expectation below is read off the shipped file. What is asserted here is
    the command and not the crossing — ``tests/test_gtf.py`` owns that — so: the
    stdout/stderr split that makes the output pipe, the record ``--json`` emits, and a
    non-zero exit naming the next action for each of the three ways it can fail.

    One test is about neither: a worm assembly answers here while ``tf-gene-list`` refuses
    the same one, because a publisher assessed worm cofactors and none has released a worm
    TF census. That asymmetry is the publishers' shape rather than a defect, and it is
    pinned so that nobody smooths it away.
    """

    def _table(self, assembly: str) -> CofactorTable:
        """The cofactor table shipped for ``assembly``'s species, whichever publisher wrote it."""
        species = metadata.assembly_metadata(assembly).species
        assert species is not None, f"{assembly} has no species in the assembly table"
        table = cofactor_table(species)
        assert table is not None, f"no cofactor table ships for {species}"
        return table

    def _registered(
        self,
        liulab_data: Path,
        *gene_ids: str,
        assembly: str = _MOUSE_ASSEMBLY,
        name: str = _MOUSE_ANNOTATION,
    ) -> None:
        """Register a GTF declaring ``gene_ids`` as ``assembly``'s ``name``, where it lives."""
        source = liulab_data / f"{assembly}.{name}.gtf"
        source.write_text(_gtf_declaring(*gene_ids))
        _register(assembly, liulab_data / "genome" / assembly, source, name)

    def _listed(self, assembly: str, count: int) -> list[str]:
        """``count`` gene ids, one per listed stem, versioned as GENCODE spells them."""
        return [f"{stem}.1" for stem in self._table(assembly).cofactor_stems[:count]]

    def test_only_the_gene_ids_reach_stdout_so_the_output_pipes(self, liulab_data: Path) -> None:
        gene_ids = self._listed(_MOUSE_ASSEMBLY, 2)
        self._registered(liulab_data, *gene_ids)

        result = runner.invoke(app, ["tf-cofactor-list", _MOUSE_ASSEMBLY])

        assert result.exit_code == 0
        assert result.stdout == "".join(f"{gene_id}\n" for gene_id in gene_ids)
        # Whose list it is must be printed and must not cost the pipe, so it goes beside
        # the ids: the heading, the publishers' own attribution, and what the crossing
        # cost — the stems the table holds that this annotation carries no gene for.
        assert f"{_MOUSE_ASSEMBLY} / {_MOUSE_ANNOTATION}" in result.stderr
        assert "Mus musculus" in result.stderr
        assert self._table(_MOUSE_ASSEMBLY).provenance.attribution() in result.stderr
        unresolved = len(self._table(_MOUSE_ASSEMBLY).cofactor_stems) - len(gene_ids)
        assert f"2 cofactors, 2 gene ids, {unresolved} stems" in result.stderr

    def test_the_annotation_may_be_named_instead_of_defaulted(self, liulab_data: Path) -> None:
        gene_ids = self._listed(_MOUSE_ASSEMBLY, 1)
        self._registered(liulab_data, *gene_ids, name="mine")

        result = runner.invoke(app, ["tf-cofactor-list", _MOUSE_ASSEMBLY, "--annotation", "mine"])

        assert result.exit_code == 0
        assert result.stdout.splitlines() == gene_ids
        assert f"{_MOUSE_ASSEMBLY} / mine" in result.stderr

    def test_json_carries_the_cofactors_the_provenance_and_the_unresolved_stems(
        self, liulab_data: Path
    ) -> None:
        table = self._table(_MOUSE_ASSEMBLY)
        stem = table.cofactor_stems[0]
        self._registered(liulab_data, f"{stem}.1")

        result = runner.invoke(app, ["tf-cofactor-list", _MOUSE_ASSEMBLY, "--json"])

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert list(payload) == [
            "assembly",
            "annotation",
            "species",
            "provenance",
            "cofactors",
            "gene_ids",
            "unresolved",
        ]
        assert (payload["assembly"], payload["annotation"]) == (_MOUSE_ASSEMBLY, _MOUSE_ANNOTATION)
        assert payload["gene_ids"] == [f"{stem}.1"]
        # One provenance record per publisher that contributed, never one flattened row.
        assert [source["pubmed_id"] for source in payload["provenance"]["sources"]] == [
            source.pubmed_id for source in table.provenance.sources
        ]
        assert payload["unresolved"]  # what this annotation carries no gene for, not dropped
        cells = dict(zip(table.columns, table.rows[table.gene_id_stems.index(stem)], strict=True))
        cofactor = payload["cofactors"][0]
        assert (cofactor["symbol"], cofactor["source"]) == (cells["symbol"], cells["source"])
        assert cofactor["is_cofactor"] is True
        assert cofactor["classifications"]["animaltfdb_category"] == cells["animaltfdb_category"]

    def test_a_worm_assembly_answers_although_tf_gene_list_refuses_the_same_one(
        self, liulab_data: Path
    ) -> None:
        # AnimalTFDB assessed worm cofactors and no publisher has censused worm
        # transcription factors, so one command answers and the other refuses for one
        # registered annotation. WormBase's ids carry no version, and a stem carrying
        # none resolves to itself, so they are registered exactly as the table spells them.
        gene_ids = list(self._table("ce11").cofactor_stems[:2])
        self._registered(liulab_data, *gene_ids, assembly="ce11", name="wormbase_ws298")

        answered = runner.invoke(app, ["tf-cofactor-list", "ce11"])
        refused = runner.invoke(app, ["tf-gene-list", "ce11"])

        assert answered.exit_code == 0
        assert answered.stdout.splitlines() == gene_ids
        assert "Caenorhabditis elegans" in answered.stderr
        assert refused.exit_code == 1
        assert "no TF census ships" in _output(refused)

    def test_an_unregistered_annotation_exits_one_naming_the_command_that_registers_it(
        self, liulab_data: Path
    ) -> None:
        result = runner.invoke(app, ["tf-cofactor-list", _MOUSE_ASSEMBLY])

        assert result.exit_code == 1
        assert result.stdout == ""
        next_action = f"genome register-annotation {_MOUSE_ASSEMBLY} {_MOUSE_ANNOTATION}"
        assert next_action in _output(result)

    def test_a_species_no_cofactor_table_ships_for_exits_one_naming_the_species_that_have_one(
        self, liulab_data: Path
    ) -> None:
        # Mouse gene ids registered for a yeast assembly: the species is the assembly's own
        # and never what the GTF happens to hold, so this is refused rather than answered.
        self._registered(
            liulab_data, *self._listed(_MOUSE_ASSEMBLY, 1), assembly="sacCer3", name="ensgene_v101"
        )

        result = runner.invoke(app, ["tf-cofactor-list", "sacCer3"])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert "no cofactor table ships" in _output(result)
        assert "Saccharomyces cerevisiae" in _output(result)
        assert "Mus musculus" in _output(result)  # …and what may be asked about instead

    def test_an_assembly_nothing_names_a_species_for_exits_one_saying_so(
        self, liulab_data: Path
    ) -> None:
        # Not the same fact as no table ships: the question was which species this is, and
        # no row answered it. The message says which shipped table could not be chosen, so
        # a cofactor question is never refused with a sentence about a census.
        self._registered(
            liulab_data, *self._listed(_MOUSE_ASSEMBLY, 1), assembly="tiny", name="mine"
        )

        result = runner.invoke(app, ["tf-cofactor-list", "tiny", "--annotation", "mine"])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert "nothing says what species 'tiny' is, so no cofactor table" in _output(result)
        assert "Mus musculus" in _output(result)


class TestXrefCommand:
    """``genome xref`` — the shell surface over an **Xref set**, driven off the fixture.

    Offline throughout, the way ``tests/test_xref.py`` is: the fake fetch serves the
    committed Alliance slice and the curated rows are pinned to it, so the command prepares
    and reads a real set under the test's own data root. What is asserted here is the
    command and not the hop — ``tests/test_xref.py`` owns that — so: the direction being
    named rather than sniffed out of the id strings, the stdout/stderr split that makes the
    output pipe, the ids that resolved to nothing staying visible in *both* renderings, and
    a non-zero exit naming the next action for each way it can fail.

    The strongest claim here is that the command holds no logic the API does not, and it is
    checked the only way that can be: ``--json`` is asserted equal to what the same two
    calls answer in Python, whole, in both directions.
    """

    def test_it_converts_foreign_ids_to_gene_id_stems(self, xref_release: FakeFetch) -> None:
        result = runner.invoke(app, ["xref", _XREF_SPECIES, "--to-stems", ENTREZ, "7157", "672"])

        assert result.exit_code == 0
        assert result.stdout.splitlines() == [
            "7157\tENSG00000141510",
            "672\tENSG00000012048",
        ]

    def test_it_converts_gene_id_stems_to_foreign_ids(self, xref_release: FakeFetch) -> None:
        result = runner.invoke(
            app, ["xref", _XREF_SPECIES, "--from-stems", HGNC, "ENSG00000141510"]
        )

        assert result.exit_code == 0
        assert result.stdout.splitlines() == ["ENSG00000141510\tHGNC:11998"]

    def test_the_direction_is_named_and_never_inferred_from_the_ids(
        self, xref_release: FakeFetch
    ) -> None:
        # One string, one namespace, two directions, two different answers. `HGNC:11998` is
        # an HGNC id and is not a **Gene id stem**, and nothing here works that out from the
        # characters: the flag says which way the hop goes, and asking the wrong way answers
        # *nothing found* rather than quietly turning around.
        toward = runner.invoke(app, ["xref", _XREF_SPECIES, "--to-stems", HGNC, "HGNC:11998"])
        away = runner.invoke(app, ["xref", _XREF_SPECIES, "--from-stems", HGNC, "HGNC:11998"])

        assert toward.exit_code == away.exit_code == 0
        assert toward.stdout.splitlines() == ["HGNC:11998\tENSG00000141510"]
        assert away.stdout.splitlines() == ["HGNC:11998\t"]

    def test_naming_no_direction_exits_two_naming_both_flags(self, xref_release: FakeFetch) -> None:
        result = runner.invoke(app, ["xref", _XREF_SPECIES, "7157"])

        assert result.exit_code == 2
        assert result.stdout == ""
        assert "--to-stems" in _output(result)
        assert "--from-stems" in _output(result)

    def test_naming_both_directions_exits_two(self, xref_release: FakeFetch) -> None:
        result = runner.invoke(
            app, ["xref", _XREF_SPECIES, "--to-stems", ENTREZ, "--from-stems", HGNC, "7157"]
        )

        assert result.exit_code == 2
        assert result.stdout == ""
        assert "exactly one" in _output(result)

    def test_json_is_the_answer_the_api_renders_toward_the_hub(
        self, xref_release: FakeFetch
    ) -> None:
        asked = ["7157", "8086", "999999999"]

        result = runner.invoke(app, ["xref", _XREF_SPECIES, "--to-stems", ENTREZ, *asked, "--json"])

        assert result.exit_code == 0
        assert (
            _json.loads(result.stdout) == XrefSet(_XREF_SPECIES).to_stems(asked, ENTREZ).as_json()
        )

    def test_json_is_the_answer_the_api_renders_away_from_the_hub(
        self, xref_release: FakeFetch
    ) -> None:
        asked = ["ENSG00000141510", "ENSG00000012048", "ENSG00000288541"]

        result = runner.invoke(
            app, ["xref", _XREF_SPECIES, "--from-stems", UNIPROT, *asked, "--json"]
        )

        assert result.exit_code == 0
        assert (
            _json.loads(result.stdout)
            == XrefSet(_XREF_SPECIES).from_stems(asked, UNIPROT).as_json()
        )

    def test_the_pairs_go_to_stdout_and_the_provenance_to_stderr_so_the_output_pipes(
        self, xref_release: FakeFetch
    ) -> None:
        result = runner.invoke(app, ["xref", _XREF_SPECIES, "--to-stems", ENTREZ, "7157"])

        assert result.exit_code == 0
        assert result.stdout == "7157\tENSG00000141510\n"
        # Which publisher said so, and when, is what a reader needs and what a pipeline must
        # not be handed: it goes beside the pairs rather than among them.
        assert f"{ALLIANCE} {_XREF_RELEASE}" in result.stderr
        assert _XREF_SPECIES in result.stderr
        assert lookup_xref(_XREF_SPECIES).url in result.stderr
        assert f"{ENTREZ} ids -> gene id stems" in result.stderr

    def test_ids_that_resolved_to_nothing_stay_in_the_human_output(
        self, xref_release: FakeFetch
    ) -> None:
        # The one thing a hand-rolled join drops. `HGNC:10041` is a real human gene the
        # Alliance lists with no Ensembl cross-reference at all, so it has no hub to reach.
        result = runner.invoke(
            app, ["xref", _XREF_SPECIES, "--to-stems", HGNC, "HGNC:11998", _XREF_NO_HUB]
        )

        assert result.exit_code == 0
        assert result.stdout.splitlines() == [
            "HGNC:11998\tENSG00000141510",
            f"{_XREF_NO_HUB}\t",
        ]
        assert "1 this release names none for" in result.stderr

    def test_ids_that_resolved_to_nothing_stay_in_the_json(self, xref_release: FakeFetch) -> None:
        result = runner.invoke(
            app,
            ["xref", _XREF_SPECIES, "--to-stems", HGNC, "HGNC:11998", _XREF_NO_HUB, "--json"],
        )

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["unresolved"] == [_XREF_NO_HUB]
        assert _XREF_NO_HUB not in payload["resolved"]

    def test_an_id_naming_two_stems_prints_both_and_nothing_picks_one(
        self, xref_release: FakeFetch
    ) -> None:
        # 6.2% of human HGNC ids name two stems in this release, so this is the ordinary
        # case rather than an edge one.
        result = runner.invoke(app, ["xref", _XREF_SPECIES, "--to-stems", ENTREZ, "8086"])

        assert result.exit_code == 0
        assert result.stdout.splitlines() == [
            "8086\tENSG00000094914",
            "8086\tENSG00000291836",
        ]

    def test_the_answer_names_the_source_and_release_that_produced_it(
        self, xref_release: FakeFetch
    ) -> None:
        result = runner.invoke(app, ["xref", _XREF_SPECIES, "--to-stems", ENTREZ, "7157", "--json"])

        payload = _json.loads(result.stdout)
        assert (payload["source"], payload["release"]) == (ALLIANCE, _XREF_RELEASE)
        assert (payload["species"], payload["namespace"]) == (_XREF_SPECIES, ENTREZ)

    def test_a_source_may_be_named(self, xref_release: FakeFetch) -> None:
        result = runner.invoke(
            app, ["xref", _XREF_SPECIES, "--source", ALLIANCE, "--to-stems", ENTREZ, "7157"]
        )

        assert result.exit_code == 0
        assert result.stdout == "7157\tENSG00000141510\n"
        assert f"{ALLIANCE} {_XREF_RELEASE}" in result.stderr

    def test_omitting_the_source_uses_the_default_xref_source(
        self, xref_release: FakeFetch
    ) -> None:
        named = runner.invoke(
            app,
            ["xref", _XREF_SPECIES, "--source", ALLIANCE, "--to-stems", ENTREZ, "7157", "--json"],
        )
        defaulted = runner.invoke(
            app, ["xref", _XREF_SPECIES, "--to-stems", ENTREZ, "7157", "--json"]
        )

        assert named.exit_code == defaulted.exit_code == 0
        assert _json.loads(defaulted.stdout) == _json.loads(named.stdout)
        assert lookup_xref(_XREF_SPECIES).default is True

    def test_a_source_no_set_exists_for_exits_one_naming_the_ones_that_do(
        self, xref_release: FakeFetch
    ) -> None:
        result = runner.invoke(
            app, ["xref", _XREF_SPECIES, "--source", "ensembl", "--to-stems", ENTREZ, "7157"]
        )

        assert result.exit_code == 1
        assert result.stdout == ""
        assert ALLIANCE in _output(result)

    def test_an_unsupported_species_exits_one_naming_the_species_that_have_a_set(
        self, xref_release: FakeFetch
    ) -> None:
        result = runner.invoke(app, ["xref", "Danio rerio", "--to-stems", ENTREZ, "7157"])

        assert result.exit_code == 1
        assert result.stdout == ""
        for species in xref_species():
            assert species in _output(result)
        # Refused before anything was fetched: a species with no Ensembl presence has no hub
        # to hang a namespace off and is unanswerable by design, not pending a download.
        assert xref_release.calls == []

    def test_a_namespace_the_source_does_not_carry_exits_one_naming_the_ones_it_does(
        self, xref_release: FakeFetch
    ) -> None:
        # The three species have three different authorities, so a human set asked for MGI
        # ids fails loudly rather than answering nothing — the failure that would otherwise
        # look like a gene list with no matches.
        result = runner.invoke(app, ["xref", _XREF_SPECIES, "--to-stems", MGI, "MGI:88276"])

        assert result.exit_code == 1
        assert result.stdout == ""
        for namespace in (ENSEMBL, ENTREZ, UNIPROT, HGNC):
            assert namespace in _output(result)

    def test_the_reverse_direction_refuses_the_same_namespace(
        self, xref_release: FakeFetch
    ) -> None:
        result = runner.invoke(app, ["xref", _XREF_SPECIES, "--from-stems", MGI, "ENSG00000141510"])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert HGNC in _output(result)

    def test_a_set_that_is_not_downloaded_exits_one_naming_the_call_for_a_login_node(
        self, monkeypatch: pytest.MonkeyPatch, xref_pinned: None
    ) -> None:
        def no_internet(url: str, dest_dir: Path, **kwargs: object) -> Path:
            raise ConnectionError("the compute node has no internet")

        monkeypatch.setattr(fetch_mod, "fetch_url", no_internet)

        result = runner.invoke(app, ["xref", _XREF_SPECIES, "--to-stems", ENTREZ, "7157"])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert xref_prepare_command(_XREF_SPECIES, ALLIANCE, _XREF_RELEASE) in _output(result)
        assert "login node" in _output(result)

    def test_a_set_left_unfinished_exits_one_naming_the_repair(
        self, xref_release: FakeFetch
    ) -> None:
        directory = xref_set_dir(_XREF_SPECIES, ALLIANCE, _XREF_RELEASE)
        directory.mkdir(parents=True)
        (directory / xref_slice_name(_XREF_SPECIES)).write_bytes(b"half a file")

        result = runner.invoke(app, ["xref", _XREF_SPECIES, "--to-stems", ENTREZ, "7157"])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert xref_prepare_command(_XREF_SPECIES, ALLIANCE, _XREF_RELEASE) in _output(result)

    def test_the_progress_display_is_suppressed_under_json(self, xref_release: FakeFetch) -> None:
        runner.invoke(app, ["xref", _XREF_SPECIES, "--to-stems", ENTREZ, "7157", "--json"])

        assert xref_release.last.progressbar is False

    def test_the_progress_display_is_drawn_without_it(self, xref_release: FakeFetch) -> None:
        runner.invoke(app, ["xref", _XREF_SPECIES, "--to-stems", ENTREZ, "7157"])

        assert xref_release.last.progressbar is True


class TestHomologsCommand:
    """``genome homologs`` — the shell surface over a **Homology set**, on the fixtures.

    Offline throughout, the way ``tests/test_homology.py`` is: the fake fetch serves the
    committed Compara subsamples and the command prepares and reads a real set under the
    test's own data root. What is asserted here is the command and not the set —
    ``tests/test_homology.py`` owns the slice, the partition and the answer — so: all three
    pairings reaching a shell, the publisher's **Homology type** surviving the render, the
    stdout/stderr split that makes the output pipe, the **Dropped partner**s and the null
    quality scores being said out loud, and a non-zero exit naming the next action for each
    way it can fail.

    The strongest claim here is that the command holds no logic the API does not, and it is
    checked the two ways that can be: ``--json`` is asserted equal to what the same call
    answers in Python, whole, and the text rows are asserted to be that JSON's own values
    in its own key order — so nothing the shell sees was assembled here.
    """

    @pytest.mark.parametrize(("species", "other", "links"), _HOMOLOGY_PAIRS)
    def test_it_answers_for_every_pairing_among_human_mouse_and_worm(
        self, fake_fetch: FakeFetch, data_dir: Path, species: str, other: str, links: int
    ) -> None:
        _serve_compara(fake_fetch, species, other)
        asked = _homology_stems(data_dir, species, other)

        result = runner.invoke(app, ["homologs", species, other, *asked])

        assert result.exit_code == 0
        assert len(result.stdout.splitlines()) == links

    def test_a_stem_prints_every_homolog_with_the_publishers_own_type_and_picks_none(
        self, fake_fetch: FakeFetch
    ) -> None:
        _serve_compara(fake_fetch, _HUMAN, _WORM)

        result = runner.invoke(app, ["homologs", _HUMAN, _WORM, _THREE_WORM_HOMOLOGS])

        assert result.exit_code == 0
        rows = [line.split("\t") for line in result.stdout.splitlines()]
        assert [row[1] for row in rows] == list(_THE_THREE_WORMS)
        # Verbatim, and the label the publisher's tree assigned rather than a count of the
        # rows in front of you: three partners here and the type still reads one2many.
        assert {row[2] for row in rows} == {"ortholog_one2many"}

    def test_orthologs_are_the_default_answer(self, fake_fetch: FakeFetch, data_dir: Path) -> None:
        _serve_compara(fake_fetch, _MOUSE, _WORM)
        asked = _homology_stems(data_dir, _MOUSE, _WORM)

        result = runner.invoke(app, ["homologs", _MOUSE, _WORM, *asked])

        assert result.exit_code == 0
        types = {line.split("\t")[2] for line in result.stdout.splitlines()}
        assert types
        assert all(kind.startswith("ortholog_") for kind in types)
        assert "orthologs" in result.stderr

    def test_paralogs_come_back_only_on_request_and_the_heading_says_which_was_asked(
        self, fake_fetch: FakeFetch, data_dir: Path
    ) -> None:
        # Release 116 publishes no cross-species paralogy for these pairs — counted over
        # the whole human dump — so what is asserted is that the flag reaches the API and
        # that the render says which question was asked, not that a paralogy row appeared.
        # Claiming one would claim something about the publisher that is not true.
        _serve_compara(fake_fetch, _HUMAN, _WORM)
        asked = _homology_stems(data_dir, _HUMAN, _WORM)

        result = runner.invoke(app, ["homologs", _HUMAN, _WORM, *asked, "--paralogs", "--json"])
        heading = runner.invoke(app, ["homologs", _HUMAN, _WORM, *asked, "--paralogs"])

        assert result.exit_code == heading.exit_code == 0
        assert (
            _json.loads(result.stdout)
            == HomologySet(_HUMAN, _WORM, progressbar=False)
            .homologs(asked, paralogs=True)
            .as_json()
        )
        assert "paralogy included" in heading.stderr

    def test_a_paralogy_link_would_be_marked_by_the_type_column_it_already_carries(
        self, fake_fetch: FakeFetch, data_dir: Path
    ) -> None:
        # *Not an ortholog* stays distinguishable from *absent* because the publisher's own
        # label is a column of every row rather than a filter applied before printing: a
        # duplication label would print in the same place `ortholog_one2one` does now.
        _serve_compara(fake_fetch, _HUMAN, _WORM)
        asked = _homology_stems(data_dir, _HUMAN, _WORM)

        result = runner.invoke(app, ["homologs", _HUMAN, _WORM, *asked, "--paralogs"])

        assert result.exit_code == 0
        assert {line.split("\t")[2] for line in result.stdout.splitlines()} == {
            "ortholog_one2one",
            "ortholog_one2many",
            "ortholog_many2many",
        }

    def test_json_is_the_answer_the_api_renders(
        self, fake_fetch: FakeFetch, data_dir: Path
    ) -> None:
        _serve_compara(fake_fetch, _MOUSE, _HUMAN)
        asked = [*_homology_stems(data_dir, _MOUSE, _HUMAN), _NO_HOMOLOG]

        result = runner.invoke(app, ["homologs", _MOUSE, _HUMAN, *asked, "--json"])

        assert result.exit_code == 0
        assert (
            _json.loads(result.stdout)
            == HomologySet(_MOUSE, _HUMAN, progressbar=False).homologs(asked).as_json()
        )

    def test_the_links_go_to_stdout_and_the_provenance_to_stderr_so_the_output_pipes(
        self, fake_fetch: FakeFetch
    ) -> None:
        _serve_compara(fake_fetch, _HUMAN, _WORM)
        row = homology_metadata(_HUMAN, _WORM, HOMOLOGY_RELEASE)
        assert row is not None

        result = runner.invoke(app, ["homologs", _HUMAN, _WORM, _ONE_WORM_HOMOLOG])

        assert result.exit_code == 0
        assert result.stdout.startswith(f"{_ONE_WORM_HOMOLOG}\t{_THE_ONE_WORM}\tortholog_one2one\t")
        assert len(result.stdout.splitlines()) == 1
        # Who asserted it, and from which release and file, is what a reader needs and what
        # a pipeline must not be handed: it goes beside the links rather than among them.
        assert row.attribution() in result.stderr
        assert f"{_HUMAN} -> {_WORM}" in result.stderr

    def test_the_text_columns_are_the_json_links_own_values_so_the_two_cannot_drift(
        self, fake_fetch: FakeFetch
    ) -> None:
        # The whole claim that the command holds no logic: every cell printed is a value
        # the API put in the answer, in the order the API writes it, with the publisher's
        # own `NULL` where it recorded nothing.
        _serve_compara(fake_fetch, _MOUSE, _HUMAN)
        asked = "ENSMUSG00000074698"

        text = runner.invoke(app, ["homologs", _MOUSE, _HUMAN, asked])
        rendered = runner.invoke(app, ["homologs", _MOUSE, _HUMAN, asked, "--json"])

        links = _json.loads(rendered.stdout)["resolved"][asked]
        assert text.stdout.splitlines() == [
            "\t".join("NULL" if value is None else str(value) for value in link.values())
            for link in links
        ]

    def test_a_stem_this_release_names_no_homolog_for_stays_in_the_human_output(
        self, fake_fetch: FakeFetch
    ) -> None:
        # The one thing a hand-rolled join drops. A stem with no link gets a row of its
        # own with every other column empty, so nothing leaves shorter than it arrived —
        # and empty is not `NULL`, which would claim a link the publisher scored nothing on.
        _serve_compara(fake_fetch, _HUMAN, _WORM)

        result = runner.invoke(app, ["homologs", _HUMAN, _WORM, _ONE_WORM_HOMOLOG, _NO_HOMOLOG])

        assert result.exit_code == 0
        rows = result.stdout.splitlines()
        assert len(rows) == 2
        assert rows[1].split("\t")[0] == _NO_HOMOLOG
        assert set(rows[1].split("\t")[1:]) == {""}
        assert "1 this release names no homolog for" in result.stderr

    def test_a_stem_that_named_nothing_stays_in_the_json(self, fake_fetch: FakeFetch) -> None:
        _serve_compara(fake_fetch, _HUMAN, _WORM)

        result = runner.invoke(
            app, ["homologs", _HUMAN, _WORM, _ONE_WORM_HOMOLOG, _NO_HOMOLOG, "--json"]
        )

        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["unresolved"] == [_NO_HOMOLOG]
        assert _NO_HOMOLOG not in payload["resolved"]

    def test_a_worm_pairing_names_the_null_quality_scores_before_a_filter_empties(
        self, fake_fetch: FakeFetch
    ) -> None:
        # Compara records neither score on any link of *either* worm pairing, so a shell
        # user about to write `awk -F'\t' '$6 > 50'` is told the column is null throughout
        # rather than discovering it when the filter comes back empty.
        _serve_compara(fake_fetch, _HUMAN, _WORM)

        result = runner.invoke(app, ["homologs", _HUMAN, _WORM, _ONE_WORM_HOMOLOG])
        payload = _json.loads(
            runner.invoke(app, ["homologs", _HUMAN, _WORM, _ONE_WORM_HOMOLOG, "--json"]).stdout
        )

        assert result.exit_code == 0
        # The line that says it, not merely the word appearing somewhere: the column list
        # above it names every column, so a warning nobody wrote would pass a looser check.
        (quality,) = [
            line for line in result.stderr.splitlines() if line.strip().startswith("quality")
        ]
        for column in QUALITY_SCORE_COLUMNS:
            assert column in quality
        assert "empties" in quality
        assert payload["null_quality_scores"] == list(QUALITY_SCORE_COLUMNS)

    def test_the_scored_pairing_says_so_rather_than_staying_silent(
        self, fake_fetch: FakeFetch
    ) -> None:
        # A pair Compara did score says that too: silence would read the same as a warning
        # nobody printed.
        _serve_compara(fake_fetch, _MOUSE, _HUMAN)

        result = runner.invoke(app, ["homologs", _MOUSE, _HUMAN, "ENSMUSG00000074698"])

        assert result.exit_code == 0
        (quality,) = [
            line for line in result.stderr.splitlines() if line.strip().startswith("quality")
        ]
        assert "carry values" in quality
        assert "empties" not in quality

    def test_the_dropped_partners_are_reported_where_an_answer_was_narrowed(
        self, fake_fetch: FakeFetch, data_dir: Path
    ) -> None:
        # Counted *and* named, both off the answer. None are dropped on release 116 — it
        # publishes no cross-species paralogy for the ortholog filter to remove — and zero
        # is printed as an answer rather than left as a silence.
        _serve_compara(fake_fetch, _HUMAN, _WORM)
        asked = _homology_stems(data_dir, _HUMAN, _WORM)

        result = runner.invoke(app, ["homologs", _HUMAN, _WORM, *asked])
        payload = _json.loads(
            runner.invoke(app, ["homologs", _HUMAN, _WORM, *asked, "--json"]).stdout
        )

        assert result.exit_code == 0
        assert "dropped partners" in result.stderr
        assert payload["dropped_partners"] == []

    def test_an_unsupported_species_exits_one_naming_the_species_that_have_a_set(
        self, fake_fetch: FakeFetch
    ) -> None:
        result = runner.invoke(app, ["homologs", "Danio rerio", _HUMAN, _ONE_WORM_HOMOLOG])

        assert result.exit_code == 1
        assert result.stdout == ""
        for species in homology_species():
            assert species in _output(result)
        # Refused before anything was fetched: nobody pinned this species, which must never
        # read as this species having no homologs.
        assert fake_fetch.calls == []

    def test_a_pair_of_one_species_exits_one_because_a_set_relates_two(
        self, fake_fetch: FakeFetch
    ) -> None:
        result = runner.invoke(app, ["homologs", _HUMAN, "homo_sapiens", _ONE_WORM_HOMOLOG])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert "two different species" in _output(result)
        assert fake_fetch.calls == []

    def test_a_release_that_is_not_pinned_exits_one_naming_the_ones_that_are(
        self, fake_fetch: FakeFetch
    ) -> None:
        result = runner.invoke(
            app, ["homologs", _HUMAN, _WORM, _ONE_WORM_HOMOLOG, "--release", "115"]
        )

        assert result.exit_code == 1
        assert result.stdout == ""
        assert HOMOLOGY_RELEASE in _output(result)
        assert fake_fetch.calls == []

    def test_a_versioned_gene_id_exits_one_naming_the_stem_to_pass(
        self, fake_fetch: FakeFetch
    ) -> None:
        # Compara writes its ids bare, so the versioned spelling would match nothing and
        # come back unresolved looking exactly like a gene it never placed in a tree.
        _serve_compara(fake_fetch, _HUMAN, _WORM)

        result = runner.invoke(app, ["homologs", _HUMAN, _WORM, f"{_ONE_WORM_HOMOLOG}.18"])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert _ONE_WORM_HOMOLOG in _output(result)

    def test_a_pair_taken_from_the_wrong_compara_file_exits_one_naming_the_other_file(
        self, fake_fetch: FakeFetch
    ) -> None:
        # The published partition, not a staged one: the real release-116 human dump holds
        # zero human/mouse rows. Serving it for that pair is what a release that had
        # re-partitioned looks like from a shell, and it must be an error naming the other
        # file rather than an empty answer that reads as *these species share no homologs*.
        fake_fetch.serve(_COMPARA_FIXTURES[_HUMAN])

        result = runner.invoke(app, ["homologs", _HUMAN, _MOUSE, "ENSG00000172150"])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert "homo_sapiens/Compara.116.protein_default.homologies.tsv.gz" in _output(result)
        assert "homology_metadata.tsv" in _output(result)

    def test_a_dump_that_is_not_comparas_exits_one_naming_the_columns_it_should_have(
        self, fake_fetch: FakeFetch
    ) -> None:
        fake_fetch.serve("tiny.gtf.gz")

        result = runner.invoke(app, ["homologs", _HUMAN, _MOUSE, "ENSG00000172150"])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert "gene_stable_id" in _output(result)

    def test_a_set_that_is_not_downloaded_exits_one_naming_the_call_for_a_login_node(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        def no_internet(url: str, dest_dir: Path, **kwargs: object) -> Path:
            raise ConnectionError("the compute node has no internet")

        monkeypatch.setattr(fetch_mod, "fetch_url", no_internet)

        result = runner.invoke(app, ["homologs", _HUMAN, _MOUSE, "ENSG00000172150"])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert homology_prepare_command(_HUMAN, _MOUSE, HOMOLOGY_RELEASE) in _output(result)
        assert "login node" in _output(result)

    def test_a_set_left_unfinished_exits_one_naming_the_repair(self, fake_fetch: FakeFetch) -> None:
        _serve_compara(fake_fetch, _HUMAN, _WORM)
        prepared = HomologySet(_HUMAN, _WORM, progressbar=False)
        (prepared.path.parent / RECORD_NAME).unlink()

        result = runner.invoke(app, ["homologs", _HUMAN, _WORM, _ONE_WORM_HOMOLOG])

        assert result.exit_code == 1
        assert result.stdout == ""
        assert "rm -rf" in _output(result)

    def test_the_progress_display_is_suppressed_under_json(self, fake_fetch: FakeFetch) -> None:
        _serve_compara(fake_fetch, _HUMAN, _WORM)

        runner.invoke(app, ["homologs", _HUMAN, _WORM, _ONE_WORM_HOMOLOG, "--json"])

        assert fake_fetch.last.progressbar is False

    def test_the_progress_display_is_drawn_without_it(self, fake_fetch: FakeFetch) -> None:
        _serve_compara(fake_fetch, _HUMAN, _WORM)

        runner.invoke(app, ["homologs", _HUMAN, _WORM, _ONE_WORM_HOMOLOG])

        assert fake_fetch.last.progressbar is True


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
