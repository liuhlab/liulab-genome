"""Tests for the ``genome assembly`` sub-app."""

from __future__ import annotations

import json as _json
import shutil
from dataclasses import dataclass
from pathlib import Path

import pytest

from genome.assembly import download as download_mod
from genome.assembly import metadata
from genome.assembly.fasta import PREPARATION_TOOLS, GenomeFiles
from genome.assembly.metadata import AssemblyMetadata
from genome.cli import app
from genome.store.completion import read_record, record_path, write_record

from .._cli import output, runner
from ..conftest import CHIMERA_COMPONENTS, COMPONENT_ANNOTATION, FakeFetch

_PREPARATION_PRESENT = all(shutil.which(t) is not None for t in PREPARATION_TOOLS)

#: sha256 of the committed ``tiny.fa``, which the fake fetch serves as any assembly.
_TINY_FA_SHA256 = "9316629bab14f9298a043f8b92e1e04a573b12d6a367ccc07c8f8040e5a13981"


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


class _OfflineTinyFasta:
    """Serve the committed ``tiny.fa.gz`` in place of any download, for a whole class.

    Shared by every class below whose tests need nothing more than that one fixture
    file offline — one definition rather than a copy of the same three lines in each.
    """

    @pytest.fixture(autouse=True)
    def _offline(self, fake_fetch: FakeFetch, offline_prepare: None) -> None:
        fake_fetch.serve("tiny.fa.gz")


class TestRegister(_OfflineTinyFasta):
    """``genome assembly register`` — prepare an assembly and say what landed.

    Offline throughout: ``fake_fetch`` serves the committed ``tiny.fa.gz`` in place of
    any download, the ``HEAD`` name check is stubbed, and the shared ``liulab_data``
    fixture puts the assembly directory under this test's own root. The assembly is
    ``tiny``, which no shipped row lists, so nothing is pinned for the fixture to
    disagree with.
    """

    def test_registers_and_reports_where_it_landed_as_text_and_json(
        self, liulab_data: Path
    ) -> None:
        result = runner.invoke(app, ["assembly", "register", "tiny"])
        assert result.exit_code == 0
        assert str(liulab_data / "genome" / "tiny") in result.stdout
        assert _TINY_FA_SHA256 in result.stdout
        assert (liulab_data / "genome" / "tiny" / "tiny.fa").is_file()

        # Already registered: answered from the record, and the same digest either way.
        json_result = runner.invoke(app, ["assembly", "register", "tiny", "--json"])
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert payload["assembly"] == "tiny"
        assert payload["directory"] == str(liulab_data / "genome" / "tiny")
        assert payload["sha256"] == _TINY_FA_SHA256
        assert sorted(payload["files"]) == [
            "tiny.2bit",
            "tiny.chrom.sizes",
            "tiny.fa",
            "tiny.fa.fai",
        ]

    def test_a_broken_directory_is_refused_force_repairs_it_and_a_source_never_asks_ucsc(
        self, liulab_data: Path, data_dir: Path
    ) -> None:
        directory = liulab_data / "genome" / "tiny"
        directory.mkdir(parents=True)
        (directory / "tiny.fa").write_text("half a genome\n")

        refused = runner.invoke(app, ["assembly", "register", "tiny"])
        assert refused.exit_code == 1
        assert "genome assembly register tiny --force" in output(refused)

        repaired = runner.invoke(app, ["assembly", "register", "tiny", "--force", "--json"])
        assert repaired.exit_code == 0
        assert _json.loads(repaired.stdout)["sha256"] == _TINY_FA_SHA256

        sourced = runner.invoke(
            app,
            [
                "assembly",
                "register",
                "tiny_from_source",
                "--source",
                str(data_dir / "tiny.fa.gz"),
                "--json",
            ],
        )
        assert sourced.exit_code == 0
        assert _json.loads(sourced.stdout)["source_url"] == str(data_dir / "tiny.fa.gz")


class TestRegisterResolvesTheName(_OfflineTinyFasta):
    """What a name means, settled by four checks in order — and the two refusals.

    A record already here, then a source the caller named, then a name whose every part
    is prepared here or listed in the shipped table, then the download that was always
    the answer. Offline throughout, and the fetch step is recorded rather than merely
    stubbed: what a refusal is asserted on is that nothing was fetched at all, since
    turning one mistyped string into a whole-genome download per part is the failure the
    gate exists to prevent.
    """

    def test_a_named_source_an_unlisted_download_and_a_bad_name_are_each_resolved_correctly(
        self, data_dir: Path, fake_fetch: FakeFetch
    ) -> None:
        # What a --component flag would have bought, bought for less: the components are
        # in the name, so typing them in the wrong order is detectable — and none of
        # this touches the network.
        mis_ordered = runner.invoke(app, ["assembly", "register", "ecHT115_ce11"])
        assert mis_ordered.exit_code == 1
        assert "`genome assembly register ce11_ecHT115`" in output(mis_ordered)

        missing = runner.invoke(app, ["assembly", "register", "ce11_ecHT115"])
        assert missing.exit_code == 1
        assert "`genome assembly register ce11`" in output(missing)
        assert "`genome assembly register ecHT115`" in output(missing)

        # --force repairs a directory; it does not answer the question of what belongs
        # in one, so it is not a bypass of the same gate.
        forced = runner.invoke(app, ["assembly", "register", "ce11_ecHT115", "--force"])
        assert forced.exit_code == 1
        assert "`genome assembly register ce11`" in output(forced)
        assert fake_fetch.calls == []

        # hg38 and mm10 are both listed, so the name alone reads as two assemblies and
        # would refuse on a machine holding neither. Saying where the bytes come from is
        # the caller answering the question, and it is believed — and still without
        # touching the network.
        source = data_dir / "tiny.fa.gz"
        result = runner.invoke(
            app, ["assembly", "register", "hg38_mm10", "--source", str(source), "--json"]
        )
        assert result.exit_code == 0
        payload = _json.loads(result.stdout)
        assert payload["source_url"] == str(source)
        assert "components" not in payload["details"]
        assert fake_fetch.calls == []

        # The whole separation between ce11_ecHT115 and a free-form local key: neither
        # `my` nor `ref` is an assembly here or in the table, so `my_ref` is one name
        # somebody chose and the download is the answer it always was (ADR-0003).
        unlisted = runner.invoke(app, ["assembly", "register", "my_ref", "--json"])
        assert unlisted.exit_code == 0
        assert fake_fetch.last.url.endswith("my_ref.fa.gz")
        assert "components" not in _json.loads(unlisted.stdout)["details"]

    def test_an_existing_record_is_rebuilt_by_force_until_it_is_lost(
        self, data_dir: Path, liulab_data: Path, fake_fetch: FakeFetch
    ) -> None:
        # The clause that stops a plain hg38_mm10 seeded years ago from silently becoming
        # a chimera: it was registered as an ordinary assembly, so that is what --force
        # registers again — and only once the record is gone does the same directory read
        # as a chimera of hg38 and mm10, neither of which this machine has.
        assert (
            runner.invoke(
                app, ["assembly", "register", "hg38_mm10", "--source", str(data_dir / "tiny.fa.gz")]
            ).exit_code
            == 0
        )

        rebuilt = runner.invoke(app, ["assembly", "register", "hg38_mm10", "--force", "--json"])
        assert rebuilt.exit_code == 0
        assert "components" not in _json.loads(rebuilt.stdout)["details"]
        assert fake_fetch.last.url.endswith("hg38_mm10.fa.gz")

        record_path(liulab_data / "genome" / "hg38_mm10").unlink()
        lost = runner.invoke(app, ["assembly", "register", "hg38_mm10", "--force"])
        assert lost.exit_code == 1
        assert "`genome assembly register hg38`" in output(lost)
        assert "`genome assembly register mm10`" in output(lost)


class TestVerify(_OfflineTinyFasta):
    """``genome assembly verify`` — re-read a FASTA and check it against the official row."""

    def test_reports_the_digest_as_text_and_json_a_mismatch_or_nothing_registered(
        self, data_dir: Path
    ) -> None:
        # sacCer3's row pins the real genome's digest; the fixture is a subsample of it,
        # so this is the mismatch a copy from a bad mirror would produce.
        mismatch = runner.invoke(
            app, ["assembly", "verify", "sacCer3", "--fasta", str(data_dir / "tiny.fa")]
        )
        assert mismatch.exit_code == 1
        assert "sha256 mismatch" in output(mismatch)

        unregistered = runner.invoke(app, ["assembly", "verify", "tiny"])
        assert unregistered.exit_code == 1
        assert "genome assembly register tiny" in output(unregistered)

        assert runner.invoke(app, ["assembly", "register", "tiny"]).exit_code == 0

        result = runner.invoke(app, ["assembly", "verify", "tiny"])
        assert result.exit_code == 0
        assert _TINY_FA_SHA256 in result.stdout

        json_result = runner.invoke(app, ["assembly", "verify", "tiny", "--json"])
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert payload["sha256"] == _TINY_FA_SHA256
        # No row lists "tiny", so what it is held to is the digest its own registration
        # recorded — the fallback, and the payload says which answered.
        assert payload["expected"] == _TINY_FA_SHA256
        assert payload["expected_from"] == "record"
        assert payload["verified"] is True
        # An assembly that is not a chimera has no components to be asked about — null,
        # rather than a status that would read as a check somebody made.
        assert payload["components"] is None


class TestWhatAVerifiedDigestWasHeldTo(_OfflineTinyFasta):
    """Three answers, three sentences — and never one wording covering two of them.

    Being held to the digest the lab pinned, being held to the one this machine last
    produced, and being held to nothing at all are different results, and a caller who
    cannot tell them apart reads the weakest as the strongest.
    """

    def test_the_record_the_row_and_nothing_at_all_each_say_a_different_sentence(
        self, monkeypatch: pytest.MonkeyPatch, data_dir: Path
    ) -> None:
        assert runner.invoke(app, ["assembly", "register", "tiny"]).exit_code == 0

        record_pinned = runner.invoke(app, ["assembly", "verify", "tiny"])
        assert record_pinned.exit_code == 0
        assert "own registration recorded" in record_pinned.stdout
        assert "not an independent pin" in record_pinned.stdout

        # Once the metadata table pins its own digest for the same assembly, that row is
        # what a re-run is held to — a stronger claim than the record it already trusted.
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
        row_pinned = runner.invoke(app, ["assembly", "verify", "tiny"])
        assert row_pinned.exit_code == 0
        assert "matches the digest the metadata table pins for it" in row_pinned.stdout

        # A FASTA handed over by hand, checked against an assembly whose row pins nothing
        # and which is not registered here: there is no digest to be held to at all.
        nothing_pinned = runner.invoke(
            app,
            ["assembly", "verify", "ce11_ecHT115", "--fasta", str(data_dir / "tiny.fa"), "--json"],
        )
        assert nothing_pinned.exit_code == 0
        payload = _json.loads(nothing_pinned.stdout)
        assert (payload["expected"], payload["expected_from"], payload["verified"]) == (
            None,
            None,
            False,
        )
        # Nothing is asked about components either: the assembly's own registration is
        # not what is being verified.
        assert payload["components"] is None

        nothing_pinned_human = runner.invoke(
            app, ["assembly", "verify", "ce11_ecHT115", "--fasta", str(data_dir / "tiny.fa")]
        )
        assert "nothing to check it against" in nothing_pinned_human.stdout
        assert "components" not in nothing_pinned_human.stdout


class TestTableRow:
    """``genome assembly table-row`` — download an assembly and print its finished table row.

    Offline throughout: ``fake_fetch`` serves the committed ``tiny.fa.gz`` in place of
    any download, and the shared ``liulab_data`` fixture puts the assembly directory
    under this test's own root. hg38 and sacCer3 are used because the shipped table pins a
    source URL for both, which also skips the network name check.
    """

    @pytest.fixture(autouse=True)
    def _offline(self, fake_fetch: FakeFetch) -> None:
        fake_fetch.serve("tiny.fa.gz")

    def test_prints_the_row_as_text_and_json_and_reports_an_existing_pin(self) -> None:
        result = runner.invoke(app, ["assembly", "table-row", "hg38"])
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

        json_result = runner.invoke(app, ["assembly", "table-row", "hg38", "--json"])
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert payload["assembly_name"] == "hg38"
        assert payload["ncbi_taxid"] == 9606
        assert payload["sha256"] == _TINY_FA_SHA256

        # sacCer3's row already pins the real genome's digest, and the fixture is a
        # subsample of it, so the two disagree. This is the command a maintainer runs
        # precisely when an upstream file has changed and the pin must be regenerated,
        # so it prints what actually arrived instead of refusing.
        mismatched = runner.invoke(app, ["assembly", "table-row", "sacCer3"])
        assert mismatched.exit_code == 0
        mismatched_row = mismatched.stdout.strip().split("\t")
        assert mismatched_row[0] == "sacCer3"
        assert mismatched_row[metadata.METADATA_FIELDS.index("sha256")] == _TINY_FA_SHA256

    def test_a_chimera_is_refused_before_anything_is_downloaded(
        self, fake_fetch: FakeFetch
    ) -> None:
        # A chimera pins nothing, so this command has no job to do for one. The refusal
        # describes the row — the name, and every other column blank — rather than
        # printing a line that would look like one it computed something for.
        result = runner.invoke(app, ["assembly", "table-row", "ce11_ecHT115"])

        assert result.exit_code == 1
        assert "ce11, ecHT115" in output(result)
        assert "no sha256" in output(result)
        assert "genome assembly verify ce11_ecHT115" in output(result)
        assert fake_fetch.calls == []


#: The repair every chimera error names. Quoted here so the tests below can assert that a
#: message carries it and then run exactly it — a message naming a command nobody can
#: follow is worse than no message.
_CHIMERA_REPAIR = "genome assembly register tinyCe_tinySc --force"


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
    """``genome assembly register <name>`` is the only build spelling, end to end.

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
            seeded = runner.invoke(
                app, ["assembly", "register", name, "--source", str(component.fasta)]
            )
            assert seeded.exit_code == 0, output(seeded)
            if annotate:
                gtf = component.gtf
                assert gtf is not None
                built = runner.invoke(
                    app, ["annotation", "register-gtf", name, str(gtf), COMPONENT_ANNOTATION]
                )
                assert built.exit_code == 0, output(built)

    def test_naming_a_chimera_builds_it_and_the_report_says_what_it_is_made_of(
        self, liulab_data: Path
    ) -> None:
        # The line that used to read `source None` is the components, and the closing one
        # names the annotation this same command registered.
        self._register_components("tinyCe", "tinySc", annotate=True)

        result = runner.invoke(app, ["assembly", "register", "tinyCe_tinySc"])

        assert result.exit_code == 0, output(result)
        assert "  components  tinyCe, tinySc" in result.stdout
        assert "source" not in result.stdout
        assert f"  annotation  {COMPONENT_ANNOTATION}+{COMPONENT_ANNOTATION}" in result.stdout
        assert (liulab_data / "genome" / "tinyCe_tinySc" / "tinyCe_tinySc.fa").is_file()

    def test_a_build_with_nothing_to_merge_the_json_payload_and_verify_all_agree(
        self, liulab_data: Path
    ) -> None:
        self._register_components("tinyCe", "tinySc")

        result = runner.invoke(app, ["assembly", "register", "tinyCe_tinySc"])
        assert result.exit_code == 0, output(result)
        assert "  annotation  none" in result.stdout

        # It already carried the components, which is why nothing was added to it.
        json_result = runner.invoke(app, ["assembly", "register", "tinyCe_tinySc", "--json"])
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert payload["source_url"] is None
        assert [entry["name"] for entry in payload["details"]["components"]] == [
            "tinyCe",
            "tinySc",
        ]

        verify_json = runner.invoke(app, ["assembly", "verify", "tinyCe_tinySc", "--json"])
        verify_payload = _json.loads(verify_json.stdout)
        verify_human = runner.invoke(app, ["assembly", "verify", "tinyCe_tinySc"])
        assert verify_payload["components"] == "unchanged"
        assert "components  unchanged" in verify_human.stdout

        # The line prints either way: a chimera whose components could not be compared
        # is unproven, and silence would be exactly what a pass looks like.
        directory = liulab_data / "genome" / "tinyCe_tinySc"
        record = read_record(directory)
        assert record is not None
        for entry in record.details["components"]:
            entry["sha256"] = None
        write_record(directory, record)

        unknown_json = runner.invoke(app, ["assembly", "verify", "tinyCe_tinySc", "--json"])
        unknown_payload = _json.loads(unknown_json.stdout)
        unknown_human = runner.invoke(app, ["assembly", "verify", "tinyCe_tinySc"])
        assert unknown_payload["components"] == "unknown"
        assert unknown_human.exit_code == 0
        assert "components  unknown" in unknown_human.stdout

    def test_the_named_repair_rebuilds_a_mismatched_or_a_lost_record_by_name(
        self, tmp_path: Path, liulab_data: Path
    ) -> None:
        # Run verbatim, not paraphrased: this command used to route to the downloader and
        # fail with "Unknown UCSC assembly", so every chimera error quoted a repair nobody
        # could follow. And the hole this closes: opening by name used to return from the
        # chimera's own record, which vouches for its files and can say nothing about the
        # ones they were copied from — only building and verifying used to ask.
        self._register_components("tinyCe", "tinySc")
        assert runner.invoke(app, ["assembly", "register", "tinyCe_tinySc"]).exit_code == 0
        corrected = str(_corrected_component(tmp_path / "corrected.fa"))
        assert (
            runner.invoke(
                app, ["assembly", "register", "tinySc", "--force", "--source", corrected]
            ).exit_code
            == 0
        )
        refused = runner.invoke(app, ["assembly", "register", "tinyCe_tinySc"])
        assert refused.exit_code == 1
        assert _CHIMERA_REPAIR in output(refused)

        repaired = runner.invoke(app, _CHIMERA_REPAIR.split()[1:])
        assert repaired.exit_code == 0, output(repaired)
        # Rebuilt, not merely re-recorded: the corrected component's bases are in it.
        fasta = (liulab_data / "genome" / "tinyCe_tinySc" / "tinyCe_tinySc.fa").read_text()
        assert "ACGTACGTAC" in fasta
        verified = runner.invoke(app, ["assembly", "verify", "tinyCe_tinySc", "--json"])
        assert _json.loads(verified.stdout)["components"] == "unchanged"

        # The residual a lost record leaves: the name is the only surviving information
        # about what this directory was, and it is enough.
        record_path(liulab_data / "genome" / "tinyCe_tinySc").unlink()
        refused_again = runner.invoke(app, ["assembly", "register", "tinyCe_tinySc"])
        assert refused_again.exit_code == 1
        assert _CHIMERA_REPAIR in output(refused_again)

        repaired_again = runner.invoke(app, _CHIMERA_REPAIR.split()[1:])
        assert repaired_again.exit_code == 0, output(repaired_again)
        assert "  components  tinyCe, tinySc" in repaired_again.stdout
