"""Tests for the ``genome annotation`` sub-app."""

from __future__ import annotations

import json as _json
import shutil
from dataclasses import replace
from pathlib import Path

import gffutils
import pytest

from genome.annotation import (
    AnnotationRegistry,
    GtfAnnotation,
    MergeSource,
    annotation_dir,
    register_merged_gtf,
)
from genome.annotation import metadata as annotation_metadata
from genome.annotation.curated import curated_gene_list
from genome.annotation.metadata import AnnotationMetadata
from genome.cli import app
from genome.store.completion import record_path

from .._cli import output, runner
from ..conftest import FakeFetch

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


def _register(assembly: str, assembly_dir: Path, gtf: Path, name: str) -> GtfAnnotation:
    """Register ``gtf`` under ``assembly_dir``, so a command has something to report on."""
    return AnnotationRegistry.locate(assembly, assembly_dir).register_path(gtf, name)


def _feature_types(database_path: Path) -> list[str]:
    """The kinds of feature a built database holds, with the connection closed behind us."""
    database = gffutils.FeatureDB(str(database_path))
    try:
        return sorted(database.featuretypes())
    finally:
        database.conn.close()


class _OfflineTinyGtf:
    """Serve the committed ``tiny.gtf.gz`` as the one row the table offers, for a class.

    The two classes below that register an annotation *by name* need exactly this and
    nothing else — the CLI takes no metadata argument by design, so the table it reads is
    what moves — and they had a copy of it each. The classes registering by path need
    neither half and do not inherit it.
    """

    @pytest.fixture(autouse=True)
    def _offline(self, fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch) -> None:
        fake_fetch.serve("tiny.gtf.gz")
        monkeypatch.setattr(annotation_metadata, "annotation_table", lambda: (_TINY_ANNOTATION,))


class TestRegisterAnnotation(_OfflineTinyGtf):
    """``genome annotation register`` — fetch, verify and build an annotation by name.

    The CLI is a thin client over the shipped table and takes no metadata argument by
    design, so the table itself is what is stood in for here: one row pointing at the
    committed ``tiny.gtf.gz`` and pinning its digest. The registration it drives — the
    fetch, the checksum, the real gffutils build, the record — is the shipped code.
    """

    def test_a_broken_directory_is_refused_force_repairs_it_and_the_repair_reports_correctly(
        self, liulab_data: Path
    ) -> None:
        directory = liulab_data / "genome" / "tiny" / "gtf" / "ensgene_v101"
        directory.mkdir(parents=True)
        (directory / "ensgene_v101.db").write_bytes(b"half a database")

        refused = runner.invoke(app, ["annotation", "register", "tiny", "ensgene_v101"])
        assert refused.exit_code == 1
        assert "genome annotation register tiny ensgene_v101 --force" in output(refused)

        repaired = runner.invoke(
            app, ["annotation", "register", "tiny", "ensgene_v101", "--force", "--json"]
        )
        assert repaired.exit_code == 0
        assert _json.loads(repaired.stdout)["sha256"] == _TINY_GTF_SHA256

        # Already registered by the repair above: answered from the same record, text
        # and json alike.
        result = runner.invoke(app, ["annotation", "register", "tiny", "ensgene_v101"])
        assert result.exit_code == 0
        assert str(directory) in result.stdout
        assert _TINY_GTF_SHA256 in result.stdout
        assert (directory / "ensgene_v101.gtf").is_file()
        assert (directory / "ensgene_v101.db").is_file()

        json_result = runner.invoke(
            app, ["annotation", "register", "tiny", "ensgene_v101", "--json"]
        )
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert payload["assembly"] == "tiny"
        assert payload["name"] == "ensgene_v101"
        assert payload["directory"] == str(directory)
        assert payload["source_url"] == _ANNOTATION_URL
        assert payload["sha256"] == _TINY_GTF_SHA256
        assert sorted(payload["files"]) == ["ensgene_v101.db", "ensgene_v101.gtf"]

        unlisted = runner.invoke(app, ["annotation", "register", "tiny", "nope"])
        assert unlisted.exit_code == 1
        assert "ensgene_v101" in output(unlisted)
        # …and the command that registers a GTF the table does not list, since that is
        # what a caller who named an unlisted annotation is most likely reaching for.
        assert "genome annotation register-gtf tiny" in output(unlisted)

    def test_the_chromosome_check_is_stood_down_and_feature_inference_is_reachable(
        self,
        fake_fetch: FakeFetch,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        liulab_data: Path,
    ) -> None:
        # The committed Ensembl-spelled GTF (I, II, III) against a UCSC-spelled
        # assembly (chrI, chrII, chrIII): refused by default, and registered anyway
        # once the caller says they have looked at the mismatch and accept it.
        fake_fetch.serve("ensembl_style.gtf")
        assembly_dir = liulab_data / "genome" / "tiny"
        assembly_dir.mkdir(parents=True)
        (assembly_dir / "tiny.chrom.sizes").write_text("chrI\t10000\nchrII\t10000\nchrIII\t10000\n")
        ensembl_row = replace(
            _TINY_ANNOTATION,
            url="https://mirror.example.invalid/annotations/ensembl_style.gtf",
            sha256=None,
        )
        monkeypatch.setattr(annotation_metadata, "annotation_table", lambda: (ensembl_row,))

        refused = runner.invoke(app, ["annotation", "register", "tiny", "ensgene_v101"])
        assert refused.exit_code == 1
        assert "chromosome" in output(refused)

        stood_down = runner.invoke(
            app,
            ["annotation", "register", "tiny", "ensgene_v101", "--no-check-chromosomes", "--json"],
        )
        assert stood_down.exit_code == 0
        details = _json.loads(stood_down.stdout)["details"]
        assert details["chromosomes_checked"] is False
        # The reason rides in `details` as the record holds it — no second spelling of it
        # for the JSON surface to drift from.
        assert details["chromosomes_unchecked_because"] == "caller-override"

        # Nothing says a listed annotation declares genes and transcripts, so the
        # inference the API exposes has to be reachable on this command as well.
        bare = tmp_path / "bare.gtf"
        bare.write_text(_BARE_GTF)
        fake_fetch.serve(bare)
        bare_row = replace(
            _TINY_ANNOTATION,
            name="bare",
            provider="somebody",
            version="1",
            url="https://mirror.example.invalid/annotations/bare.gtf",
            sha256=None,
        )
        monkeypatch.setattr(annotation_metadata, "annotation_table", lambda: (bare_row,))

        inferred = runner.invoke(
            app, ["annotation", "register", "tiny", "bare", "--infer-genes", "--infer-transcripts"]
        )
        assert inferred.exit_code == 0
        database = liulab_data / "genome" / "tiny" / "gtf" / "bare" / "bare.db"
        assert _feature_types(database) == ["exon", "gene", "transcript"]


class TestRegisterGtf:
    """``genome annotation register-gtf`` — register a GTF the annotation table does not list.

    The by-path way in, from a shell: no table row, no download, no checksum to compare
    against — the caller says where the file is. Offline by construction, since the GTF
    is a local one, and the shared ``liulab_data`` fixture puts the assembly directory
    under this test's own root.
    """

    def test_registers_reports_as_text_json_is_listed_and_a_missing_gtf_is_refused(
        self, tmp_path: Path, data_dir: Path, liulab_data: Path
    ) -> None:
        source = data_dir / "tiny.gtf"
        directory = liulab_data / "genome" / "tiny" / "gtf" / "mine"

        result = runner.invoke(app, ["annotation", "register-gtf", "tiny", str(source), "mine"])
        assert result.exit_code == 0
        assert str(directory) in result.stdout
        assert str(source) in result.stdout
        assert (directory / "mine.gtf").is_file()
        assert (directory / "mine.db").is_file()

        json_result = runner.invoke(
            app, ["annotation", "register-gtf", "tiny", str(source), "mine", "--json"]
        )
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert payload["assembly"] == "tiny"
        assert payload["name"] == "mine"
        assert payload["directory"] == str(directory)
        assert payload["source_url"] == str(source)
        assert payload["sha256"] == _TINY_GTF_SHA256
        assert sorted(payload["files"]) == ["mine.db", "mine.gtf"]

        listed = runner.invoke(app, ["annotation", "list", "tiny", "--json"])
        assert listed.exit_code == 0
        rows = _json.loads(listed.stdout)["annotations"]
        assert [(row["name"], row["offered"], row["registered"]) for row in rows] == [
            ("mine", False, True)
        ]

        missing = runner.invoke(
            app, ["annotation", "register-gtf", "tiny", str(tmp_path / "nope.gtf"), "mine"]
        )
        assert missing.exit_code == 1
        assert "GTF file not found" in output(missing)

    def test_a_broken_directory_is_refused_and_force_repairs_it(
        self, data_dir: Path, liulab_data: Path
    ) -> None:
        source = data_dir / "tiny.gtf"
        directory = liulab_data / "genome" / "tiny" / "gtf" / "broken"
        directory.mkdir(parents=True)
        (directory / "broken.db").write_bytes(b"half a database")

        refused = runner.invoke(app, ["annotation", "register-gtf", "tiny", str(source), "broken"])
        assert refused.exit_code == 1
        assert f"genome annotation register-gtf tiny {source} broken --force" in output(refused)

        repaired = runner.invoke(
            app, ["annotation", "register-gtf", "tiny", str(source), "broken", "--force", "--json"]
        )
        assert repaired.exit_code == 0
        assert _json.loads(repaired.stdout)["sha256"] == _TINY_GTF_SHA256

    def test_the_chromosome_check_is_stood_down_and_feature_inference_is_reachable(
        self, tmp_path: Path, data_dir: Path, liulab_data: Path
    ) -> None:
        # The committed Ensembl-spelled GTF (I, II, III) against a UCSC-spelled assembly
        # (chrI, chrII, chrIII): the assembly's chrom.sizes is found from its name, so
        # this way in checks the names too — and stands the check down when asked.
        source = data_dir / "ensembl_style.gtf"
        assembly_dir = liulab_data / "genome" / "tiny"
        assembly_dir.mkdir(parents=True)
        (assembly_dir / "tiny.chrom.sizes").write_text("chrI\t10000\nchrII\t10000\nchrIII\t10000\n")

        refused = runner.invoke(app, ["annotation", "register-gtf", "tiny", str(source), "mine"])
        assert refused.exit_code == 1
        assert "chromosome" in output(refused)

        stood_down = runner.invoke(
            app,
            [
                "annotation",
                "register-gtf",
                "tiny",
                str(source),
                "mine",
                "--no-check-chromosomes",
                "--json",
            ],
        )
        assert stood_down.exit_code == 0
        details = _json.loads(stood_down.stdout)["details"]
        assert details["chromosomes_checked"] is False
        assert details["chromosomes_unchecked_because"] == "caller-override"

        # Without the flags the database holds exons and nothing else — genes and
        # transcripts are what a caller registers an annotation for.
        bare = tmp_path / "bare.gtf"
        bare.write_text(_BARE_GTF)
        gtf_root = liulab_data / "genome" / "tiny" / "gtf"
        assert (
            runner.invoke(app, ["annotation", "register-gtf", "tiny", str(bare), "exons"]).exit_code
            == 0
        )
        assert _feature_types(gtf_root / "exons" / "exons.db") == ["exon"]

        inferred = runner.invoke(
            app,
            [
                "annotation",
                "register-gtf",
                "tiny",
                str(bare),
                "genes",
                "--infer-genes",
                "--infer-transcripts",
            ],
        )
        assert inferred.exit_code == 0
        assert _feature_types(gtf_root / "genes" / "genes.db") == ["exon", "gene", "transcript"]


class TestWhatARegistrationSaysAboutTheChromosomes(_OfflineTinyGtf):
    """Both registration commands say which of four things happened to the name check.

    ``--no-check-chromosomes`` used to be answered with "register the assembly first",
    which the caller may well have done already: the record could not tell *nothing to
    check against* from *the caller stood the check down*, so the surface picked one and
    was wrong half the time. Each state now has its own sentence, and the one that had
    advice to give is the only one that gives any.
    """

    #: The advice that belongs to exactly one of the four states.
    _ADVICE = "register the assembly first"

    @staticmethod
    def _prepare_assembly(liulab_data: Path) -> None:
        """Put the assembly's ``chrom.sizes`` where the check looks for it."""
        assembly_dir = liulab_data / "genome" / "tiny"
        assembly_dir.mkdir(parents=True, exist_ok=True)
        (assembly_dir / "tiny.chrom.sizes").write_text("chrI\t10000\nchrII\t10000\nchrIII\t10000\n")

    def test_a_check_that_ran_or_had_nothing_to_run_against_says_which(
        self, data_dir: Path, liulab_data: Path
    ) -> None:
        # No chrom.sizes yet: the assembly is not registered here, and registering it is
        # exactly what would let the names be verified — the one state that advises.
        nothing_to_check = runner.invoke(app, ["annotation", "register", "tiny", "ensgene_v101"])
        assert nothing_to_check.exit_code == 0
        assert "chromosomes not checked" in nothing_to_check.stdout
        assert self._ADVICE in nothing_to_check.stdout
        shutil.rmtree(liulab_data / "genome" / "tiny" / "gtf" / "ensgene_v101")

        # A check that ran is reported by both commands rather than left to silence,
        # which would read exactly the same as a surface that had nothing good to say.
        self._prepare_assembly(liulab_data)
        by_name = runner.invoke(app, ["annotation", "register", "tiny", "ensgene_v101"])
        by_path = runner.invoke(
            app, ["annotation", "register-gtf", "tiny", str(data_dir / "tiny.gtf"), "mine"]
        )
        for result in (by_name, by_path):
            assert result.exit_code == 0
            assert "chromosomes checked" in result.stdout
            assert self._ADVICE not in result.stdout

    def test_standing_the_check_down_or_an_old_record_is_advised_nothing_wrong(
        self, data_dir: Path, liulab_data: Path
    ) -> None:
        # The bug this fixes: the assembly is registered, the caller turned the check off
        # deliberately, and being told to register the assembly first is wrong.
        self._prepare_assembly(liulab_data)
        by_name = runner.invoke(
            app, ["annotation", "register", "tiny", "ensgene_v101", "--no-check-chromosomes"]
        )
        by_path = runner.invoke(
            app,
            [
                "annotation",
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

        # A record written before the reason existed reports it as unknown: the record
        # returned is the one already on disk, whose bare `false` stands for either
        # reason. Neither may be claimed, and neither raises.
        path = record_path(annotation_dir(liulab_data / "genome" / "tiny", "ensgene_v101"))
        written = _json.loads(path.read_text())
        written["details"] = {"chromosomes_checked": False}
        path.write_text(_json.dumps(written))

        unknown = runner.invoke(app, ["annotation", "register", "tiny", "ensgene_v101"])
        assert unknown.exit_code == 0
        assert "does not say why" in unknown.stdout
        assert self._ADVICE not in unknown.stdout
        assert "stood down" not in unknown.stdout


class TestAnnotations:
    """``genome annotation list`` — what the tables offer, set against what is registered here.

    The shipped table answers here rather than one stood up for the test: reporting it
    is this command's whole job. hg38 is the assembly, whose row offers one annotation
    and flags it as the default.
    """

    def test_nothing_registered_reports_as_text_and_json_and_an_unlisted_assembly_is_empty(
        self, liulab_data: Path
    ) -> None:
        result = runner.invoke(app, ["annotation", "list", "hg38"])
        assert result.exit_code == 0
        assert "gencode_v50" in result.stdout
        assert "offered, not registered" in result.stdout
        assert "genome annotation register hg38 gencode_v50" in result.stdout
        # Nothing was prepared to answer the question — the assembly is not even there.
        assert not (liulab_data / "genome" / "hg38").exists()

        json_result = runner.invoke(app, ["annotation", "list", "hg38", "--json"])
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert payload["assembly"] == "hg38"
        assert payload["directory"] == str(liulab_data / "genome" / "hg38")
        assert payload["default_annotation"] == "gencode_v50"
        assert [
            (row["name"], row["offered"], row["registered"]) for row in payload["annotations"]
        ] == [("gencode_v50", True, False)]

        unlisted = runner.invoke(app, ["annotation", "list", "tiny"])
        assert unlisted.exit_code == 0
        assert "tiny" in unlisted.stdout

    def test_it_sets_registered_against_offered_and_a_broken_offered_one_names_its_repair(
        self, data_dir: Path, liulab_data: Path
    ) -> None:
        assembly_dir = liulab_data / "genome" / "hg38"
        _register("hg38", assembly_dir, data_dir / "tiny.gtf", "mine")

        result = runner.invoke(app, ["annotation", "list", "hg38", "--json"])
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

        # It used to read `offered, not registered` — indistinguishable from one nobody
        # had ever fetched — and the closing line sent the reader to a command that
        # would itself raise and demand --force.
        _register("hg38", assembly_dir, data_dir / "tiny.gtf", "gencode_v50")
        record_path(annotation_dir(assembly_dir, "gencode_v50")).unlink()

        broken = runner.invoke(app, ["annotation", "list", "hg38"])
        assert broken.exit_code == 0
        assert "offered, not registered" not in broken.stdout
        assert "broken" in broken.stdout
        assert "genome annotation register hg38 gencode_v50 --force" in broken.stdout
        default_line = next(
            line for line in broken.stdout.splitlines() if line.startswith("default:")
        )
        assert "--force" in default_line

    def test_a_broken_unlisted_annotation_is_listed_and_json_carries_its_broken_state(
        self, data_dir: Path, liulab_data: Path
    ) -> None:
        assembly_dir = liulab_data / "genome" / "hg38"
        _register("hg38", assembly_dir, data_dir / "tiny.gtf", "healthy")
        annotation = _register("hg38", assembly_dir, data_dir / "tiny.gtf", "mine")
        annotation.db.write_bytes(b"truncated")

        result = runner.invoke(app, ["annotation", "list", "hg38"])
        assert result.exit_code == 0
        assert "mine" in result.stdout
        assert "broken" in result.stdout
        assert (
            f"genome annotation register-gtf hg38 {data_dir / 'tiny.gtf'} mine --force"
            in result.stdout
        )

        json_result = runner.invoke(app, ["annotation", "list", "hg38", "--json"])
        assert json_result.exit_code == 0
        rows = {row["name"]: row for row in _json.loads(json_result.stdout)["annotations"]}
        assert rows["mine"]["broken"] is True
        assert rows["mine"]["registered"] is False
        assert rows["mine"]["repair"].endswith("mine --force")
        assert "mine.db" in rows["mine"]["problem"]
        # One broken annotation costs neither the exit code nor the ones beside it.
        assert rows["healthy"]["broken"] is False
        assert rows["healthy"]["registered"] is True
        assert rows["gencode_v50"]["broken"] is False


class TestGeneCategoryCommands:
    """``genome annotation gene-list`` and ``genome annotation gene-categories``.

    A category's genes, and which categories exist.

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

    def test_only_the_gene_ids_reach_stdout_annotation_may_be_named_and_json_keeps_sources_apart(
        self, liulab_data: Path, data_dir: Path
    ) -> None:
        self._registered(liulab_data, data_dir)
        category = self._declared("ensgene_v101")[0]

        result = runner.invoke(app, ["annotation", "gene-list", "sacCer3", category])
        assert result.exit_code == 0
        assert result.stdout == "".join(
            f"{gene_id}\n" for gene_id in self._ids("ensgene_v101", category)
        )
        # The attribution is worth printing and must not cost the pipe, so it goes beside it.
        assert category in result.stderr
        assert "ensgene_v101" in result.stderr

        named = runner.invoke(
            app, ["annotation", "gene-list", "sacCer3", category, "--annotation", "ensgene_v101"]
        )
        assert named.exit_code == 0
        assert named.stdout.splitlines() == list(self._ids("ensgene_v101", category))

        json_result = runner.invoke(app, ["annotation", "gene-list", "sacCer3", category, "--json"])
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert list(payload) == ["assembly", "annotation", "category", "gene_ids", "sources"]
        assert payload["assembly"] == "sacCer3"
        assert payload["gene_ids"] == list(self._ids("ensgene_v101", category))
        assert [list(source) for source in payload["sources"]] == [
            ["component", "annotation", "description", "source", "gene_ids"]
        ]
        assert payload["sources"][0]["component"] is None

    def test_gene_categories_prints_rows_json_and_a_merged_annotations_per_component_split(
        self, liulab_data: Path, data_dir: Path, tmp_path: Path
    ) -> None:
        self._registered(liulab_data, data_dir)
        declared = self._declared("ensgene_v101")

        result = runner.invoke(app, ["annotation", "gene-categories", "sacCer3"])
        assert result.exit_code == 0
        lines = result.stdout.splitlines()
        assert lines[0] == "categories for sacCer3 / ensgene_v101"
        assert [line.split()[0] for line in lines[1:]] == list(declared)
        assert [line.split()[1] for line in lines[1:]] == [
            str(len(self._ids("ensgene_v101", category))) for category in declared
        ]

        json_result = runner.invoke(app, ["annotation", "gene-categories", "sacCer3", "--json"])
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert [entry["category"] for entry in payload] == list(declared)
        assert all(entry["gene_ids"] for entry in payload)

        # The case #111 was opened over: worm rRNA and its food's arrive as one category
        # and must stay distinguishable inside it.
        self._merged(liulab_data, tmp_path)
        merged = runner.invoke(app, ["annotation", "gene-categories", "ce11_ecHT115"])
        assert merged.exit_code == 0
        merged_lines = merged.stdout.splitlines()
        assert (
            merged_lines[0] == "categories for ce11_ecHT115 / wormbase_ws298+refseq_rs_2025_06_26"
        )
        assert any("(ce11: " in line and "ecHT115: " in line for line in merged_lines[1:])

        shared = next(
            category
            for category in self._declared("wormbase_ws298")
            if category in self._declared("refseq_rs_2025_06_26")
        )
        gene_list = runner.invoke(
            app, ["annotation", "gene-list", "ce11_ecHT115", shared, "--json"]
        )
        assert gene_list.exit_code == 0
        merged_payload = _json.loads(gene_list.stdout)
        assert [source["component"] for source in merged_payload["sources"]] == ["ce11", "ecHT115"]
        assert merged_payload["gene_ids"] == [
            *self._ids("wormbase_ws298", shared),
            *self._ids("refseq_rs_2025_06_26", shared),
        ]

    def test_missing_curated_list_unregistered_annotation_or_unknown_category_exits_one(
        self, liulab_data: Path, data_dir: Path
    ) -> None:
        unregistered = runner.invoke(app, ["annotation", "gene-categories", "sacCer3"])
        assert unregistered.exit_code == 1
        assert "genome annotation register sacCer3 ensgene_v101" in output(unregistered)

        self._registered(liulab_data, data_dir)
        bad_category = runner.invoke(
            app, ["annotation", "gene-list", "sacCer3", "no_such_category"]
        )
        assert bad_category.exit_code == 1
        assert bad_category.stdout == ""
        assert "no_such_category" in output(bad_category)
        for category in self._declared("ensgene_v101"):
            assert category in output(bad_category)

        # Not an empty list of genes and not exit 0: the caller has to be able to tell
        # *nothing is known here* from *there are none of these genes*.
        _register("tiny", liulab_data / "genome" / "tiny", data_dir / "tiny.gtf", "mine")
        no_list = runner.invoke(app, ["annotation", "gene-list", "tiny", "rRNA"])
        assert no_list.exit_code == 1
        assert no_list.stdout == ""
        assert "no curated gene list ships" in output(no_list)
        assert "ensgene_v101" in output(no_list)  # …and which annotations do have one


class TestWhatTheseSurfacesPrintUnchanged:
    """``list`` and ``register-gtf`` print what they printed, asserted line for line.

    Both are pinned to their whole output rather than to a phrase inside it, because
    "nothing at all" is the claim: a line added to either of them fails here. The
    ``--json`` half is pinned the same way and for a stronger reason: a script parses it
    positionally as often as by key, so a reordered key is a break nobody would see.
    Only the spelling that reaches them moved; ``doctor``, whose whole output is pinned the
    same way, is asserted beside the other root commands in ``tests/test_cli.py``.
    """

    def test_annotations_json_is_the_same_keys_in_the_same_order(self) -> None:
        result = runner.invoke(app, ["annotation", "list", "hg38", "--json"])

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
            app,
            ["annotation", "register-gtf", "tiny", str(data_dir / "tiny.gtf"), "mine", "--json"],
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
        result = runner.invoke(app, ["annotation", "list", "hg38"])

        assert result.exit_code == 0
        assert result.stdout == (
            f"annotations for hg38 in {tmp_path / 'genome' / 'hg38'}\n"
            f"  gencode_v50  offered, not registered  GENCODE v50\n"
            f"default: gencode_v50 — not registered here; register it with "
            f"`genome annotation register hg38 gencode_v50`\n"
        )
