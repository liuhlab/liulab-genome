"""Tests for genome.io.gtf — registering annotations and building gffutils databases.

Databases are built **for real** from the committed ``tiny.gtf`` fixtures: the build is
most of what registering an annotation is, so substituting gffutils would test nothing.
That needs only gffutils (pure Python + SQLite), never the native bioinformatics
binaries, so nothing here is gated on a tool skip.

Nothing touches the network either: the shared ``fake_fetch`` fixture replaces the
package's one fetch step with a copy out of ``tests/data``, and the annotation table is
injected as an in-memory :class:`AnnotationMetadata` record rather than faked as a TSV —
the shipped table is tested in test_metadata.
"""

from __future__ import annotations

import gzip
import json
from dataclasses import asdict
from pathlib import Path

import gffutils
import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.gene_list import (
    GeneCategoryNotDeclaredError,
    GeneListAssemblyMismatchError,
    NoGeneCategoriesError,
    curated_gene_list,
)
from genome.io import gtf as gtf_module
from genome.io.completion import (
    RegistrationMismatchError,
    UnfinishedRegistrationError,
    read_record,
    record_path,
    work_dir,
)
from genome.io.gtf import (
    AnnotationNotRegisteredError,
    AnnotationRegistry,
    ChromosomeMismatchError,
    GtfAnnotation,
    MergeSource,
    NoCofactorTableError,
    NoGeneFeaturesError,
    NoTFCensusError,
    UnknownSpeciesError,
    _reject_unknown_chromosomes,
    annotation_dir,
    annotation_status,
    default_annotation,
    discard_merged_annotation,
    gene_list,
    gene_lists,
    list_annotations,
    list_broken_annotations,
    register_annotation,
    register_gtf,
    register_merged_gtf,
    tf_cofactor_list,
    tf_gene_list,
)
from genome.io.registration import AssemblyDir
from genome.io.results import chromosome_check_summary
from genome.io.utils import ChecksumMismatchError
from genome.metadata import AnnotationMetadata, assembly_metadata
from genome.tf.cofactor import UNIFORM_COLUMNS as COFACTOR_UNIFORM_COLUMNS
from genome.tf.cofactor import CofactorTable, cofactor_table
from genome.tf.gene import UNIFORM_COLUMNS, TFGeneTable, tf_gene_table

from .conftest import FakeFetch
from .test_source import _module_level_imports

# A minimal but valid GTF: one gene with a transcript and an exon. Standard
# gene/transcript features are declared, so the default no-inference path applies.
_GTF = (
    "\n".join(
        [
            'chrI\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
            'chrI\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
            'chrI\ttest\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
        ]
    )
    + "\n"
)

# A bare exon-level GTF: exon lines and nothing else, which is what gene/transcript
# inference exists for. Built with inference off it yields a database of exons alone.
_BARE_GTF = (
    "\n".join(
        [
            'chrI\ttest\texon\t1\t50\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
            'chrI\ttest\texon\t60\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
        ]
    )
    + "\n"
)

#: sha256 of the committed ``tiny.gtf`` — the *unpacked* bytes ``tiny.gtf.gz`` yields.
_TINY_GTF_SHA256 = "255f43bd9abef76424d1c2d89a40cccc1a36215409bbc8f32dcead49ca3baf5e"

#: A URL that is nothing like any provider's, so using it can only come from a row.
_PINNED_URL = "https://mirror.example.invalid/annotations/tiny.gtf.gz"

#: The name the fixture annotation is registered under throughout.
_NAME = "ensgene_v101"


def _row(
    *,
    name: str = _NAME,
    url: str = _PINNED_URL,
    sha256: str | None = None,
    default: bool = True,
) -> AnnotationMetadata:
    """An in-memory annotation row for the ``tiny`` assembly."""
    return AnnotationMetadata(
        assembly="tiny",
        name=name,
        provider="UCSC",
        version="ensGene.v101",
        url=url,
        sha256=sha256,
        default=default,
    )


def _register_by_name(
    assembly_dir: Path,
    assembly: str,
    name: str,
    *,
    force: bool = False,
    progressbar: bool = True,
    metadata: AnnotationMetadata | None = None,
    check_chromosomes: bool = True,
) -> GtfAnnotation:
    """Register ``name`` through a registry bound to ``assembly_dir``.

    A registry addressed by directory rather than opened, which is what these tests want:
    the assembly's ``chrom.sizes`` is found under the directory exactly as an opened
    assembly's is.
    """
    return AnnotationRegistry(AssemblyDir(assembly=assembly, path=assembly_dir)).register(
        name,
        force=force,
        progressbar=progressbar,
        metadata=metadata,
        check_chromosomes=check_chromosomes,
    )


def _register_by_path(
    assembly_dir: Path,
    gtf: str | Path,
    name: str,
    *,
    assembly: str = "tiny",
    chrom_sizes: str | Path | None = None,
    force: bool = False,
    check_chromosomes: bool = True,
    disable_infer_genes: bool = True,
    disable_infer_transcripts: bool = True,
) -> GtfAnnotation:
    """Register the GTF at ``gtf`` under ``assembly_dir`` through a registry.

    The setup half of most of this module: a registry addressed by directory rather than
    opened, answering with the paths, which is what a test asserting on files wants. Left
    to itself it checks against ``<assembly_dir>/<assembly>.chrom.sizes`` — absent in most
    of these directories, so nothing is checked — and ``chrom_sizes`` names another file
    where a test has written one somewhere else.
    """
    return AnnotationRegistry(
        AssemblyDir(assembly=assembly, path=assembly_dir), chrom_sizes=chrom_sizes
    ).register_path(
        gtf,
        name,
        force=force,
        check_chromosomes=check_chromosomes,
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )


def _feature_types(database_path: Path) -> list[str]:
    """The kinds of feature a built database holds, with the connection closed behind us."""
    database = gffutils.FeatureDB(str(database_path))
    try:
        return sorted(database.featuretypes())
    finally:
        database.conn.close()


def _write_chrom_sizes(assembly_dir: Path, *names: str, assembly: str = "tiny") -> Path:
    """Write ``<assembly>.chrom.sizes`` naming ``names``, where an assembly's own sits."""
    assembly_dir.mkdir(parents=True, exist_ok=True)
    path = assembly_dir / f"{assembly}.chrom.sizes"
    path.write_text("".join(f"{name}\t10000\n" for name in names))
    return path


class TestRegisterByPath:
    """``AnnotationRegistry.register_path`` — the way in for a GTF the table does not list.

    Placing, decompressing and rebuilding, asserted on the paths the registry answers
    with. What a broken directory tells the caller to run is asserted once, over in
    :class:`TestAnnotationRegistry`, where the assembly name that composes it lives.
    """

    def test_a_plain_gtf_is_copied_built_and_recorded(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly = tmp_path / "asm"

        annotation = _register_by_path(assembly, src, "WS298")

        assert annotation.gtf == annotation_dir(assembly, "WS298") / "WS298.gtf"
        assert annotation.gtf.read_text() == _GTF
        assert annotation.db.is_file()
        assert list(list_annotations(assembly)) == ["WS298"]

        record = read_record(annotation_dir(assembly, "WS298"))
        assert record is not None
        assert record.kind == "annotation"
        assert record.name == "WS298"
        assert sorted(record.files) == ["WS298.db", "WS298.gtf"]
        assert record.source_url == str(src)

    def test_a_gzipped_gtf_is_decompressed(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf.gz"
        with gzip.open(src, "wt") as handle:
            handle.write(_GTF)
        assembly = tmp_path / "asm"

        annotation = _register_by_path(assembly, src, "WS298")

        # Stored as a plain .gtf with decompressed contents, and the db builds.
        assert annotation.gtf.suffix == ".gtf"
        assert annotation.gtf.read_text() == _GTF
        assert annotation.db.is_file()
        assert list(list_annotations(assembly)) == ["WS298"]

    def test_a_missing_source_says_what_to_pass(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="GTF file not found"):
            _register_by_path(tmp_path / "asm", tmp_path / "nope.gtf", "X")

    def test_reregistering_a_valid_one_returns_it_without_rebuilding(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly = tmp_path / "asm"
        first = _register_by_path(assembly, src, "WS298")
        built_at = first.db.stat().st_mtime_ns

        # Silently: `filterwarnings = ["error"]` fails the test on any warning at all.
        second = _register_by_path(assembly, src, "WS298")

        assert second == first
        assert second.db.stat().st_mtime_ns == built_at

    def test_force_rebuilds(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly = tmp_path / "asm"
        _register_by_path(assembly, src, "WS298")

        annotation = _register_by_path(assembly, src, "WS298", force=True)

        assert annotation.db.is_file()
        assert list(list_annotations(assembly)) == ["WS298"]


class TestRegisterByName:
    """``AnnotationRegistry.register`` — naming an annotation is enough to have it on disk."""

    @pytest.fixture(autouse=True)
    def _serve_the_gtf(self, fake_fetch: FakeFetch) -> FakeFetch:
        fake_fetch.serve("tiny.gtf.gz")
        return fake_fetch

    def test_it_fetches_verifies_builds_and_records(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        annotation = _register_by_name(
            tmp_path,
            "tiny",
            _NAME,
            progressbar=False,
            metadata=_row(sha256=_TINY_GTF_SHA256),
        )

        assert fake_fetch.last.url == _PINNED_URL
        assert annotation.gtf == annotation_dir(tmp_path, _NAME) / f"{_NAME}.gtf"
        assert annotation.gtf.is_file()
        assert annotation.db.is_file()

        record = read_record(annotation_dir(tmp_path, _NAME))
        assert record is not None
        assert record.kind == "annotation"
        assert record.name == _NAME
        assert record.source_url == _PINNED_URL
        assert record.sha256 == _TINY_GTF_SHA256
        assert sorted(record.files) == [f"{_NAME}.db", f"{_NAME}.gtf"]
        # The archive went with the working area once the record was written.
        assert not work_dir(annotation_dir(tmp_path, _NAME)).exists()

    def test_the_database_it_builds_answers_queries(self, tmp_path: Path) -> None:
        annotation = _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=_row())

        database = gffutils.FeatureDB(str(annotation.db))
        try:
            transcripts = list(database.features_of_type("transcript"))
            assert len(transcripts) == 18
            assert {feature.seqid for feature in transcripts} == {"chrI", "chrII", "chrIII"}
        finally:
            database.conn.close()

    def test_a_wrong_checksum_raises_naming_both_digests(self, tmp_path: Path) -> None:
        wrong = "0" * 64

        with pytest.raises(ChecksumMismatchError) as excinfo:
            _register_by_name(
                tmp_path, "tiny", _NAME, progressbar=False, metadata=_row(sha256=wrong)
            )

        assert wrong in str(excinfo.value)
        assert _TINY_GTF_SHA256 in str(excinfo.value)
        # Nothing that could not be vouched for reached the annotation's own files.
        directory = annotation_dir(tmp_path, _NAME)
        assert not (directory / f"{_NAME}.gtf").exists()
        assert read_record(directory) is None

    def test_a_row_that_pins_no_digest_records_whatever_arrived(self, tmp_path: Path) -> None:
        _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=_row())

        record = read_record(annotation_dir(tmp_path, _NAME))
        assert record is not None
        assert record.sha256 == _TINY_GTF_SHA256

    def test_an_uncompressed_url_is_placed_as_it_arrives(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        fake_fetch.serve("tiny.gtf")
        url = "https://mirror.example.invalid/annotations/tiny.gtf"

        annotation = _register_by_name(
            tmp_path, "tiny", _NAME, progressbar=False, metadata=_row(url=url)
        )

        assert annotation.gtf.read_text().startswith("chrII\tensGene.v101\ttranscript")

    def test_a_name_no_row_lists_says_what_is_offered(self, tmp_path: Path) -> None:
        # Against the shipped table, which lists exactly one annotation for sacCer3.
        with pytest.raises(ValueError, match="no annotation named 'nope'") as excinfo:
            _register_by_name(tmp_path, "sacCer3", "nope", progressbar=False)

        assert "ensgene_v101" in str(excinfo.value)
        assert "register-gtf" in str(excinfo.value)  # the way in for one no row lists
        assert "register_path" in str(excinfo.value)  # ...and the same from Python

    def test_reregistering_a_valid_one_is_a_silent_no_op(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        row = _row(sha256=_TINY_GTF_SHA256)
        first = _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)
        built_at = first.db.stat().st_mtime_ns

        # Silently: `filterwarnings = ["error"]` fails the test on any warning at all.
        second = _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)

        assert second == first
        assert second.db.stat().st_mtime_ns == built_at  # not rebuilt
        assert len(fake_fetch.calls) == 1  # nothing fetched twice

    def test_a_half_built_annotation_is_reported_as_broken(self, tmp_path: Path) -> None:
        # A gffutils build killed part-way: a database file, and no record.
        directory = annotation_dir(tmp_path, _NAME)
        directory.mkdir(parents=True)
        (directory / f"{_NAME}.db").write_bytes(b"half a database")

        with pytest.raises(UnfinishedRegistrationError) as excinfo:
            _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=_row())

        assert f"genome register-annotation tiny {_NAME} --force" in str(excinfo.value)

    def test_a_record_that_disagrees_with_disk_raises(self, tmp_path: Path) -> None:
        row = _row(sha256=_TINY_GTF_SHA256)
        annotation = _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)
        annotation.db.write_bytes(b"truncated")

        with pytest.raises(RegistrationMismatchError, match="disagrees with its"):
            _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)

    def test_force_repairs_what_the_error_named(self, tmp_path: Path) -> None:
        directory = annotation_dir(tmp_path, _NAME)
        directory.mkdir(parents=True)
        (directory / f"{_NAME}.db").write_bytes(b"half a database")

        annotation = _register_by_name(
            tmp_path, "tiny", _NAME, progressbar=False, force=True, metadata=_row()
        )

        assert read_record(directory) is not None
        assert annotation.db.stat().st_size > len(b"half a database")

    def test_force_keeps_a_gtf_whose_digest_still_matches(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        row = _row(sha256=_TINY_GTF_SHA256)
        _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)
        record_path(annotation_dir(tmp_path, _NAME)).unlink()

        _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, force=True, metadata=row)

        assert len(fake_fetch.calls) == 1  # the GTF on disk proved itself; nothing refetched

    def test_force_refetches_when_the_row_pins_nothing_to_prove_it_against(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        row = _row()
        _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)

        _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, force=True, metadata=row)

        assert len(fake_fetch.calls) == 2

    def test_the_record_carries_gffutils_rather_than_a_tool_version(self, tmp_path: Path) -> None:
        # gffutils is a Python library, not an External tool resolved on PATH, so its
        # version is provenance in details and never a tool version.
        _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=_row())

        record = read_record(annotation_dir(tmp_path, _NAME))
        assert record is not None
        assert record.tool_versions == {}
        assert record.details["gffutils_version"] == gffutils.__version__
        assert record.details["provider"] == "UCSC"
        assert record.details["version"] == "ensGene.v101"


class TestListAnnotations:
    """Registered means a record that agrees with disk — never a database file."""

    def test_a_database_without_a_record_is_not_registered(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        _register_by_path(tmp_path, src, "finished")
        halfway = annotation_dir(tmp_path, "halfway")
        halfway.mkdir(parents=True)
        (halfway / "halfway.db").write_bytes(b"half a database")

        assert list(list_annotations(tmp_path)) == ["finished"]

    def test_a_record_that_disagrees_with_disk_is_not_registered(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        annotation = _register_by_path(tmp_path, src, "WS298")
        assert list(list_annotations(tmp_path)) == ["WS298"]

        annotation.db.write_bytes(b"truncated")

        assert list_annotations(tmp_path) == {}


class TestDiscardMergedAnnotation:
    """What a chimera build removes when its next build merges under another name."""

    def test_a_merged_annotation_goes_and_takes_an_emptied_gtf_tree_with_it(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly = tmp_path / "asm"
        chrom_sizes = _write_chrom_sizes(assembly, "chrI__tiny")
        register_merged_gtf(
            assembly,
            "a+b",
            [MergeSource("tiny", "a", src)],
            separator="__",
            chrom_sizes=chrom_sizes,
        )

        assert discard_merged_annotation(assembly, "a+b") is True

        assert list_annotations(assembly) == {}
        # A chimera that merges nothing carries no gtf/ tree, and one whose last derived
        # annotation has just gone is in exactly that state.
        assert not (assembly / "gtf").exists()

    def test_an_annotation_a_caller_registered_by_hand_is_never_removed(
        self, tmp_path: Path
    ) -> None:
        # The name comes from a previous build's record, and a name is not ownership: only
        # a record showing a merge wrote it is.
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        _register_by_path(tmp_path, src, "a+b")

        assert discard_merged_annotation(tmp_path, "a+b") is False
        assert list(list_annotations(tmp_path)) == ["a+b"]

    def test_a_name_nothing_is_registered_under_is_not_an_error(self, tmp_path: Path) -> None:
        assert discard_merged_annotation(tmp_path, "a+b") is False


class TestListBrokenAnnotations:
    """The complement of ``list_annotations``: what is on disk and cannot be trusted.

    Every directory under ``gtf/`` is registered, broken, or not begun; this reports
    the middle one, which is otherwise invisible to anything that lists.
    """

    def test_a_directory_holding_files_with_no_record_is_broken(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        _register_by_path(tmp_path, src, "mine")
        record_path(annotation_dir(tmp_path, "mine")).unlink()

        broken = list_broken_annotations(tmp_path, "tiny")

        assert list(broken) == ["mine"]
        assert broken["mine"].directory == annotation_dir(tmp_path, "mine")
        assert "holds files but no .completion.json" in broken["mine"].problem

    def test_a_record_that_disagrees_with_disk_is_broken(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        annotation = _register_by_path(tmp_path, src, "mine")
        annotation.db.write_bytes(b"truncated")

        broken = list_broken_annotations(tmp_path, "tiny")

        assert list(broken) == ["mine"]
        assert "mine.db" in broken["mine"].problem

    def test_a_finished_annotation_is_not_broken(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        _register_by_path(tmp_path, src, "mine")

        assert list_broken_annotations(tmp_path, "tiny") == {}

    def test_an_empty_directory_is_not_broken(self, tmp_path: Path) -> None:
        # ADR-0007: an absent or empty directory is a fresh registration, not a broken
        # one. A run interrupted before it downloaded anything must not be reported.
        annotation_dir(tmp_path, "mine").mkdir(parents=True)
        work_dir(annotation_dir(tmp_path, "mine")).mkdir()

        assert list_broken_annotations(tmp_path, "tiny") == {}

    def test_an_assembly_with_no_annotations_at_all_answers_empty(self, tmp_path: Path) -> None:
        assert list_broken_annotations(tmp_path, "sacCer3") == {}

    def test_a_broken_one_leaves_the_others_listed(self, tmp_path: Path) -> None:
        # The invariant: one broken annotation must not stop the rest being found, and
        # nothing here raises — reporting is the whole point.
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        _register_by_path(tmp_path, src, "healthy")
        _register_by_path(tmp_path, src, "damaged")
        record_path(annotation_dir(tmp_path, "damaged")).unlink()

        assert list(list_annotations(tmp_path)) == ["healthy"]
        assert list(list_broken_annotations(tmp_path, "tiny")) == ["damaged"]

    def test_a_name_the_table_offers_is_repaired_by_name(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        _register_by_path(tmp_path, src, "ensgene_v101")
        record_path(annotation_dir(tmp_path, "ensgene_v101")).unlink()

        broken = list_broken_annotations(tmp_path, "sacCer3")

        assert broken["ensgene_v101"].repair == (
            "genome register-annotation sacCer3 ensgene_v101 --force"
        )
        assert broken["ensgene_v101"].repair in broken["ensgene_v101"].problem

    def test_an_unlisted_one_is_repaired_from_the_path_its_record_remembers(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        annotation = _register_by_path(tmp_path, src, "mine")
        annotation.db.write_bytes(b"truncated")

        broken = list_broken_annotations(tmp_path, "tiny")

        assert broken["mine"].repair == f"genome register-gtf tiny {src} mine --force"

    def test_an_unlisted_one_whose_source_is_unknowable_says_so(self, tmp_path: Path) -> None:
        # No record survives to say which GTF it was built from, so there is no path to
        # print: the command is named with the one thing it still needs filled in,
        # rather than a path that would not run.
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        _register_by_path(tmp_path, src, "mine")
        record_path(annotation_dir(tmp_path, "mine")).unlink()

        broken = list_broken_annotations(tmp_path, "tiny")

        assert broken["mine"].repair == "genome register-gtf tiny <path> mine --force"

    def test_an_unlisted_one_whose_source_is_gone_is_not_named_as_a_command(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        annotation = _register_by_path(tmp_path, src, "mine")
        annotation.db.write_bytes(b"truncated")
        src.unlink()

        broken = list_broken_annotations(tmp_path, "tiny")

        assert str(src) not in broken["mine"].repair
        assert broken["mine"].repair == "genome register-gtf tiny <path> mine --force"


class TestDefaultAnnotation:
    """The one rule that decides a default, wherever the question is asked from."""

    def test_the_flagged_row_decides_it(self) -> None:
        assert default_annotation([_row()], []) == _NAME

    def test_the_flag_wins_over_the_sole_registered_annotation(self) -> None:
        # Everyone in the lab reaches for the same one, whatever this machine happens
        # to hold.
        assert default_annotation([_row()], ["something_else"]) == _NAME

    def test_an_explicit_choice_wins_over_the_flag(self) -> None:
        assert default_annotation([_row()], [_NAME], explicit="something_else") == "something_else"

    def test_nothing_flagged_falls_back_to_the_sole_registered_annotation(self) -> None:
        assert default_annotation([_row(default=False)], ["only_one"]) == "only_one"

    def test_nothing_flagged_and_no_sole_annotation_leaves_no_default(self) -> None:
        assert default_annotation([], []) is None
        assert default_annotation([], ["one", "two"]) is None


class TestAnnotationStatus:
    """What an assembly's table offers, set against what is registered on this machine."""

    def test_it_reports_what_is_offered_with_nothing_registered(self, liulab_data: Path) -> None:
        # The case it most needs to serve: a fresh machine, where the answer is
        # entirely the shipped table's.
        payload = annotation_status("sacCer3")

        assert payload.assembly == "sacCer3"
        assert payload.directory == liulab_data / "genome" / "sacCer3"
        assert payload.default_annotation == "ensgene_v101"
        rows = payload.annotations
        assert [(r.name, r.offered, r.registered) for r in rows] == [("ensgene_v101", True, False)]
        assert rows[0].provider == "UCSC"
        assert rows[0].path is None

    def test_the_payload_it_serializes_is_the_rows_under_their_own_names(
        self, liulab_data: Path
    ) -> None:
        # `--json` is this report rendered, so a row's fields and the payload's keys are
        # one spelling: a surface reads attributes and never names a key of its own.
        payload = annotation_status("sacCer3")

        assert payload.as_json() == {
            "assembly": "sacCer3",
            "directory": str(liulab_data / "genome" / "sacCer3"),
            "default_annotation": "ensgene_v101",
            "annotations": [asdict(row) for row in payload.annotations],
        }

    def test_the_default_annotations_own_row_is_reachable_without_a_search(self) -> None:
        # What the closing line of `genome annotations` needs: the default's own state,
        # so "not registered here" and "broken here" are told apart by the report itself.
        offered = annotation_status("sacCer3")
        nothing = annotation_status("tiny")

        default = offered.default_row
        assert default is offered.annotations[0]
        assert default is not None
        assert not default.registered
        # No default decided, so no row is about one — the state a fresh unlisted
        # assembly is in, and not an error.
        assert nothing.default_annotation is None
        assert nothing.default_row is None

    def test_it_creates_nothing_and_fetches_nothing(
        self, fake_fetch: FakeFetch, liulab_data: Path
    ) -> None:
        annotation_status("hg38")

        assert fake_fetch.calls == []
        assert not (liulab_data / "genome" / "hg38").exists()

    def test_a_registered_annotation_the_table_offers_is_reported_as_both(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly_dir = tmp_path / "asm"
        annotation = _register_by_path(assembly_dir, src, "ensgene_v101")

        payload = annotation_status("sacCer3", cache_dir=assembly_dir)

        rows = payload.annotations
        assert [(r.name, r.offered, r.registered) for r in rows] == [("ensgene_v101", True, True)]
        assert rows[0].path == str(annotation.gtf)

    def test_a_registered_annotation_no_row_lists_is_reported_too(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly_dir = tmp_path / "asm"
        _register_by_path(assembly_dir, src, "mine")

        payload = annotation_status("tiny", cache_dir=assembly_dir)

        rows = payload.annotations
        assert [(r.name, r.offered, r.registered) for r in rows] == [("mine", False, True)]
        assert rows[0].provider is None
        assert rows[0].broken is False
        assert rows[0].problem is None
        assert payload.default_annotation == "mine"  # nothing flagged, and it is alone

    def test_a_broken_offered_annotation_is_reported_as_broken_not_as_absent(
        self, tmp_path: Path
    ) -> None:
        # The bug this closes: half-registered and never-fetched looked identical here.
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly_dir = tmp_path / "asm"
        _register_by_path(assembly_dir, src, "ensgene_v101")
        record_path(annotation_dir(assembly_dir, "ensgene_v101")).unlink()

        payload = annotation_status("sacCer3", cache_dir=assembly_dir)

        rows = payload.annotations
        assert [(r.name, r.offered, r.registered, r.broken) for r in rows] == [
            ("ensgene_v101", True, False, True)
        ]
        assert rows[0].repair == "genome register-annotation sacCer3 ensgene_v101 --force"
        assert "holds files but no .completion.json" in str(rows[0].problem)
        assert rows[0].path is None

    def test_a_broken_unlisted_annotation_is_reported_at_all(self, tmp_path: Path) -> None:
        # No row lists it and no record vouches for it, so nothing used to mention it.
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly_dir = tmp_path / "asm"
        annotation = _register_by_path(assembly_dir, src, "mine")
        annotation.db.write_bytes(b"truncated")

        payload = annotation_status("tiny", cache_dir=assembly_dir)

        rows = payload.annotations
        assert [(r.name, r.offered, r.registered, r.broken) for r in rows] == [
            ("mine", False, False, True)
        ]
        assert rows[0].repair == f"genome register-gtf tiny {src} mine --force"

    def test_one_broken_annotation_does_not_hide_the_others(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly_dir = tmp_path / "asm"
        _register_by_path(assembly_dir, src, "healthy")
        _register_by_path(assembly_dir, src, "damaged")
        record_path(annotation_dir(assembly_dir, "damaged")).unlink()

        payload = annotation_status("tiny", cache_dir=assembly_dir)

        rows = payload.annotations
        assert [(r.name, r.registered, r.broken) for r in rows] == [
            ("damaged", False, True),
            ("healthy", True, False),
        ]
        # A broken one is not a registered one, so it never becomes the sole-registered
        # default either.
        assert payload.default_annotation == "healthy"


class TestAnnotationRegistry:
    """One assembly's annotations, and the four states each of them can be in.

    The registry is the one place *registered / broken / offered / not begun* is
    assembled: everything that used to rebuild that four-way state — a genome opening,
    the status report, the error a name nobody registered earns — asks this instead.
    """

    def _registered(self, tmp_path: Path, name: str) -> Path:
        """Register the fixture GTF under ``name`` and return the source it came from."""
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(tmp_path, source, name)
        return source

    def test_the_four_states_are_settled_in_one_construction(self, tmp_path: Path) -> None:
        self._registered(tmp_path, "ensgene_v101")
        self._registered(tmp_path, "damaged")
        record_path(annotation_dir(tmp_path, "damaged")).unlink()

        registry = AnnotationRegistry.locate("sacCer3", tmp_path)

        assert registry.registered == ["ensgene_v101"]
        assert [entry.name for entry in registry.broken] == ["damaged"]
        assert [record.name for record in registry.offered] == ["ensgene_v101"]
        assert registry.default == "ensgene_v101"

    def test_it_reads_the_directory_it_was_pointed_at(self, tmp_path: Path) -> None:
        # The assembly dir travels with the registry rather than being re-derived from
        # the data root at each question.
        elsewhere = tmp_path / "elsewhere"
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(elsewhere, source, "mine")

        assert AnnotationRegistry.locate("tiny", elsewhere).registered == ["mine"]
        assert AnnotationRegistry.locate("tiny", tmp_path).registered == []

    def test_the_path_of_a_registered_annotation_is_its_gtf(self, tmp_path: Path) -> None:
        self._registered(tmp_path, "mine")

        registry = AnnotationRegistry.locate("tiny", tmp_path)

        assert registry.path("mine") == annotation_dir(tmp_path, "mine") / "mine.gtf"

    def test_a_name_nothing_knows_says_what_is_registered_and_what_is_offered(
        self, tmp_path: Path
    ) -> None:
        registry = AnnotationRegistry.locate("sacCer3", tmp_path)

        with pytest.raises(AnnotationNotRegisteredError) as excinfo:
            registry.path("no_such_annotation")

        message = str(excinfo.value)
        assert "no_such_annotation" in message
        assert "ensgene_v101" in message  # what the table does offer

    def test_the_path_of_a_broken_one_names_its_repair_in_one_hop(self, tmp_path: Path) -> None:
        # Not the command that would itself raise and demand --force: the one that works.
        self._registered(tmp_path, "ensgene_v101")
        record_path(annotation_dir(tmp_path, "ensgene_v101")).unlink()

        registry = AnnotationRegistry.locate("sacCer3", tmp_path)

        with pytest.raises(AnnotationNotRegisteredError) as excinfo:
            registry.path("ensgene_v101")

        assert "genome register-annotation sacCer3 ensgene_v101 --force" in str(excinfo.value)

    def test_registering_by_path_adopts_it_without_reading_the_disk_again(
        self, tmp_path: Path
    ) -> None:
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        registry = AnnotationRegistry.locate("tiny", tmp_path)
        assert registry.registered == []

        annotation = registry.register_path(source, "mine")

        assert registry.registered == ["mine"]
        assert registry.path("mine") == annotation.gtf
        assert registry.default == "mine"  # nothing flagged, and it is alone

    def test_registering_by_name_fetches_and_adopts(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        fake_fetch.serve("tiny.gtf.gz")
        registry = AnnotationRegistry.locate("tiny", tmp_path)

        annotation = registry.register(_NAME, progressbar=False, metadata=_row())

        assert registry.registered == [_NAME]
        assert registry.path(_NAME) == annotation.gtf

    def test_a_default_already_decided_is_never_displaced_by_a_registration(
        self, tmp_path: Path
    ) -> None:
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        registry = AnnotationRegistry.locate("sacCer3", tmp_path)
        assert registry.default == "ensgene_v101"  # the table's flag

        registry.register_path(source, "mine")

        assert registry.default == "ensgene_v101"

    def test_an_explicit_default_wins_and_need_not_be_registered(self, tmp_path: Path) -> None:
        registry = AnnotationRegistry.locate("sacCer3", tmp_path, default="mine")

        assert registry.default == "mine"
        assert registry.registered == []

    def test_registering_over_a_broken_directory_stops_reporting_it(self, tmp_path: Path) -> None:
        source = self._registered(tmp_path, "mine")
        record_path(annotation_dir(tmp_path, "mine")).unlink()
        registry = AnnotationRegistry.locate("tiny", tmp_path)
        assert [entry.name for entry in registry.broken] == ["mine"]

        registry.register_path(source, "mine", force=True)

        assert registry.broken == []
        assert registry.registered == ["mine"]

    def test_addressed_by_assembly_name_a_broken_directory_names_a_command(
        self, tmp_path: Path
    ) -> None:
        # A registry always knows its assembly, so the repair it names is a command a
        # shell can run rather than a Python call with the assembly left to guess at.
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        directory = annotation_dir(tmp_path, "WS298")
        directory.mkdir(parents=True)
        (directory / "WS298.db").write_bytes(b"half a database")

        with pytest.raises(UnfinishedRegistrationError) as excinfo:
            AnnotationRegistry.locate("tiny", tmp_path).register_path(source, "WS298")

        assert f"genome register-gtf tiny {source} WS298 --force" in str(excinfo.value)

    def test_it_finds_the_assemblys_chrom_sizes_without_being_told(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        _write_chrom_sizes(tmp_path, "chrI", "chrII", "chrIII")

        with pytest.raises(ChromosomeMismatchError):
            AnnotationRegistry.locate("tiny", tmp_path).register_path(
                data_dir / "ensembl_style.gtf", _NAME
            )

    def test_a_chrom_sizes_it_is_handed_is_the_one_checked_against(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        # What a Genome opened somewhere of its own passes: the chrom.sizes it actually
        # prepared, rather than the one the layout would name for its assembly.
        sizes = _write_chrom_sizes(tmp_path, "chrI", "chrII", "chrIII", assembly="elsewhere")

        with pytest.raises(ChromosomeMismatchError):
            AnnotationRegistry(
                AssemblyDir.locate("tiny", tmp_path), chrom_sizes=sizes
            ).register_path(data_dir / "ensembl_style.gtf", _NAME)

    def test_its_status_is_what_the_status_report_answers(self, tmp_path: Path) -> None:
        self._registered(tmp_path, "healthy")
        self._registered(tmp_path, "damaged")
        record_path(annotation_dir(tmp_path, "damaged")).unlink()

        assert AnnotationRegistry.locate("tiny", tmp_path).status() == annotation_status(
            "tiny", cache_dir=tmp_path
        )

    def test_nothing_is_created_by_asking(self, liulab_data: Path) -> None:
        registry = AnnotationRegistry.locate("sacCer3", liulab_data / "genome" / "sacCer3")

        assert registry.registered == []
        assert registry.broken == []
        assert registry.default == "ensgene_v101"
        assert not (liulab_data / "genome").exists()


class TestChromosomeNames:
    """A GTF's sequence names must be the assembly's, and are checked before the build.

    The mismatch case is the committed ``ensembl_style.gtf`` — ``tiny.gtf``'s own 85
    features with the ``chr`` prefix stripped — against a ``chrI``/``chrII``/``chrIII``
    assembly: the UCSC-versus-Ensembl case in real bytes.
    """

    #: How the fixture assembly spells its three sequences.
    _UCSC = ("chrI", "chrII", "chrIII")

    def test_an_ensembl_spelled_gtf_is_refused_naming_every_offender(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        sizes = _write_chrom_sizes(tmp_path, *self._UCSC)

        with pytest.raises(ChromosomeMismatchError) as excinfo:
            _register_by_path(tmp_path, data_dir / "ensembl_style.gtf", _NAME, chrom_sizes=sizes)

        assert excinfo.value.missing == ("I", "II", "III")
        message = str(excinfo.value)
        assert "I, II, III" in message
        assert "chrI" in message  # what the assembly spells them as
        assert "check_chromosomes=False" in message

    def test_a_refused_gtf_costs_neither_a_database_nor_a_directory(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        # The check runs before the build, and before anything is placed: the annotation
        # directory is left exactly as it was found.
        sizes = _write_chrom_sizes(tmp_path, *self._UCSC)

        with pytest.raises(ChromosomeMismatchError):
            _register_by_path(tmp_path, data_dir / "ensembl_style.gtf", _NAME, chrom_sizes=sizes)

        assert not annotation_dir(tmp_path, _NAME).exists()
        assert list(tmp_path.rglob("*.db")) == []

    def test_refused_by_name_it_stays_refused_rather_than_reading_as_interrupted(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        fake_fetch.serve("ensembl_style.gtf")
        _write_chrom_sizes(tmp_path, *self._UCSC)
        row = _row(url="https://mirror.example.invalid/annotations/ensembl_style.gtf")

        with pytest.raises(ChromosomeMismatchError):
            _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)

        directory = annotation_dir(tmp_path, _NAME)
        assert not (directory / f"{_NAME}.gtf").exists()  # never placed
        assert list(tmp_path.rglob("*.db")) == []  # never paid for the build
        assert read_record(directory) is None

        # Running it again reports the same problem, not an interrupted registration.
        with pytest.raises(ChromosomeMismatchError):
            _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)

    def test_sequences_the_annotation_never_mentions_are_not_an_error(self, tmp_path: Path) -> None:
        # Strict one way only: the GTF names chrI alone, the assembly carries five.
        sizes = _write_chrom_sizes(tmp_path, *self._UCSC, "chrM", "scaffold_17")
        source = tmp_path / "one-chromosome.gtf"
        source.write_text(_GTF)

        annotation = _register_by_path(tmp_path, source, "WS298", chrom_sizes=sizes)

        assert annotation.db.is_file()
        record = read_record(annotation_dir(tmp_path, "WS298"))
        assert record is not None
        assert record.details["chromosomes_checked"] is True

    def test_a_gzipped_source_is_checked_without_being_unpacked_first(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        source = tmp_path / "ensembl_style.gtf.gz"
        with gzip.open(source, "wt") as handle:
            handle.write((data_dir / "ensembl_style.gtf").read_text())
        sizes = _write_chrom_sizes(tmp_path, *self._UCSC)

        with pytest.raises(ChromosomeMismatchError):
            _register_by_path(tmp_path, source, _NAME, chrom_sizes=sizes)

        assert not annotation_dir(tmp_path, _NAME).exists()

    def test_a_wholesale_mismatch_lists_ten_names_and_counts_the_rest(self, tmp_path: Path) -> None:
        offenders = [f"scaffold_{n}" for n in range(25)]
        source = tmp_path / "many.gtf"
        source.write_text(
            "".join(f'{name}\ttest\texon\t1\t100\t.\t+\t.\tgene_id "g1";\n' for name in offenders)
        )
        sizes = _write_chrom_sizes(tmp_path, *self._UCSC)

        with pytest.raises(ChromosomeMismatchError) as excinfo:
            _register_by_path(tmp_path, source, _NAME, chrom_sizes=sizes)

        assert len(excinfo.value.missing) == 25  # every one of them is on the exception
        assert "(and 15 more)" in str(excinfo.value)  # ten of them are in the message

    def test_comment_lines_are_not_taken_for_chromosomes(self, tmp_path: Path) -> None:
        sizes = _write_chrom_sizes(tmp_path, *self._UCSC)
        source = tmp_path / "commented.gtf"
        source.write_text("##description: a header\n#!genome-build tiny\n" + _GTF)

        assert _register_by_path(tmp_path, source, "WS298", chrom_sizes=sizes).db.is_file()

    def test_the_override_registers_a_mismatched_gtf_anyway(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        sizes = _write_chrom_sizes(tmp_path, *self._UCSC)

        annotation = _register_by_path(
            tmp_path,
            data_dir / "ensembl_style.gtf",
            _NAME,
            chrom_sizes=sizes,
            check_chromosomes=False,
        )

        database = gffutils.FeatureDB(str(annotation.db))
        try:
            assert {feature.seqid for feature in database.features_of_type("transcript")} == {
                "I",
                "II",
                "III",
            }
        finally:
            database.conn.close()
        record = read_record(annotation_dir(tmp_path, _NAME))
        assert record is not None
        assert record.details["chromosomes_checked"] is False
        assert record.details["chromosomes_unchecked_because"] == "caller-override"

    def test_the_override_registers_by_name_too(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        fake_fetch.serve("ensembl_style.gtf")
        _write_chrom_sizes(tmp_path, *self._UCSC)
        row = _row(url="https://mirror.example.invalid/annotations/ensembl_style.gtf")

        annotation = _register_by_name(
            tmp_path, "tiny", _NAME, progressbar=False, metadata=row, check_chromosomes=False
        )

        assert annotation.db.is_file()
        record = read_record(annotation_dir(tmp_path, _NAME))
        assert record is not None
        assert record.details["chromosomes_checked"] is False
        assert record.details["chromosomes_unchecked_because"] == "caller-override"

    def test_a_matching_gtf_registered_by_name_records_that_it_was_checked(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        fake_fetch.serve("tiny.gtf.gz")
        _write_chrom_sizes(tmp_path, *self._UCSC)

        _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=_row())

        record = read_record(annotation_dir(tmp_path, _NAME))
        assert record is not None
        assert record.details["chromosomes_checked"] is True
        # A check that ran and did not raise passed, so there is no reason beside it.
        assert record.details["chromosomes_unchecked_because"] is None

    def test_without_a_chrom_sizes_there_is_nothing_to_check_and_the_record_says_so(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        # An annotation registered before its assembly was prepared: no chrom.sizes
        # exists, so the names cannot be checked. The record says they were not, and
        # says it was for want of that file rather than because anyone asked to skip it.
        fake_fetch.serve("ensembl_style.gtf")
        row = _row(url="https://mirror.example.invalid/annotations/ensembl_style.gtf")

        annotation = _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)

        assert annotation.db.is_file()
        record = read_record(annotation_dir(tmp_path, _NAME))
        assert record is not None
        assert record.details["chromosomes_checked"] is False
        assert record.details["chromosomes_unchecked_because"] == "no-chrom-sizes"


@pytest.fixture(scope="session")
def gtf_scratch(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A reusable GTF path for the property test (session-scoped: hypothesis-safe)."""
    return tmp_path_factory.mktemp("chromosomes") / "generated.gtf"


#: Sequence names that look like the ones references actually use, kept free of
#: whitespace so every generated line is a well-formed GTF record.
_chrom_name = st.text(alphabet="chrIVXM_.0123456789", min_size=1, max_size=12)


@given(
    in_gtf=st.lists(_chrom_name, max_size=8),
    in_assembly=st.lists(_chrom_name, max_size=8),
)
def test_the_check_names_every_offender_never_a_subset(
    gtf_scratch: Path, in_gtf: list[str], in_assembly: list[str]
) -> None:
    gtf_scratch.write_text(
        "".join(f'{name}\ttest\texon\t1\t100\t.\t+\t.\tgene_id "g1";\n' for name in in_gtf)
    )
    known = frozenset(in_assembly)
    missing = set(in_gtf) - known

    if not missing:
        _reject_unknown_chromosomes(gtf_scratch, known, name=_NAME)
        return
    with pytest.raises(ChromosomeMismatchError) as excinfo:
        _reject_unknown_chromosomes(gtf_scratch, known, name=_NAME)
    assert excinfo.value.missing == tuple(sorted(missing))


class TestRegisterAnnotation:
    """``register_annotation`` — the same by assembly name, answering with the record."""

    def test_it_returns_the_record_plus_where_it_landed(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        fake_fetch.serve("tiny.gtf.gz")

        payload = register_annotation(
            "tiny",
            _NAME,
            cache_dir=tmp_path,
            progressbar=False,
            metadata=_row(sha256=_TINY_GTF_SHA256),
        )

        assert payload.record.kind == "annotation"
        assert payload.name == _NAME
        assert payload.assembly == "tiny"
        assert payload.directory == annotation_dir(tmp_path, _NAME)
        assert payload.source_url == _PINNED_URL
        assert payload.sha256 == _TINY_GTF_SHA256
        directory = annotation_dir(tmp_path, _NAME)
        assert payload.record.files == {
            name: (directory / name).stat().st_size for name in (f"{_NAME}.gtf", f"{_NAME}.db")
        }
        assert payload.file_names == [f"{_NAME}.db", f"{_NAME}.gtf"]

    def test_the_payload_it_serializes_is_the_record_plus_where_it_landed(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        # The `--json` payload is the completion record under its own on-disk key names,
        # with the two facts a record does not hold about itself. A type wraps those
        # names; it never renames them, because lab directories are read by both.
        fake_fetch.serve("tiny.gtf.gz")

        payload = register_annotation(
            "tiny", _NAME, cache_dir=tmp_path, progressbar=False, metadata=_row()
        )

        assert payload.as_json() == {
            **asdict(payload.record),
            "assembly": "tiny",
            "directory": str(annotation_dir(tmp_path, _NAME)),
        }
        assert list(payload.as_json())[-2:] == ["assembly", "directory"]

    def test_what_the_chromosome_check_settled_is_read_off_the_record(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        # The sentence belongs to the record and to the API that reads it, so a surface
        # printing it names none of the record's fields. Nothing is registered as the
        # assembly here, so there was no chrom.sizes to check against.
        fake_fetch.serve("tiny.gtf.gz")

        payload = register_annotation(
            "tiny", _NAME, cache_dir=tmp_path, progressbar=False, metadata=_row()
        )

        assert payload.chromosome_check == chromosome_check_summary(payload.record.details)
        assert "nothing to check against" in payload.chromosome_check

    def test_the_chromosome_check_reaches_this_way_in_too(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        fake_fetch.serve("ensembl_style.gtf")
        _write_chrom_sizes(tmp_path, "chrI", "chrII", "chrIII")
        row = _row(url="https://mirror.example.invalid/annotations/ensembl_style.gtf")

        with pytest.raises(ChromosomeMismatchError):
            register_annotation("tiny", _NAME, cache_dir=tmp_path, progressbar=False, metadata=row)

        payload = register_annotation(
            "tiny",
            _NAME,
            cache_dir=tmp_path,
            progressbar=False,
            metadata=row,
            check_chromosomes=False,
        )

        assert payload.record.details == {
            "provider": "UCSC",
            "version": "ensGene.v101",
            "gffutils_version": gffutils.__version__,
            "chromosomes_checked": False,
            "chromosomes_unchecked_because": "caller-override",
        }

    def test_it_files_the_annotation_under_the_assembly_data_dir(
        self, fake_fetch: FakeFetch, liulab_data: Path
    ) -> None:
        fake_fetch.serve("tiny.gtf.gz")

        payload = register_annotation("tiny", _NAME, progressbar=False, metadata=_row())

        assert payload.directory == liulab_data / "genome" / "tiny" / "gtf" / _NAME

    def test_the_inference_knobs_reach_the_database_build(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        # A bare exon-level GTF: with inference left off the database holds exons and
        # nothing else, so a caller who asks for genes and transcripts must be able to
        # get them from this way in too.
        bare = tmp_path / "bare.gtf"
        bare.write_text(_BARE_GTF)
        fake_fetch.serve(bare)
        row = _row(url="https://mirror.example.invalid/annotations/bare.gtf")

        default = register_annotation(
            "tiny", "exons_only", cache_dir=tmp_path, progressbar=False, metadata=row
        )
        inferred = register_annotation(
            "tiny",
            "with_genes",
            cache_dir=tmp_path,
            progressbar=False,
            metadata=row,
            disable_infer_genes=False,
            disable_infer_transcripts=False,
        )

        assert default.name == "exons_only"
        assert _feature_types(annotation_dir(tmp_path, "exons_only") / "exons_only.db") == ["exon"]
        assert inferred.name == "with_genes"
        assert _feature_types(annotation_dir(tmp_path, "with_genes") / "with_genes.db") == [
            "exon",
            "gene",
            "transcript",
        ]


class TestRegisterGtf:
    """``register_gtf`` — a GTF no row lists, addressed by assembly name."""

    def test_it_returns_the_record_plus_where_it_landed(
        self, tmp_path: Path, liulab_data: Path
    ) -> None:
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)

        payload = register_gtf("tiny", source, "WS298")

        directory = liulab_data / "genome" / "tiny" / "gtf" / "WS298"
        assert payload.record.kind == "annotation"
        assert payload.name == "WS298"
        assert payload.assembly == "tiny"
        assert payload.directory == directory
        assert payload.source_url == str(source)
        assert payload.record.files == {
            name: (directory / name).stat().st_size for name in ("WS298.gtf", "WS298.db")
        }

    def test_cache_dir_overrides_which_assembly_directory_it_is_filed_under(
        self, tmp_path: Path
    ) -> None:
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        elsewhere = tmp_path / "elsewhere"

        payload = register_gtf("tiny", source, "WS298", cache_dir=elsewhere)

        assert payload.directory == annotation_dir(elsewhere, "WS298")

    def test_it_finds_the_assembly_chrom_sizes_without_being_told(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        # Naming the assembly is what says where its chrom.sizes is, so an
        # Ensembl-spelled GTF is refused rather than silently registered unchecked.
        _write_chrom_sizes(tmp_path, "chrI", "chrII", "chrIII")

        with pytest.raises(ChromosomeMismatchError):
            register_gtf("tiny", data_dir / "ensembl_style.gtf", _NAME, cache_dir=tmp_path)

    def test_the_override_registers_the_mismatch_anyway_and_the_record_says_so(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        _write_chrom_sizes(tmp_path, "chrI", "chrII", "chrIII")

        payload = register_gtf(
            "tiny",
            data_dir / "ensembl_style.gtf",
            _NAME,
            cache_dir=tmp_path,
            check_chromosomes=False,
        )

        assert payload.record.details == {
            "gffutils_version": gffutils.__version__,
            "chromosomes_checked": False,
            "chromosomes_unchecked_because": "caller-override",
        }

    def test_the_inference_knobs_reach_the_database_build(self, tmp_path: Path) -> None:
        source = tmp_path / "bare.gtf"
        source.write_text(_BARE_GTF)

        register_gtf("tiny", source, "exons_only", cache_dir=tmp_path)
        register_gtf(
            "tiny",
            source,
            "with_genes",
            cache_dir=tmp_path,
            disable_infer_genes=False,
            disable_infer_transcripts=False,
        )

        assert _feature_types(annotation_dir(tmp_path, "exons_only") / "exons_only.db") == ["exon"]
        assert _feature_types(annotation_dir(tmp_path, "with_genes") / "with_genes.db") == [
            "exon",
            "gene",
            "transcript",
        ]

    def test_a_missing_source_says_what_to_pass(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="GTF file not found"):
            register_gtf("tiny", tmp_path / "nope.gtf", "WS298", cache_dir=tmp_path)

    def test_a_broken_directory_names_the_command_that_repairs_it(self, tmp_path: Path) -> None:
        # Addressed by assembly name, the repair is a command a shell can run rather
        # than a Python call with the assembly left to guess at.
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        directory = annotation_dir(tmp_path, "WS298")
        directory.mkdir(parents=True)
        (directory / "WS298.db").write_bytes(b"half a database")

        with pytest.raises(UnfinishedRegistrationError) as excinfo:
            register_gtf("tiny", source, "WS298", cache_dir=tmp_path)

        assert f"genome register-gtf tiny {source} WS298 --force" in str(excinfo.value)

    def test_force_repairs_what_the_error_named(self, tmp_path: Path) -> None:
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        directory = annotation_dir(tmp_path, "WS298")
        directory.mkdir(parents=True)
        (directory / "WS298.db").write_bytes(b"half a database")

        payload = register_gtf("tiny", source, "WS298", cache_dir=tmp_path, force=True)

        assert payload.name == "WS298"
        assert list(list_annotations(tmp_path)) == ["WS298"]


# ---------------------------------------------------------------------------------------
# Which genes are in a category, and what an annotation that cannot say answers
# ---------------------------------------------------------------------------------------

#: An annotation whose curated gene list ships today, and the assembly the table files it
#: under. Named rather than derived: registering it under *another* assembly is what the
#: guard below is about, so the pairing has to be written down somewhere the test controls.
_CURATED = "gencode_v50"
_CURATED_ASSEMBLY = "hg38"

#: The two contributors a **Merged annotation** of worm and its food is made of, and the
#: components whose sequences their features sit on.
_WORM, _WORM_COMPONENT = "wormbase_ws298", "ce11"
_FOOD, _FOOD_COMPONENT = "refseq_rs_2025_06_26", "ecHT115"
_CHIMERA = f"{_WORM_COMPONENT}_{_FOOD_COMPONENT}"


def _declared(annotation: str) -> tuple[str, ...]:
    """The categories a shipped curated list declares, in file order — never a fixed set.

    Which categories exist is the curated list's to say and differs per annotation, so
    every assertion below reads them off the shipped file rather than naming them.
    """
    listed = curated_gene_list(annotation)
    assert listed is not None, f"no curated gene list ships for {annotation}"
    return tuple(listed.categories)


def _curated_ids(annotation: str, category: str) -> tuple[str, ...]:
    """The gene ids a shipped curated list puts in ``category``."""
    listed = curated_gene_list(annotation)
    assert listed is not None
    return listed.categories[category].gene_ids


def _register_merged(
    assembly_dir: Path,
    name: str,
    source: Path,
    *,
    assembly: str = _CHIMERA,
    food: str = _FOOD,
) -> GtfAnnotation:
    """Register a merged annotation of the worm/food pair under ``name``, from one GTF.

    A real :func:`register_merged_gtf`, not a hand-written record: what a merge writes
    into ``details`` is what the gene-category path reads back, so standing that in would
    test the stand-in. The fixture GTF is merged twice under two component names, which is
    all the record needs to carry — nothing here reads the features.

    ``food`` names the second contributor's annotation, so a caller may merge in one no
    curated list ships for and test what a contributor that cannot contribute does.
    """
    chrom_sizes = _write_chrom_sizes(
        assembly_dir,
        *(f"chrI__{component}" for component in (_WORM_COMPONENT, _FOOD_COMPONENT)),
        assembly=assembly,
    )
    return register_merged_gtf(
        assembly_dir,
        name,
        [
            MergeSource(_WORM_COMPONENT, _WORM, source),
            MergeSource(_FOOD_COMPONENT, food, source),
        ],
        separator="__",
        chrom_sizes=chrom_sizes,
    )


class TestGeneList:
    """``AnnotationRegistry.gene_list`` — the genes one annotation puts in one category.

    The shipped curated lists answer here, since which categories exist is data and no
    fixture may pretend otherwise. What is asserted is structure: who contributed, what
    order the ids arrive in, and — the point of the whole surface — that an annotation
    which cannot answer raises rather than answering with nothing.
    """

    def _registry(self, tmp_path: Path, *, assembly: str = _CURATED_ASSEMBLY) -> AnnotationRegistry:
        """Register the fixture GTF under a curated annotation's name and open the registry."""
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(tmp_path, source, _CURATED, assembly=assembly)
        return AnnotationRegistry.locate(assembly, tmp_path)

    def test_a_plain_annotation_answers_as_one_source_belonging_to_no_component(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path)
        category = _declared(_CURATED)[0]

        answer = registry.gene_list(category, _CURATED)

        assert (answer.assembly, answer.annotation, answer.category) == (
            _CURATED_ASSEMBLY,
            _CURATED,
            category,
        )
        # One contributor, and `component` is None rather than the assembly's own name:
        # attribution is what a merge needs, and there is nothing here to attribute.
        assert [source.component for source in answer.sources] == [None]
        assert [source.annotation for source in answer.sources] == [_CURATED]

    def test_the_ids_are_the_curated_lists_own_in_the_order_it_lists_them(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path)
        category = _declared(_CURATED)[0]

        answer = registry.gene_list(category, _CURATED)

        assert answer.gene_ids == list(_curated_ids(_CURATED, category))
        assert answer.gene_ids  # never an empty answer: a declared category has genes

    def test_a_source_carries_what_membership_means_and_where_it_came_from(
        self, tmp_path: Path
    ) -> None:
        # The two sentences travel with the ids, because they are what says whether these
        # ids mean what the caller's metric needs.
        answer = self._registry(tmp_path).gene_list(_declared(_CURATED)[0], _CURATED)

        assert answer.sources[0].description.strip()
        assert answer.sources[0].source.strip()

    def test_a_category_the_list_does_not_declare_raises_and_lists_the_ones_it_does(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path)

        with pytest.raises(GeneCategoryNotDeclaredError) as excinfo:
            registry.gene_list("no_such_category", _CURATED)

        message = str(excinfo.value)
        assert "no_such_category" in message
        for category in _declared(_CURATED):
            assert category in message

    def test_an_annotation_nothing_ships_a_list_for_raises_the_other_absence(
        self, tmp_path: Path
    ) -> None:
        # The fact #111 exists for: *no categories are declared* and *this category is not
        # declared* are different, and neither is an empty answer.
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(tmp_path, source, "mine")
        registry = AnnotationRegistry.locate("tiny", tmp_path)

        with pytest.raises(NoGeneCategoriesError) as excinfo:
            registry.gene_list("rRNA", "mine")

        message = str(excinfo.value)
        assert "mine" in message
        assert _CURATED in message  # …and which annotations do declare categories

    def test_the_two_absences_are_told_apart_by_type_and_caught_together(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path)
        source = tmp_path / "ann.gtf"
        _register_by_path(tmp_path, source, "mine", assembly=_CURATED_ASSEMBLY)
        both = AnnotationRegistry.locate(_CURATED_ASSEMBLY, tmp_path)

        with pytest.raises(LookupError) as declared:
            registry.gene_list("no_such_category", _CURATED)
        with pytest.raises(LookupError) as nothing:
            both.gene_list("no_such_category", "mine")

        assert isinstance(declared.value, GeneCategoryNotDeclaredError)
        assert isinstance(nothing.value, NoGeneCategoriesError)
        assert not isinstance(declared.value, NoGeneCategoriesError)
        assert not isinstance(nothing.value, GeneCategoryNotDeclaredError)

    def test_an_unregistered_name_earns_the_error_it_already_had(self, tmp_path: Path) -> None:
        # Callers pass names, never paths, and a name nothing registered is resolved by the
        # same `path` every other question goes through — so the message is the same one.
        registry = AnnotationRegistry.locate(_CURATED_ASSEMBLY, tmp_path)

        with pytest.raises(AnnotationNotRegisteredError) as excinfo:
            registry.gene_list("rRNA", _CURATED)

        assert f"genome register-annotation {_CURATED_ASSEMBLY} {_CURATED}" in str(excinfo.value)

    def test_naming_no_annotation_asks_the_default_one(self, tmp_path: Path) -> None:
        registry = self._registry(tmp_path)
        category = _declared(_CURATED)[0]

        assert registry.default == _CURATED
        assert registry.gene_list(category).annotation == _CURATED

    def test_no_default_and_no_name_says_which_argument_chooses_one(self, tmp_path: Path) -> None:
        registry = AnnotationRegistry.locate("tiny", tmp_path)

        with pytest.raises(ValueError, match="annotation") as excinfo:
            registry.gene_list("rRNA")

        assert "default_gtf" in str(excinfo.value)

    def test_a_list_curated_against_another_assembly_refuses_to_answer(
        self, tmp_path: Path
    ) -> None:
        # A name is unique only within its assembly, so a list found by name alone is not
        # yet known to be about this reference — and answering would hand back another
        # species' genes under this one's name.
        registry = self._registry(tmp_path, assembly="tiny")

        with pytest.raises(GeneListAssemblyMismatchError) as excinfo:
            registry.gene_list(_declared(_CURATED)[0], _CURATED)

        message = str(excinfo.value)
        assert "tiny" in message
        assert _CURATED_ASSEMBLY in message

    def test_gene_lists_returns_every_declared_category_in_file_order(self, tmp_path: Path) -> None:
        answers = self._registry(tmp_path).gene_lists(_CURATED)

        assert [answer.category for answer in answers] == list(_declared(_CURATED))
        assert all(answer.gene_ids for answer in answers)

    def test_gene_lists_raises_rather_than_answering_with_an_empty_tuple(
        self, tmp_path: Path
    ) -> None:
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(tmp_path, source, "mine")

        with pytest.raises(NoGeneCategoriesError):
            AnnotationRegistry.locate("tiny", tmp_path).gene_lists("mine")


class TestMergedGeneList:
    """A **Merged annotation**'s genes, attributed to the component each came from.

    The case #111 was opened over: a chimera of a worm and the bacterium it eats, whose
    ribosomal RNA is both species' and must not arrive as one number.
    """

    def _registry(self, tmp_path: Path, name: str = f"{_WORM}+{_FOOD}") -> AnnotationRegistry:
        """Register a merged annotation under ``name`` and open the chimera's registry."""
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_merged(tmp_path, name, source)
        return AnnotationRegistry.locate(_CHIMERA, tmp_path)

    def test_each_contributor_answers_for_its_own_component(self, tmp_path: Path) -> None:
        registry = self._registry(tmp_path)
        shared = next(category for category in _declared(_WORM) if category in _declared(_FOOD))

        answer = registry.gene_list(shared, f"{_WORM}+{_FOOD}")

        assert [(source.component, source.annotation) for source in answer.sources] == [
            (_WORM_COMPONENT, _WORM),
            (_FOOD_COMPONENT, _FOOD),
        ]
        assert answer.annotation == f"{_WORM}+{_FOOD}"
        assert answer.assembly == _CHIMERA

    def test_the_ids_are_the_contributions_concatenated_and_never_de_duplicated(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path)
        shared = next(category for category in _declared(_WORM) if category in _declared(_FOOD))

        answer = registry.gene_list(shared, f"{_WORM}+{_FOOD}")

        assert answer.gene_ids == [
            *_curated_ids(_WORM, shared),
            *_curated_ids(_FOOD, shared),
        ]

    def test_a_contributor_that_cannot_answer_is_left_out_rather_than_raised_over(
        self, tmp_path: Path
    ) -> None:
        # One component's annotation is curated and the other's is not, so only one of the
        # two can contribute. That is an omission and not a failure: the answer carries the
        # contributor that can answer, and stays attributed to its component.
        gtf = tmp_path / "ann.gtf"
        gtf.write_text(_GTF)
        _register_merged(tmp_path, "merged", gtf, food="nobody_curated_this")
        registry = AnnotationRegistry.locate(_CHIMERA, tmp_path)
        category = _declared(_WORM)[0]

        answer = registry.gene_list(category, "merged")

        assert [contributed.component for contributed in answer.sources] == [_WORM_COMPONENT]
        assert answer.gene_ids == list(_curated_ids(_WORM, category))

    def test_a_category_no_contributor_declares_raises(self, tmp_path: Path) -> None:
        registry = self._registry(tmp_path)

        with pytest.raises(GeneCategoryNotDeclaredError):
            registry.gene_list("no_such_category", f"{_WORM}+{_FOOD}")

    def test_who_contributed_comes_from_the_record_and_not_from_splitting_the_name(
        self, tmp_path: Path
    ) -> None:
        # The name of a merge is the +-join of what went in, but it cannot say which
        # component each half came from — and that is exactly what attribution needs.
        registry = self._registry(tmp_path, name="merged")
        shared = next(category for category in _declared(_WORM) if category in _declared(_FOOD))

        answer = registry.gene_list(shared, "merged")

        assert [source.component for source in answer.sources] == [
            _WORM_COMPONENT,
            _FOOD_COMPONENT,
        ]

    def test_gene_lists_is_the_union_of_what_the_contributors_declare(self, tmp_path: Path) -> None:
        registry = self._registry(tmp_path)
        union = list(dict.fromkeys([*_declared(_WORM), *_declared(_FOOD)]))

        assert [answer.category for answer in registry.gene_lists()] == union


# ---------------------------------------------------------------------------------------
# Which of an annotation's own gene ids a Gene id stem names
# ---------------------------------------------------------------------------------------

#: A pseudoautosomal gene as GENCODE spells it on a lift: one stem, two gene ids, the
#: second the Y copy. Real ids — PLCXD1 and TP53 — because the shape of the collision is
#: the whole point and an invented id would not have it.
_PAR_X = "ENSG00000182378.14"
_PAR_Y = "ENSG00000182378.14_PAR_Y"
_ALONE = "ENSG00000141510.18"
_PAR_STEM = "ENSG00000182378"
_ALONE_STEM = "ENSG00000141510"

#: A stem no fixture below carries a gene for — what the census holds and the annotation
#: does not.
_ABSENT_STEM = "ENSG00000288541"


def _gtf_declaring(*gene_ids: str) -> str:
    """A GTF declaring one gene, transcript and exon per id, in the shape GENCODE has.

    The ids are the only thing that varies: nothing reading it looks at a coordinate, so
    every gene sits on the same interval rather than inviting anyone to.
    """
    lines: list[str] = []
    for gene_id in gene_ids:
        attributes = f'gene_id "{gene_id}"; transcript_id "{gene_id}_t";'
        lines.extend(
            f"chrI\ttest\t{feature}\t1\t100\t.\t+\t.\t{attributes}"
            for feature in ("gene", "transcript", "exon")
        )
    return "\n".join(lines) + "\n"


def _registry_declaring(
    tmp_path: Path, *gene_ids: str, assembly: str, name: str = "mine"
) -> AnnotationRegistry:
    """Register a GTF declaring ``gene_ids`` under ``name`` and open the registry.

    The setup both shipped-table surfaces share — a census against an annotation and a
    cofactor table against one differ in which file is read and in nothing about the
    annotation, so they register the same way. ``assembly`` is required: which species an
    assembly is is exactly what those tests are about, and a default would hide it.
    """
    tmp_path.mkdir(parents=True, exist_ok=True)
    source = tmp_path / f"{name}.gtf"
    source.write_text(_gtf_declaring(*gene_ids))
    _register_by_path(tmp_path, source, name, assembly=assembly)
    return AnnotationRegistry.locate(assembly, tmp_path)


class TestResolveGeneIds:
    """``AnnotationRegistry.resolve_gene_ids`` — a **Gene id stem** against real gene ids.

    The first surface here to open the **Annotation database**, so every test registers a
    fixture GTF for real and asks the registry, exactly as the gene-category tests do. What
    is asserted is what a caller holding a published table keyed by unversioned ids gets
    back: both of the ids a stem names rather than one of them, and the stems this
    annotation has nothing for said out loud.
    """

    def _registry(self, tmp_path: Path, gtf: str, name: str = "mine") -> AnnotationRegistry:
        """Register ``gtf`` under ``name`` for the ``tiny`` assembly and open the registry."""
        source = tmp_path / f"{name}.gtf"
        source.write_text(gtf)
        _register_by_path(tmp_path, source, name)
        return AnnotationRegistry.locate("tiny", tmp_path)

    def test_a_stem_resolves_to_the_versioned_id_the_annotation_spells_it_with(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path, _gtf_declaring(_ALONE))

        answer = registry.resolve_gene_ids([_ALONE_STEM], "mine")

        assert answer.resolved == {_ALONE_STEM: (_ALONE,)}
        assert answer.unresolved == ()
        assert (answer.assembly, answer.annotation) == ("tiny", "mine")

    def test_a_stem_naming_two_gene_ids_answers_with_both_rather_than_choosing(
        self, tmp_path: Path
    ) -> None:
        # The `gencode_v50lift37` case: nine stems name two genes each, eight of them a
        # pseudoautosomal pair. Answering with the first would hand back the X copy of a Y
        # gene without ever saying a choice had been made.
        registry = self._registry(tmp_path, _gtf_declaring(_PAR_Y, _PAR_X, _ALONE))

        answer = registry.resolve_gene_ids([_PAR_STEM], "mine")

        assert answer.resolved[_PAR_STEM] == (_PAR_X, _PAR_Y)

    def test_a_stem_this_annotation_carries_no_gene_for_is_reported_and_not_dropped(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path, _gtf_declaring(_ALONE))

        answer = registry.resolve_gene_ids([_ALONE_STEM, _ABSENT_STEM], "mine")

        assert answer.unresolved == (_ABSENT_STEM,)
        assert _ABSENT_STEM not in answer.resolved
        # …and what did resolve is still there, so an absence costs the caller nothing else.
        assert answer.gene_ids == [_ALONE]

    def test_an_annotation_whose_ids_carry_no_version_resolves_each_stem_to_itself(
        self, tmp_path: Path
    ) -> None:
        # WormBase and SGD never versioned a gene id, and an Ensembl-shaped assumption
        # would leave both unresolvable. An id with no version is its own stem.
        registry = self._registry(tmp_path, _GTF)

        answer = registry.resolve_gene_ids(["g1", "g2"], "mine")

        assert answer.resolved == {"g1": ("g1",)}
        assert answer.unresolved == ("g2",)

    def test_the_answer_keeps_the_order_the_stems_were_asked_about_and_asks_repeats_once(
        self, tmp_path: Path
    ) -> None:
        # A caller passing a few thousand at once reads its own list against the answer.
        registry = self._registry(tmp_path, _gtf_declaring(_ALONE, _PAR_X))

        answer = registry.resolve_gene_ids(
            [_PAR_STEM, _ABSENT_STEM, _ALONE_STEM, _PAR_STEM], "mine"
        )

        assert list(answer.resolved) == [_PAR_STEM, _ALONE_STEM]
        assert answer.gene_ids == [_PAR_X, _ALONE]
        assert answer.unresolved == (_ABSENT_STEM,)

    def test_asking_about_no_stems_answers_emptily_rather_than_about_every_gene(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path, _gtf_declaring(_ALONE))

        answer = registry.resolve_gene_ids([], "mine")

        assert (dict(answer.resolved), answer.unresolved) == ({}, ())

    def test_an_annotation_declaring_no_genes_says_so_rather_than_failing_every_stem(
        self, tmp_path: Path
    ) -> None:
        # An exon-level GTF registers as exons alone, since reconstructing the features
        # above them is off by default. Every stem would come back unresolved, which reads
        # as *this annotation has none of your genes* and is a different fact entirely.
        registry = self._registry(tmp_path, _BARE_GTF)

        with pytest.raises(NoGeneFeaturesError) as excinfo:
            registry.resolve_gene_ids(["g1"], "mine")

        message = str(excinfo.value)
        assert "mine" in message
        assert "--infer-genes" in message  # …and the argument that rebuilds it with genes

    def test_an_unregistered_name_earns_the_error_it_already_had(self, tmp_path: Path) -> None:
        # Resolved through the same lookup every other question goes through, so the
        # message and its repair are the ones that surface already has.
        registry = AnnotationRegistry.locate(_CURATED_ASSEMBLY, tmp_path)

        with pytest.raises(AnnotationNotRegisteredError) as excinfo:
            registry.resolve_gene_ids([_ALONE_STEM], _CURATED)

        assert f"genome register-annotation {_CURATED_ASSEMBLY} {_CURATED}" in str(excinfo.value)

    def test_naming_no_annotation_asks_the_default_one(self, tmp_path: Path) -> None:
        registry = self._registry(tmp_path, _gtf_declaring(_ALONE))

        assert registry.default == "mine"
        assert registry.resolve_gene_ids([_ALONE_STEM]).annotation == "mine"

    def test_a_merged_annotation_answers_from_its_own_database_for_either_component(
        self, tmp_path: Path
    ) -> None:
        # A merge rewrites seqnames and never a `gene_id`, so its one database holds both
        # components' genes under the ids their own annotations gave them — and resolution
        # reads that database, exactly as it reads any other annotation's.
        worm = tmp_path / "worm.gtf"
        worm.write_text(_gtf_declaring("WBGene00004512"))
        food = tmp_path / "food.gtf"
        food.write_text(_gtf_declaring(_ALONE))
        chrom_sizes = _write_chrom_sizes(
            tmp_path,
            *(f"chrI__{component}" for component in (_WORM_COMPONENT, _FOOD_COMPONENT)),
            assembly=_CHIMERA,
        )
        register_merged_gtf(
            tmp_path,
            "merged",
            [
                MergeSource(_WORM_COMPONENT, _WORM, worm),
                MergeSource(_FOOD_COMPONENT, _FOOD, food),
            ],
            separator="__",
            chrom_sizes=chrom_sizes,
        )
        registry = AnnotationRegistry.locate(_CHIMERA, tmp_path)

        answer = registry.resolve_gene_ids(["WBGene00004512", _ALONE_STEM], "merged")

        assert answer.gene_ids == ["WBGene00004512", _ALONE]
        assert answer.unresolved == ()


# ---------------------------------------------------------------------------------------
# Which of an annotation's genes a published census judges transcription factors
# ---------------------------------------------------------------------------------------

#: An assembly whose species has a census, one whose species has none, and one nothing
#: names a species for at all — the chimera, which is more than one species by
#: construction. Named rather than derived: which species a census ships for is exactly
#: what these tests are about, so the pairing is written down where the test controls it.
_CENSUSED_ASSEMBLY = _CURATED_ASSEMBLY
_UNCENSUSED_ASSEMBLY = "ce11"

#: ZBED1 as GENCODE spells it on a pseudoautosomal region: one **Gene id stem**, two gene
#: ids, and a gene Lambert judges a transcription factor — so the collision a resolver must
#: not silently pick from is a **TF gene**'s own, rather than an invented one.
_TF_PAR_X = "ENSG00000214717.13"
_TF_PAR_Y = "ENSG00000214717.13_PAR_Y"
_TF_PAR_STEM = "ENSG00000214717"


def _census() -> TFGeneTable:
    """The census shipped for the species of the assembly the tests below register against.

    Read off the shipped file rather than named, exactly as the curated gene lists are:
    which genes are transcription factors is the census's judgement, and no fixture may
    stand in for it.
    """
    species = assembly_metadata(_CENSUSED_ASSEMBLY).species
    assert species is not None, f"{_CENSUSED_ASSEMBLY} has no species in the assembly table"
    census = tf_gene_table(species)
    assert census is not None, f"no census ships for {species}"
    return census


def _census_row(stem: str) -> dict[str, str | None]:
    """The census's own row for one **Gene id stem**, keyed by its own column names."""
    census = _census()
    row = census.rows[census.gene_id_stems.index(stem)]
    return dict(zip(census.columns, row, strict=True))


def _rejected_stem() -> str:
    """A stem the census assessed and judged *not* to be a transcription factor."""
    census = _census()
    positive = set(census.assessed_positive)
    return next(stem for stem in census.gene_id_stems if stem not in positive)


def _stem_assessed(assessment: str) -> str:
    """A stem the census judged a TF under one of its own **TF assessment** grades."""
    return next(
        stem
        for stem in _census().assessed_positive
        if _census_row(stem)["tf_assessment"] == assessment
    )


def _versioned(stem: str) -> str:
    """One gene id an annotation might spell that stem with — a version it never carries."""
    return f"{stem}.1"


class TestTFGeneList:
    """``AnnotationRegistry.tf_gene_list`` — where a published census meets an annotation.

    The shipped census answers here, as the shipped curated lists do for gene categories:
    what a census holds is its publisher's business and no fixture pretends otherwise. What
    is asserted is the crossing — the census's **Gene id stem**s arriving as this
    annotation's own gene ids, what rode back unresolved, and that an assembly no census
    can answer for raises rather than answering with nothing.
    """

    def _registry(
        self,
        tmp_path: Path,
        *gene_ids: str,
        assembly: str = _CENSUSED_ASSEMBLY,
        name: str = "mine",
    ) -> AnnotationRegistry:
        """Register a GTF declaring ``gene_ids`` under ``name`` and open the registry."""
        return _registry_declaring(tmp_path, *gene_ids, assembly=assembly, name=name)

    def test_the_census_arrives_in_this_annotations_own_gene_ids(self, tmp_path: Path) -> None:
        # The whole point of the surface: the answer joins to a counts matrix keyed by this
        # annotation's ids, with no normalisation left for the caller.
        stems = _census().assessed_positive[:2]
        gene_ids = [_versioned(stem) for stem in stems]
        registry = self._registry(tmp_path, *gene_ids)

        answer = registry.tf_gene_list("mine")

        assert (answer.assembly, answer.annotation) == (_CENSUSED_ASSEMBLY, "mine")
        assert [gene.gene_id_stem for gene in answer.genes] == list(stems)
        assert answer.gene_ids == gene_ids

    def test_only_the_genes_the_census_judged_transcription_factors_are_carried(
        self, tmp_path: Path
    ) -> None:
        # The common case is not 2,765 rows to filter down to 1,639.
        positive, rejected = _census().assessed_positive[0], _rejected_stem()
        registry = self._registry(tmp_path, _versioned(positive), _versioned(rejected))

        answer = registry.tf_gene_list("mine")

        assert [gene.gene_id_stem for gene in answer.genes] == [positive]
        assert [gene.is_tf for gene in answer.genes] == [True]

    def test_a_caller_can_widen_to_the_genes_the_census_assessed_and_rejected(
        self, tmp_path: Path
    ) -> None:
        # …and widening carries the verdict rather than dropping it: a gene the census
        # assessed and turned down arrives saying so, which is not the same fact as a gene
        # it never looked at, and that one is absent from both answers.
        positive, rejected = _census().assessed_positive[0], _rejected_stem()
        registry = self._registry(tmp_path, _versioned(positive), _versioned(rejected))

        answer = registry.tf_gene_list("mine", include_rejected=True)

        assert {gene.gene_id_stem: gene.is_tf for gene in answer.genes} == {
            positive: True,
            rejected: False,
        }

    def test_a_gene_carries_the_censuss_family_and_every_judgement_it_recorded(
        self, tmp_path: Path
    ) -> None:
        census = _census()
        stem = census.assessed_positive[0]
        registry = self._registry(tmp_path, _versioned(stem))

        gene = registry.tf_gene_list("mine").genes[0]

        cells = _census_row(stem)
        assert gene.symbol == cells["symbol"]
        assert gene.dbd_family == cells["dbd_family"]
        # Everything the publisher recorded beyond the uniform four, under its own names:
        # the assessment, the binding mode, the motif status, the KRAB flag and the votes.
        assert dict(gene.judgements) == {
            name: cells[name] for name in census.columns[len(UNIFORM_COLUMNS) :]
        }

    def test_a_caller_can_tighten_on_the_assessment_the_census_recorded(
        self, tmp_path: Path
    ) -> None:
        # The **TF assessment** is graded, and tightening to `Known motif` or loosening to
        # include `Inferred motif` is a re-filter on what the answer already carries rather
        # than a second flag this package invents.
        known, inferred = _stem_assessed("Known motif"), _stem_assessed("Inferred motif")
        registry = self._registry(tmp_path, _versioned(known), _versioned(inferred))

        answer = registry.tf_gene_list("mine")

        assert {gene.gene_id_stem for gene in answer.genes} == {known, inferred}
        assert [
            gene.gene_id_stem
            for gene in answer.genes
            if gene.judgements["tf_assessment"] == "Known motif"
        ] == [known]

    def test_the_verdict_travels_with_the_census_that_reached_it(self, tmp_path: Path) -> None:
        registry = self._registry(tmp_path, _versioned(_census().assessed_positive[0]))

        answer = registry.tf_gene_list("mine")

        assert answer.provenance == _census().provenance
        assert answer.provenance.publisher
        assert answer.provenance.version
        assert answer.provenance.pubmed_id
        assert answer.species == assembly_metadata(_CENSUSED_ASSEMBLY).species

    def test_stems_this_annotation_carries_no_gene_for_ride_back_on_the_answer(
        self, tmp_path: Path
    ) -> None:
        census = _census()
        carried, absent = census.assessed_positive[0], census.assessed_positive[1]
        registry = self._registry(tmp_path, _versioned(carried))

        answer = registry.tf_gene_list("mine")

        assert [gene.gene_id_stem for gene in answer.genes] == [carried]
        assert absent in answer.unresolved
        # Every stem the census judged a transcription factor is accounted for one way or
        # it holds and this annotation does not is visible rather than dropped.
        assert len(answer.genes) + len(answer.unresolved) == len(census.assessed_positive)

    def test_a_stem_naming_two_gene_ids_answers_with_both(self, tmp_path: Path) -> None:
        assert _TF_PAR_STEM in _census().assessed_positive
        registry = self._registry(tmp_path, _TF_PAR_Y, _TF_PAR_X)

        answer = registry.tf_gene_list("mine")

        assert [gene.gene_ids for gene in answer.genes] == [(_TF_PAR_X, _TF_PAR_Y)]
        assert answer.gene_ids == [_TF_PAR_X, _TF_PAR_Y]

    def test_the_species_follows_the_assembly_and_never_the_ids_the_gtf_holds(
        self, tmp_path: Path
    ) -> None:
        # Human gene ids registered for a worm assembly. Asking for one species'
        # transcription factors while holding another's assembly is not expressible, so this
        # is answered about the assembly's own species and never about what is in the GTF.
        registry = self._registry(tmp_path, _TF_PAR_X, assembly=_UNCENSUSED_ASSEMBLY)

        with pytest.raises(NoTFCensusError) as excinfo:
            registry.tf_gene_list("mine")

        message = str(excinfo.value)
        assert str(assembly_metadata(_UNCENSUSED_ASSEMBLY).species) in message
        assert str(assembly_metadata(_CENSUSED_ASSEMBLY).species) in message

    @pytest.mark.parametrize("assembly", [_CHIMERA, "tiny"])
    def test_an_assembly_nothing_names_a_species_for_says_so_rather_than_guessing(
        self, tmp_path: Path, assembly: str
    ) -> None:
        registry = self._registry(tmp_path, _TF_PAR_X, assembly=assembly)

        with pytest.raises(UnknownSpeciesError) as excinfo:
            registry.tf_gene_list("mine")

        message = str(excinfo.value)
        assert assembly in message
        assert str(assembly_metadata(_CENSUSED_ASSEMBLY).species) in message

    def test_the_two_absences_are_told_apart_by_type_and_caught_together(
        self, tmp_path: Path
    ) -> None:
        # As the curated gene lists' pair already are: *no census ships for this species*
        # and *nothing says what species this is* are different answers, both lookups, and
        # neither is an empty collection.
        uncensused = self._registry(tmp_path / "worm", _TF_PAR_X, assembly=_UNCENSUSED_ASSEMBLY)
        unnamed = self._registry(tmp_path / "chimera", _TF_PAR_X, assembly=_CHIMERA)

        with pytest.raises(LookupError) as no_census:
            uncensused.tf_gene_list("mine")
        with pytest.raises(LookupError) as no_species:
            unnamed.tf_gene_list("mine")

        assert isinstance(no_census.value, NoTFCensusError)
        assert isinstance(no_species.value, UnknownSpeciesError)
        assert not isinstance(no_census.value, UnknownSpeciesError)
        assert not isinstance(no_species.value, NoTFCensusError)

    def test_an_unregistered_name_earns_the_error_it_already_had(self, tmp_path: Path) -> None:
        registry = AnnotationRegistry.locate(_CENSUSED_ASSEMBLY, tmp_path)

        with pytest.raises(AnnotationNotRegisteredError) as excinfo:
            registry.tf_gene_list(_CURATED)

        assert f"genome register-annotation {_CENSUSED_ASSEMBLY} {_CURATED}" in str(excinfo.value)

    def test_naming_no_annotation_asks_the_default_one(self, tmp_path: Path) -> None:
        registry = self._registry(tmp_path, _TF_PAR_X, name=_CURATED)

        assert registry.default == _CURATED
        assert registry.tf_gene_list().annotation == _CURATED

    def test_the_json_record_carries_the_genes_the_provenance_and_the_unresolved_stems(
        self, tmp_path: Path
    ) -> None:
        # What ``--json`` has to be able to emit: the genes with their **TF assessment** and
        # **DBD family**, the census's provenance, and the stems that resolved to nothing.
        stem = _census().assessed_positive[0]
        registry = self._registry(tmp_path, _versioned(stem))

        payload = registry.tf_gene_list("mine").as_json()

        assert payload["assembly"] == _CENSUSED_ASSEMBLY
        assert payload["gene_ids"] == [_versioned(stem)]
        assert payload["provenance"]["pubmed_id"] == _census().provenance.pubmed_id
        assert payload["unresolved"]
        gene = payload["genes"][0]
        assert (gene["gene_id_stem"], gene["gene_ids"]) == (stem, [_versioned(stem)])
        assert gene["dbd_family"] == _census_row(stem)["dbd_family"]
        assert gene["judgements"]["tf_assessment"] == _census_row(stem)["tf_assessment"]
        assert json.loads(json.dumps(payload)) == payload  # serializes as it stands

    def test_it_answers_for_an_assembly_named_rather_than_opened(self, tmp_path: Path) -> None:
        stem = _census().assessed_positive[0]
        self._registry(tmp_path, _versioned(stem))

        answer = tf_gene_list(_CENSUSED_ASSEMBLY, annotation="mine", cache_dir=tmp_path)

        assert [gene.gene_id_stem for gene in answer.genes] == [stem]
        assert answer.provenance == _census().provenance


# ---------------------------------------------------------------------------------------
# Which of an annotation's genes a publisher lists as transcription cofactors
# ---------------------------------------------------------------------------------------

#: An assembly whose species has a cofactor table, one whose species has none, and the
#: worm — which has a table although no publisher has released a worm TF census, so it is
#: the one assembly the two halves answer differently for. Named rather than derived, for
#: the reason the census pairing above is: which species ships what is what these tests
#: are about. Yeast rather than human is the untabled one on purpose: human has no table
#: only until one is built for it, and a test pinned to that would go green by accident.
_TABLED_ASSEMBLY = "mm39"
_UNTABLED_ASSEMBLY = "sacCer3"
_WORM_ASSEMBLY = "ce11"


def _cofactors(assembly: str = _TABLED_ASSEMBLY) -> CofactorTable:
    """The **Cofactor table** shipped for that assembly's species.

    Read off the shipped file rather than named, exactly as the census is: which genes a
    publisher lists as cofactors is the publisher's judgement, and no fixture may stand
    in for it.
    """
    species = assembly_metadata(assembly).species
    assert species is not None, f"{assembly} has no species in the assembly table"
    table = cofactor_table(species)
    assert table is not None, f"no cofactor table ships for {species}"
    return table


def _cofactor_row(stem: str, assembly: str = _TABLED_ASSEMBLY) -> dict[str, str | None]:
    """The table's own row for one **Gene id stem**, keyed by its own column names."""
    table = _cofactors(assembly)
    row = table.rows[table.gene_id_stems.index(stem)]
    return dict(zip(table.columns, row, strict=True))


class TestTFCofactorList:
    """``AnnotationRegistry.tf_cofactor_list`` — a shipped cofactor table meets an annotation.

    The counterpart of :class:`TestTFGeneList`, registering its fixture annotations the
    same way and asserting the same crossing: the table's **Gene id stem**s arriving as
    this annotation's own gene ids, what rode back unresolved, and that an assembly no
    published table can answer for raises rather than answering with nothing. The shipped
    table answers throughout; what it holds is its publisher's business.

    Nothing here prepares an assembly, fetches anything or reads the **Data dir**: every
    test registers one small GTF into ``tmp_path`` and asks, and the suite's autouse
    network guard is what says the answer never left the wheel.
    """

    def _registry(
        self,
        tmp_path: Path,
        *gene_ids: str,
        assembly: str = _TABLED_ASSEMBLY,
        name: str = "mine",
    ) -> AnnotationRegistry:
        """Register a GTF declaring ``gene_ids`` under ``name`` and open the registry."""
        return _registry_declaring(tmp_path, *gene_ids, assembly=assembly, name=name)

    def test_the_table_arrives_in_this_annotations_own_gene_ids(self, tmp_path: Path) -> None:
        # The whole point of the surface: the answer joins to a counts matrix keyed by this
        # annotation's ids, with no normalisation left for the caller.
        stems = _cofactors().cofactor_stems[:2]
        gene_ids = [_versioned(stem) for stem in stems]
        registry = self._registry(tmp_path, *gene_ids)

        answer = registry.tf_cofactor_list("mine")

        assert (answer.assembly, answer.annotation) == (_TABLED_ASSEMBLY, "mine")
        assert [entry.gene_id_stem for entry in answer.cofactors] == list(stems)
        assert answer.gene_ids == gene_ids

    def test_an_entry_carries_the_uniform_columns_and_the_publishers_own_beside_them(
        self, tmp_path: Path
    ) -> None:
        table = _cofactors()
        stem = table.cofactor_stems[0]
        registry = self._registry(tmp_path, _versioned(stem))

        entry = registry.tf_cofactor_list("mine").cofactors[0]

        cells = _cofactor_row(stem)
        assert (entry.symbol, entry.source) == (cells["symbol"], cells["source"])
        assert entry.is_cofactor is True
        # Everything the publisher classified it with, under that publisher's own
        # namespaced name — the AnimalTFDB family and the category joined onto it.
        assert dict(entry.classifications) == {
            name: cells[name] for name in table.columns[len(COFACTOR_UNIFORM_COLUMNS) :]
        }
        assert "animaltfdb_family" in entry.classifications

    def test_membership_travels_with_the_publishers_that_listed_the_gene(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path, _versioned(_cofactors().cofactor_stems[0]))

        answer = registry.tf_cofactor_list("mine")

        assert answer.provenance == _cofactors().provenance
        assert answer.provenance.sources
        assert all(source.publisher and source.pubmed_id for source in answer.provenance.sources)
        assert answer.species == assembly_metadata(_TABLED_ASSEMBLY).species

    def test_stems_this_annotation_carries_no_gene_for_ride_back_on_the_answer(
        self, tmp_path: Path
    ) -> None:
        table = _cofactors()
        carried, absent = table.cofactor_stems[0], table.cofactor_stems[1]
        registry = self._registry(tmp_path, _versioned(carried))

        answer = registry.tf_cofactor_list("mine")

        assert [entry.gene_id_stem for entry in answer.cofactors] == [carried]
        assert absent in answer.unresolved
        # Every stem the table lists is accounted for one way or the other, so what the
        # publisher holds and this annotation does not is visible rather than dropped.
        assert len(answer.cofactors) + len(answer.unresolved) == len(table.cofactor_stems)

    def test_a_stem_naming_two_gene_ids_answers_with_both(self, tmp_path: Path) -> None:
        # The collision is what is under test rather than the biology: two gene ids that
        # reduce to one stem, in the shape GENCODE's pseudoautosomal copies have, so a
        # resolver taking the first would hand back one of them without saying it chose.
        stem = _cofactors().cofactor_stems[0]
        first, second = _versioned(stem), f"{_versioned(stem)}_PAR_Y"
        registry = self._registry(tmp_path, second, first)

        answer = registry.tf_cofactor_list("mine")

        assert [entry.gene_ids for entry in answer.cofactors] == [(first, second)]
        assert answer.gene_ids == [first, second]

    def test_the_species_follows_the_assembly_and_never_the_ids_the_gtf_holds(
        self, tmp_path: Path
    ) -> None:
        # Mouse gene ids registered for a yeast assembly. Asking for one species' cofactors
        # while holding another's assembly is not expressible, so this is answered about the
        # assembly's own species and never about what is in the GTF.
        stem = _cofactors().cofactor_stems[0]
        registry = self._registry(tmp_path, _versioned(stem), assembly=_UNTABLED_ASSEMBLY)

        with pytest.raises(NoCofactorTableError) as excinfo:
            registry.tf_cofactor_list("mine")

        message = str(excinfo.value)
        assert str(assembly_metadata(_UNTABLED_ASSEMBLY).species) in message
        assert str(assembly_metadata(_TABLED_ASSEMBLY).species) in message

    @pytest.mark.parametrize("assembly", [_CHIMERA, "tiny"])
    def test_an_assembly_nothing_names_a_species_for_says_so_rather_than_guessing(
        self, tmp_path: Path, assembly: str
    ) -> None:
        stem = _cofactors().cofactor_stems[0]
        registry = self._registry(tmp_path, _versioned(stem), assembly=assembly)

        with pytest.raises(UnknownSpeciesError) as excinfo:
            registry.tf_cofactor_list("mine")

        message = str(excinfo.value)
        assert assembly in message
        assert "cofactor table" in message
        assert str(assembly_metadata(_TABLED_ASSEMBLY).species) in message

    def test_the_two_absences_are_told_apart_by_type_and_caught_together(
        self, tmp_path: Path
    ) -> None:
        # As the census half's pair is: *nobody has published a table for this species* and
        # *nothing says what species this is* are different answers, both lookups, and
        # neither is an empty collection.
        stem = _cofactors().cofactor_stems[0]
        untabled = self._registry(tmp_path / "yeast", _versioned(stem), assembly=_UNTABLED_ASSEMBLY)
        unnamed = self._registry(tmp_path / "chimera", _versioned(stem), assembly=_CHIMERA)

        with pytest.raises(LookupError) as no_table:
            untabled.tf_cofactor_list("mine")
        with pytest.raises(LookupError) as no_species:
            unnamed.tf_cofactor_list("mine")

        assert isinstance(no_table.value, NoCofactorTableError)
        assert isinstance(no_species.value, UnknownSpeciesError)
        assert not isinstance(no_table.value, UnknownSpeciesError)
        assert not isinstance(no_species.value, NoCofactorTableError)

    def test_the_worm_answers_here_although_the_census_half_raises_for_it(
        self, tmp_path: Path
    ) -> None:
        # The asymmetry, pinned: AnimalTFDB assessed worm cofactors and nobody has released
        # a worm TF census, so one assembly gets two different answers. That is the
        # publishers' shape and not a defect, and a test says so where it would otherwise
        # be filed as one.
        stem = _cofactors(_WORM_ASSEMBLY).cofactor_stems[0]
        registry = self._registry(tmp_path, _versioned(stem), assembly=_WORM_ASSEMBLY)

        answer = registry.tf_cofactor_list("mine")

        assert [entry.gene_id_stem for entry in answer.cofactors] == [stem]
        with pytest.raises(NoTFCensusError):
            registry.tf_gene_list("mine")

    def test_an_unregistered_name_earns_the_error_it_already_had(self, tmp_path: Path) -> None:
        registry = AnnotationRegistry.locate(_TABLED_ASSEMBLY, tmp_path)

        with pytest.raises(AnnotationNotRegisteredError) as excinfo:
            registry.tf_cofactor_list("gencode_vM39")

        assert f"genome register-annotation {_TABLED_ASSEMBLY} gencode_vM39" in str(excinfo.value)

    def test_naming_no_annotation_asks_the_default_one(self, tmp_path: Path) -> None:
        stem = _cofactors().cofactor_stems[0]
        registry = self._registry(tmp_path, _versioned(stem), name="gencode_vM39")

        assert registry.default == "gencode_vM39"
        assert registry.tf_cofactor_list().annotation == "gencode_vM39"

    def test_the_json_record_carries_the_cofactors_the_provenance_and_the_unresolved_stems(
        self, tmp_path: Path
    ) -> None:
        # What ``--json`` has to be able to emit, in the shape the TF gene list's answer
        # already uses: the entries with the publisher's own classification of each, the
        # provenance to cite, and the stems that resolved to nothing.
        stem = _cofactors().cofactor_stems[0]
        registry = self._registry(tmp_path, _versioned(stem))

        payload = registry.tf_cofactor_list("mine").as_json()

        assert payload["assembly"] == _TABLED_ASSEMBLY
        assert payload["gene_ids"] == [_versioned(stem)]
        assert payload["provenance"]["sources"][0]["pubmed_id"] == (
            _cofactors().provenance.sources[0].pubmed_id
        )
        assert payload["unresolved"]
        entry = payload["cofactors"][0]
        assert (entry["gene_id_stem"], entry["gene_ids"]) == (stem, [_versioned(stem)])
        assert entry["source"] == _cofactor_row(stem)["source"]
        assert (
            entry["classifications"]["animaltfdb_family"]
            == (_cofactor_row(stem)["animaltfdb_family"])
        )
        assert json.loads(json.dumps(payload)) == payload  # serializes as it stands

    def test_it_answers_for_an_assembly_named_rather_than_opened(self, tmp_path: Path) -> None:
        # One code path, asserted whole rather than sampled: the module-level function and
        # the method reach the same answer, so a shell surface over either says one thing.
        stem = _cofactors().cofactor_stems[0]
        registry = self._registry(tmp_path, _versioned(stem))

        answer = tf_cofactor_list(_TABLED_ASSEMBLY, annotation="mine", cache_dir=tmp_path)

        assert [entry.gene_id_stem for entry in answer.cofactors] == [stem]
        assert answer == registry.tf_cofactor_list("mine")


class TestAddressedByAssembly:
    """``gene_list`` and ``gene_lists`` addressed by assembly name rather than opened.

    A registry for the length of the call, exactly as ``annotation_status`` is, so there
    is no second code path to keep in step.
    """

    def test_it_answers_for_an_assembly_named_rather_than_opened(self, tmp_path: Path) -> None:
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(tmp_path, source, _CURATED, assembly=_CURATED_ASSEMBLY)
        category = _declared(_CURATED)[0]

        answer = gene_list(_CURATED_ASSEMBLY, category, cache_dir=tmp_path)

        assert answer.annotation == _CURATED
        assert [entry.category for entry in gene_lists(_CURATED_ASSEMBLY, cache_dir=tmp_path)] == (
            list(_declared(_CURATED))
        )

    def test_an_assembly_with_nothing_registered_raises_rather_than_answering(
        self, tmp_path: Path
    ) -> None:
        with pytest.raises(AnnotationNotRegisteredError):
            gene_list(_CURATED_ASSEMBLY, "rRNA", annotation=_CURATED, cache_dir=tmp_path)


# ---------------------------------------------------------------------------------------
# The edge this module must not grow back
# ---------------------------------------------------------------------------------------


def test_registering_an_annotation_imports_nothing_that_downloads_an_assembly() -> None:
    # The cycle, asserted closed. `io.download` imports `io.source` at the top of the
    # file, `io.chimera` imports `io.gtf`, and `io.gtf` used to import `io.download`
    # back — once for the annotations subdirectory name, which `io.registration` defines
    # and the downloader merely re-exports, and once for the package's one fetch step,
    # which `io.fetch` now holds. Each was a single line, and each grows back the moment
    # somebody reaches for a name that happens to be importable from the downloader.
    # Were the edge to return, asking what a chimera is made of would drag the whole
    # annotation build stack in behind it, and the downloader's module-level import of
    # the resolution would go back behind a deferred one.
    forbidden = {"genome.io.download", "genome.io.chimera", "genome.genome"}

    assert _module_level_imports(gtf_module) & forbidden == set()


def test_the_annotation_fetch_is_the_packages_one_fetch_step() -> None:
    # The positive half of the guard above: the edge is gone because the fetch moved to a
    # module of its own, not because this one started spelling a download itself.
    assert "genome.io.fetch" in _module_level_imports(gtf_module)
