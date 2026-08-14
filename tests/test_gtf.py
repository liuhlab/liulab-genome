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
from pathlib import Path

import gffutils
import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.io.completion import (
    RegistrationMismatchError,
    UnfinishedRegistrationError,
    read_record,
    record_path,
    work_dir,
)
from genome.io.gtf import (
    ChromosomeMismatchError,
    _reject_unknown_chromosomes,
    annotation_dir,
    annotation_status,
    chromosome_check_summary,
    default_annotation,
    fetch_annotation,
    list_annotations,
    list_broken_annotations,
    register_annotation,
    register_annotation_by_path,
    register_gtf,
)
from genome.io.utils import ChecksumMismatchError
from genome.metadata import AnnotationMetadata

from .conftest import FakeFetch

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
    """``register_gtf`` — the escape hatch for a GTF the table does not list."""

    def test_a_plain_gtf_is_copied_built_and_recorded(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly = tmp_path / "asm"

        annotation = register_gtf(assembly, src, "WS298")

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

        annotation = register_gtf(assembly, src, "WS298")

        # Stored as a plain .gtf with decompressed contents, and the db builds.
        assert annotation.gtf.suffix == ".gtf"
        assert annotation.gtf.read_text() == _GTF
        assert annotation.db.is_file()
        assert list(list_annotations(assembly)) == ["WS298"]

    def test_a_missing_source_says_what_to_pass(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="GTF file not found"):
            register_gtf(tmp_path / "asm", tmp_path / "nope.gtf", "X")

    def test_reregistering_a_valid_one_returns_it_without_rebuilding(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly = tmp_path / "asm"
        first = register_gtf(assembly, src, "WS298")
        built_at = first.db.stat().st_mtime_ns

        # Silently: `filterwarnings = ["error"]` fails the test on any warning at all.
        second = register_gtf(assembly, src, "WS298")

        assert second == first
        assert second.db.stat().st_mtime_ns == built_at

    def test_force_rebuilds(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly = tmp_path / "asm"
        register_gtf(assembly, src, "WS298")

        annotation = register_gtf(assembly, src, "WS298", force=True)

        assert annotation.db.is_file()
        assert list(list_annotations(assembly)) == ["WS298"]

    def test_a_directory_without_a_record_raises_naming_its_repair(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly = tmp_path / "asm"
        directory = annotation_dir(assembly, "WS298")
        directory.mkdir(parents=True)
        (directory / "WS298.db").write_bytes(b"half a database")

        with pytest.raises(UnfinishedRegistrationError, match="force=True"):
            register_gtf(assembly, src, "WS298")


class TestRegisterByName:
    """``fetch_annotation`` — naming an annotation is enough to have it on disk."""

    @pytest.fixture(autouse=True)
    def _serve_the_gtf(self, fake_fetch: FakeFetch) -> FakeFetch:
        fake_fetch.serve("tiny.gtf.gz")
        return fake_fetch

    def test_it_fetches_verifies_builds_and_records(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        annotation = fetch_annotation(
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
        annotation = fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=_row())

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
            fetch_annotation(
                tmp_path, "tiny", _NAME, progressbar=False, metadata=_row(sha256=wrong)
            )

        assert wrong in str(excinfo.value)
        assert _TINY_GTF_SHA256 in str(excinfo.value)
        # Nothing that could not be vouched for reached the annotation's own files.
        directory = annotation_dir(tmp_path, _NAME)
        assert not (directory / f"{_NAME}.gtf").exists()
        assert read_record(directory) is None

    def test_a_row_that_pins_no_digest_records_whatever_arrived(self, tmp_path: Path) -> None:
        fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=_row())

        record = read_record(annotation_dir(tmp_path, _NAME))
        assert record is not None
        assert record.sha256 == _TINY_GTF_SHA256

    def test_an_uncompressed_url_is_placed_as_it_arrives(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        fake_fetch.serve("tiny.gtf")
        url = "https://mirror.example.invalid/annotations/tiny.gtf"

        annotation = fetch_annotation(
            tmp_path, "tiny", _NAME, progressbar=False, metadata=_row(url=url)
        )

        assert annotation.gtf.read_text().startswith("chrII\tensGene.v101\ttranscript")

    def test_a_name_no_row_lists_says_what_is_offered(self, tmp_path: Path) -> None:
        # Against the shipped table, which lists exactly one annotation for sacCer3.
        with pytest.raises(ValueError, match="no annotation named 'nope'") as excinfo:
            fetch_annotation(tmp_path, "sacCer3", "nope", progressbar=False)

        assert "ensgene_v101" in str(excinfo.value)
        assert "register_gtf" in str(excinfo.value)

    def test_reregistering_a_valid_one_is_a_silent_no_op(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        row = _row(sha256=_TINY_GTF_SHA256)
        first = fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)
        built_at = first.db.stat().st_mtime_ns

        # Silently: `filterwarnings = ["error"]` fails the test on any warning at all.
        second = fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)

        assert second == first
        assert second.db.stat().st_mtime_ns == built_at  # not rebuilt
        assert len(fake_fetch.calls) == 1  # nothing fetched twice

    def test_a_half_built_annotation_is_reported_as_broken(self, tmp_path: Path) -> None:
        # A gffutils build killed part-way: a database file, and no record.
        directory = annotation_dir(tmp_path, _NAME)
        directory.mkdir(parents=True)
        (directory / f"{_NAME}.db").write_bytes(b"half a database")

        with pytest.raises(UnfinishedRegistrationError) as excinfo:
            fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=_row())

        assert f"genome register-annotation tiny {_NAME} --force" in str(excinfo.value)

    def test_a_record_that_disagrees_with_disk_raises(self, tmp_path: Path) -> None:
        row = _row(sha256=_TINY_GTF_SHA256)
        annotation = fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)
        annotation.db.write_bytes(b"truncated")

        with pytest.raises(RegistrationMismatchError, match="disagrees with its"):
            fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)

    def test_force_repairs_what_the_error_named(self, tmp_path: Path) -> None:
        directory = annotation_dir(tmp_path, _NAME)
        directory.mkdir(parents=True)
        (directory / f"{_NAME}.db").write_bytes(b"half a database")

        annotation = fetch_annotation(
            tmp_path, "tiny", _NAME, progressbar=False, force=True, metadata=_row()
        )

        assert read_record(directory) is not None
        assert annotation.db.stat().st_size > len(b"half a database")

    def test_force_keeps_a_gtf_whose_digest_still_matches(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        row = _row(sha256=_TINY_GTF_SHA256)
        fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)
        record_path(annotation_dir(tmp_path, _NAME)).unlink()

        fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, force=True, metadata=row)

        assert len(fake_fetch.calls) == 1  # the GTF on disk proved itself; nothing refetched

    def test_force_refetches_when_the_row_pins_nothing_to_prove_it_against(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        row = _row()
        fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)

        fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, force=True, metadata=row)

        assert len(fake_fetch.calls) == 2

    def test_the_record_carries_gffutils_rather_than_a_tool_version(self, tmp_path: Path) -> None:
        # gffutils is a Python library, not an External tool resolved on PATH, so its
        # version is provenance in details and never a tool version.
        fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=_row())

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
        register_gtf(tmp_path, src, "finished")
        halfway = annotation_dir(tmp_path, "halfway")
        halfway.mkdir(parents=True)
        (halfway / "halfway.db").write_bytes(b"half a database")

        assert list(list_annotations(tmp_path)) == ["finished"]

    def test_a_record_that_disagrees_with_disk_is_not_registered(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        annotation = register_gtf(tmp_path, src, "WS298")
        assert list(list_annotations(tmp_path)) == ["WS298"]

        annotation.db.write_bytes(b"truncated")

        assert list_annotations(tmp_path) == {}


class TestListBrokenAnnotations:
    """The complement of ``list_annotations``: what is on disk and cannot be trusted.

    Every directory under ``gtf/`` is registered, broken, or not begun; this reports
    the middle one, which is otherwise invisible to anything that lists.
    """

    def test_a_directory_holding_files_with_no_record_is_broken(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        register_gtf(tmp_path, src, "mine")
        record_path(annotation_dir(tmp_path, "mine")).unlink()

        broken = list_broken_annotations(tmp_path, "tiny")

        assert list(broken) == ["mine"]
        assert broken["mine"].directory == annotation_dir(tmp_path, "mine")
        assert "holds files but no .completion.json" in broken["mine"].problem

    def test_a_record_that_disagrees_with_disk_is_broken(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        annotation = register_gtf(tmp_path, src, "mine")
        annotation.db.write_bytes(b"truncated")

        broken = list_broken_annotations(tmp_path, "tiny")

        assert list(broken) == ["mine"]
        assert "mine.db" in broken["mine"].problem

    def test_a_finished_annotation_is_not_broken(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        register_gtf(tmp_path, src, "mine")

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
        register_gtf(tmp_path, src, "healthy")
        register_gtf(tmp_path, src, "damaged")
        record_path(annotation_dir(tmp_path, "damaged")).unlink()

        assert list(list_annotations(tmp_path)) == ["healthy"]
        assert list(list_broken_annotations(tmp_path, "tiny")) == ["damaged"]

    def test_a_name_the_table_offers_is_repaired_by_name(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        register_gtf(tmp_path, src, "ensgene_v101")
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
        annotation = register_gtf(tmp_path, src, "mine")
        annotation.db.write_bytes(b"truncated")

        broken = list_broken_annotations(tmp_path, "tiny")

        assert broken["mine"].repair == f"genome register-gtf tiny {src} mine --force"

    def test_an_unlisted_one_whose_source_is_unknowable_says_so(self, tmp_path: Path) -> None:
        # No record survives to say which GTF it was built from, so there is no path to
        # print: the command is named with the one thing it still needs filled in,
        # rather than a path that would not run.
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        register_gtf(tmp_path, src, "mine")
        record_path(annotation_dir(tmp_path, "mine")).unlink()

        broken = list_broken_annotations(tmp_path, "tiny")

        assert broken["mine"].repair == "genome register-gtf tiny <path> mine --force"

    def test_an_unlisted_one_whose_source_is_gone_is_not_named_as_a_command(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        annotation = register_gtf(tmp_path, src, "mine")
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

    def test_it_reports_what_is_offered_with_nothing_registered(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # The case it most needs to serve: a fresh machine, where the answer is
        # entirely the shipped table's.
        monkeypatch.setenv("LIULAB_DATA", str(tmp_path))

        payload = annotation_status("sacCer3")

        assert payload["assembly"] == "sacCer3"
        assert payload["directory"] == str(tmp_path / "genome" / "sacCer3")
        assert payload["default_annotation"] == "ensgene_v101"
        rows = payload["annotations"]
        assert isinstance(rows, list)
        assert [(r["name"], r["offered"], r["registered"]) for r in rows] == [
            ("ensgene_v101", True, False)
        ]
        assert rows[0]["provider"] == "UCSC"
        assert rows[0]["path"] is None

    def test_it_creates_nothing_and_fetches_nothing(
        self, fake_fetch: FakeFetch, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("LIULAB_DATA", str(tmp_path))

        annotation_status("hg38")

        assert fake_fetch.calls == []
        assert not (tmp_path / "genome" / "hg38").exists()

    def test_a_registered_annotation_the_table_offers_is_reported_as_both(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly_dir = tmp_path / "asm"
        annotation = register_gtf(assembly_dir, src, "ensgene_v101")

        payload = annotation_status("sacCer3", cache_dir=assembly_dir)

        rows = payload["annotations"]
        assert isinstance(rows, list)
        assert [(r["name"], r["offered"], r["registered"]) for r in rows] == [
            ("ensgene_v101", True, True)
        ]
        assert rows[0]["path"] == str(annotation.gtf)

    def test_a_registered_annotation_no_row_lists_is_reported_too(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly_dir = tmp_path / "asm"
        register_gtf(assembly_dir, src, "mine")

        payload = annotation_status("tiny", cache_dir=assembly_dir)

        rows = payload["annotations"]
        assert isinstance(rows, list)
        assert [(r["name"], r["offered"], r["registered"]) for r in rows] == [("mine", False, True)]
        assert rows[0]["provider"] is None
        assert rows[0]["broken"] is False
        assert rows[0]["problem"] is None
        assert payload["default_annotation"] == "mine"  # nothing flagged, and it is alone

    def test_a_broken_offered_annotation_is_reported_as_broken_not_as_absent(
        self, tmp_path: Path
    ) -> None:
        # The bug this closes: half-registered and never-fetched looked identical here.
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly_dir = tmp_path / "asm"
        register_gtf(assembly_dir, src, "ensgene_v101")
        record_path(annotation_dir(assembly_dir, "ensgene_v101")).unlink()

        payload = annotation_status("sacCer3", cache_dir=assembly_dir)

        rows = payload["annotations"]
        assert isinstance(rows, list)
        assert [(r["name"], r["offered"], r["registered"], r["broken"]) for r in rows] == [
            ("ensgene_v101", True, False, True)
        ]
        assert rows[0]["repair"] == "genome register-annotation sacCer3 ensgene_v101 --force"
        assert "holds files but no .completion.json" in str(rows[0]["problem"])
        assert rows[0]["path"] is None

    def test_a_broken_unlisted_annotation_is_reported_at_all(self, tmp_path: Path) -> None:
        # No row lists it and no record vouches for it, so nothing used to mention it.
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly_dir = tmp_path / "asm"
        annotation = register_gtf(assembly_dir, src, "mine")
        annotation.db.write_bytes(b"truncated")

        payload = annotation_status("tiny", cache_dir=assembly_dir)

        rows = payload["annotations"]
        assert isinstance(rows, list)
        assert [(r["name"], r["offered"], r["registered"], r["broken"]) for r in rows] == [
            ("mine", False, False, True)
        ]
        assert rows[0]["repair"] == f"genome register-gtf tiny {src} mine --force"

    def test_one_broken_annotation_does_not_hide_the_others(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly_dir = tmp_path / "asm"
        register_gtf(assembly_dir, src, "healthy")
        register_gtf(assembly_dir, src, "damaged")
        record_path(annotation_dir(assembly_dir, "damaged")).unlink()

        payload = annotation_status("tiny", cache_dir=assembly_dir)

        rows = payload["annotations"]
        assert isinstance(rows, list)
        assert [(r["name"], r["registered"], r["broken"]) for r in rows] == [
            ("damaged", False, True),
            ("healthy", True, False),
        ]
        # A broken one is not a registered one, so it never becomes the sole-registered
        # default either.
        assert payload["default_annotation"] == "healthy"


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
            register_gtf(tmp_path, data_dir / "ensembl_style.gtf", _NAME, chrom_sizes=sizes)

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
            register_gtf(tmp_path, data_dir / "ensembl_style.gtf", _NAME, chrom_sizes=sizes)

        assert not annotation_dir(tmp_path, _NAME).exists()
        assert list(tmp_path.rglob("*.db")) == []

    def test_refused_by_name_it_stays_refused_rather_than_reading_as_interrupted(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        fake_fetch.serve("ensembl_style.gtf")
        _write_chrom_sizes(tmp_path, *self._UCSC)
        row = _row(url="https://mirror.example.invalid/annotations/ensembl_style.gtf")

        with pytest.raises(ChromosomeMismatchError):
            fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)

        directory = annotation_dir(tmp_path, _NAME)
        assert not (directory / f"{_NAME}.gtf").exists()  # never placed
        assert list(tmp_path.rglob("*.db")) == []  # never paid for the build
        assert read_record(directory) is None

        # Running it again reports the same problem, not an interrupted registration.
        with pytest.raises(ChromosomeMismatchError):
            fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)

    def test_sequences_the_annotation_never_mentions_are_not_an_error(self, tmp_path: Path) -> None:
        # Strict one way only: the GTF names chrI alone, the assembly carries five.
        sizes = _write_chrom_sizes(tmp_path, *self._UCSC, "chrM", "scaffold_17")
        source = tmp_path / "one-chromosome.gtf"
        source.write_text(_GTF)

        annotation = register_gtf(tmp_path, source, "WS298", chrom_sizes=sizes)

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
            register_gtf(tmp_path, source, _NAME, chrom_sizes=sizes)

        assert not annotation_dir(tmp_path, _NAME).exists()

    def test_a_wholesale_mismatch_lists_ten_names_and_counts_the_rest(self, tmp_path: Path) -> None:
        offenders = [f"scaffold_{n}" for n in range(25)]
        source = tmp_path / "many.gtf"
        source.write_text(
            "".join(f'{name}\ttest\texon\t1\t100\t.\t+\t.\tgene_id "g1";\n' for name in offenders)
        )
        sizes = _write_chrom_sizes(tmp_path, *self._UCSC)

        with pytest.raises(ChromosomeMismatchError) as excinfo:
            register_gtf(tmp_path, source, _NAME, chrom_sizes=sizes)

        assert len(excinfo.value.missing) == 25  # every one of them is on the exception
        assert "(and 15 more)" in str(excinfo.value)  # ten of them are in the message

    def test_comment_lines_are_not_taken_for_chromosomes(self, tmp_path: Path) -> None:
        sizes = _write_chrom_sizes(tmp_path, *self._UCSC)
        source = tmp_path / "commented.gtf"
        source.write_text("##description: a header\n#!genome-build tiny\n" + _GTF)

        assert register_gtf(tmp_path, source, "WS298", chrom_sizes=sizes).db.is_file()

    def test_the_override_registers_a_mismatched_gtf_anyway(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        sizes = _write_chrom_sizes(tmp_path, *self._UCSC)

        annotation = register_gtf(
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

        annotation = fetch_annotation(
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

        fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=_row())

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

        annotation = fetch_annotation(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)

        assert annotation.db.is_file()
        record = read_record(annotation_dir(tmp_path, _NAME))
        assert record is not None
        assert record.details["chromosomes_checked"] is False
        assert record.details["chromosomes_unchecked_because"] == "no-chrom-sizes"


class TestReadingBackWhatWasChecked:
    """``chromosome_check_summary`` — one sentence per state, and never the wrong one.

    The states differ in what a reader should do about them, which is why they are told
    apart at all: an annotation registered before its assembly is waiting for the
    assembly, and one whose check the caller stood down is waiting for nothing.
    """

    _ADVICE = "register the assembly first"

    def test_a_check_that_ran_says_so_rather_than_saying_nothing(self) -> None:
        # Silence is not how a pass is reported: a surface printing nothing about the
        # check reads exactly like one printing that it passed.
        summary = chromosome_check_summary(
            {"chromosomes_checked": True, "chromosomes_unchecked_because": None}
        )

        assert "chromosomes checked" in summary
        assert self._ADVICE not in summary

    def test_nothing_to_check_against_is_the_one_state_that_advises(self) -> None:
        summary = chromosome_check_summary(
            {"chromosomes_checked": False, "chromosomes_unchecked_because": "no-chrom-sizes"}
        )

        assert "chromosomes not checked" in summary
        assert self._ADVICE in summary

    def test_an_override_is_never_told_to_register_the_assembly(self) -> None:
        # The bug this fixes: the assembly may well be registered, and the caller turned
        # the check off on purpose. What is left to say is what the record does not vouch
        # for, not what to do about it.
        summary = chromosome_check_summary(
            {"chromosomes_checked": False, "chromosomes_unchecked_because": "caller-override"}
        )

        assert "stood down" in summary
        assert self._ADVICE not in summary

    def test_every_state_reads_as_its_own_sentence(self) -> None:
        summaries = {
            chromosome_check_summary(details)
            for details in (
                {"chromosomes_checked": True, "chromosomes_unchecked_because": None},
                {"chromosomes_checked": False, "chromosomes_unchecked_because": "no-chrom-sizes"},
                {"chromosomes_checked": False, "chromosomes_unchecked_because": "caller-override"},
                {"chromosomes_checked": False},
            )
        }

        assert len(summaries) == 4

    def test_a_record_written_before_the_reason_existed_reads_as_unknown(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        # The real back-compatibility case, on a record that is on disk: an older version
        # wrote the bare bool, and which of the two reasons it stood for is not knowable.
        # It must read as neither, and reading it must not raise.
        register_gtf(tmp_path, data_dir / "tiny.gtf", _NAME)
        path = record_path(annotation_dir(tmp_path, _NAME))
        written = json.loads(path.read_text())
        written["details"] = {"chromosomes_checked": False}
        path.write_text(json.dumps(written))

        record = read_record(annotation_dir(tmp_path, _NAME))

        assert record is not None
        assert record.details == {"chromosomes_checked": False}
        summary = chromosome_check_summary(record.details)
        assert "does not say why" in summary
        assert self._ADVICE not in summary  # nor is it claimed to be the override
        assert "stood down" not in summary

    def test_a_reason_this_version_has_never_heard_of_reads_as_unknown_too(self) -> None:
        # Forward as well as backward: a record from a later version claiming some third
        # reason is one this version cannot report, which is the same as not knowing.
        summary = chromosome_check_summary(
            {"chromosomes_checked": False, "chromosomes_unchecked_because": "some-later-reason"}
        )

        assert "does not say why" in summary


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

        assert payload["kind"] == "annotation"
        assert payload["name"] == _NAME
        assert payload["assembly"] == "tiny"
        assert payload["directory"] == str(annotation_dir(tmp_path, _NAME))
        assert payload["source_url"] == _PINNED_URL
        assert payload["sha256"] == _TINY_GTF_SHA256
        directory = annotation_dir(tmp_path, _NAME)
        assert payload["files"] == {
            name: (directory / name).stat().st_size for name in (f"{_NAME}.gtf", f"{_NAME}.db")
        }

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

        assert payload["details"] == {
            "provider": "UCSC",
            "version": "ensGene.v101",
            "gffutils_version": gffutils.__version__,
            "chromosomes_checked": False,
            "chromosomes_unchecked_because": "caller-override",
        }

    def test_it_files_the_annotation_under_the_assembly_data_dir(
        self, fake_fetch: FakeFetch, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        fake_fetch.serve("tiny.gtf.gz")
        monkeypatch.setenv("LIULAB_DATA", str(tmp_path))

        payload = register_annotation("tiny", _NAME, progressbar=False, metadata=_row())

        assert payload["directory"] == str(tmp_path / "genome" / "tiny" / "gtf" / _NAME)

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

        assert default["name"] == "exons_only"
        assert _feature_types(annotation_dir(tmp_path, "exons_only") / "exons_only.db") == ["exon"]
        assert inferred["name"] == "with_genes"
        assert _feature_types(annotation_dir(tmp_path, "with_genes") / "with_genes.db") == [
            "exon",
            "gene",
            "transcript",
        ]


class TestRegisterAnnotationByPath:
    """``register_annotation_by_path`` — a GTF no row lists, addressed by assembly name."""

    def test_it_returns_the_record_plus_where_it_landed(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("LIULAB_DATA", str(tmp_path))
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)

        payload = register_annotation_by_path("tiny", source, "WS298")

        directory = tmp_path / "genome" / "tiny" / "gtf" / "WS298"
        assert payload["kind"] == "annotation"
        assert payload["name"] == "WS298"
        assert payload["assembly"] == "tiny"
        assert payload["directory"] == str(directory)
        assert payload["source_url"] == str(source)
        assert payload["files"] == {
            name: (directory / name).stat().st_size for name in ("WS298.gtf", "WS298.db")
        }

    def test_cache_dir_overrides_which_assembly_directory_it_is_filed_under(
        self, tmp_path: Path
    ) -> None:
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        elsewhere = tmp_path / "elsewhere"

        payload = register_annotation_by_path("tiny", source, "WS298", cache_dir=elsewhere)

        assert payload["directory"] == str(annotation_dir(elsewhere, "WS298"))

    def test_it_finds_the_assembly_chrom_sizes_without_being_told(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        # What the by-directory form cannot do: given the assembly's name it knows where
        # its chrom.sizes is, so an Ensembl-spelled GTF is refused rather than silently
        # registered unchecked.
        _write_chrom_sizes(tmp_path, "chrI", "chrII", "chrIII")

        with pytest.raises(ChromosomeMismatchError):
            register_annotation_by_path(
                "tiny", data_dir / "ensembl_style.gtf", _NAME, cache_dir=tmp_path
            )

    def test_the_override_registers_the_mismatch_anyway_and_the_record_says_so(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        _write_chrom_sizes(tmp_path, "chrI", "chrII", "chrIII")

        payload = register_annotation_by_path(
            "tiny",
            data_dir / "ensembl_style.gtf",
            _NAME,
            cache_dir=tmp_path,
            check_chromosomes=False,
        )

        assert payload["details"] == {
            "gffutils_version": gffutils.__version__,
            "chromosomes_checked": False,
            "chromosomes_unchecked_because": "caller-override",
        }

    def test_the_inference_knobs_reach_the_database_build(self, tmp_path: Path) -> None:
        source = tmp_path / "bare.gtf"
        source.write_text(_BARE_GTF)

        register_annotation_by_path("tiny", source, "exons_only", cache_dir=tmp_path)
        register_annotation_by_path(
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
            register_annotation_by_path("tiny", tmp_path / "nope.gtf", "WS298", cache_dir=tmp_path)

    def test_a_broken_directory_names_the_command_that_repairs_it(self, tmp_path: Path) -> None:
        # Addressed by assembly name, the repair is a command a shell can run — the
        # by-directory form knows no assembly name and names the Python call instead.
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        directory = annotation_dir(tmp_path, "WS298")
        directory.mkdir(parents=True)
        (directory / "WS298.db").write_bytes(b"half a database")

        with pytest.raises(UnfinishedRegistrationError) as excinfo:
            register_annotation_by_path("tiny", source, "WS298", cache_dir=tmp_path)

        assert f"genome register-gtf tiny {source} WS298 --force" in str(excinfo.value)

    def test_force_repairs_what_the_error_named(self, tmp_path: Path) -> None:
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        directory = annotation_dir(tmp_path, "WS298")
        directory.mkdir(parents=True)
        (directory / "WS298.db").write_bytes(b"half a database")

        payload = register_annotation_by_path(
            "tiny", source, "WS298", cache_dir=tmp_path, force=True
        )

        assert payload["name"] == "WS298"
        assert list(list_annotations(tmp_path)) == ["WS298"]
