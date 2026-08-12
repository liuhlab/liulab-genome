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
from pathlib import Path

import gffutils
import pytest

from genome.io.completion import (
    RegistrationMismatchError,
    UnfinishedRegistrationError,
    read_record,
    record_path,
    work_dir,
)
from genome.io.gtf import (
    annotation_dir,
    fetch_annotation,
    list_annotations,
    register_annotation,
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

#: sha256 of the committed ``tiny.gtf`` — the *unpacked* bytes ``tiny.gtf.gz`` yields.
_TINY_GTF_SHA256 = "255f43bd9abef76424d1c2d89a40cccc1a36215409bbc8f32dcead49ca3baf5e"

#: A URL that is nothing like any provider's, so using it can only come from a row.
_PINNED_URL = "https://mirror.example.invalid/annotations/tiny.gtf.gz"

#: The name the fixture annotation is registered under throughout.
_NAME = "ensgene_v101"


def _row(*, url: str = _PINNED_URL, sha256: str | None = None) -> AnnotationMetadata:
    """An in-memory annotation row for the ``tiny`` assembly."""
    return AnnotationMetadata(
        assembly="tiny",
        name=_NAME,
        provider="UCSC",
        version="ensGene.v101",
        url=url,
        sha256=sha256,
        default=True,
    )


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

    def test_it_files_the_annotation_under_the_assembly_data_dir(
        self, fake_fetch: FakeFetch, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        fake_fetch.serve("tiny.gtf.gz")
        monkeypatch.setenv("LIULAB_DATA", str(tmp_path))

        payload = register_annotation("tiny", _NAME, progressbar=False, metadata=_row())

        assert payload["directory"] == str(tmp_path / "genome" / "tiny" / "gtf" / _NAME)
