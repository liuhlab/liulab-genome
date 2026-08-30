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
    CompletionRecord,
    RegistrationMismatchError,
    UnfinishedRegistrationError,
    read_record,
    record_path,
    work_dir,
)
from genome.io.gtf import (
    UNCHECKED_CALLER_OVERRIDE,
    UNCHECKED_NO_CHROM_SIZES,
    AnnotationNotRegisteredError,
    AnnotationRegistry,
    AnnotationStatus,
    AnnotationStatusRow,
    ChromosomeMismatchError,
    GtfAnnotation,
    MergeSource,
    NoGeneFeaturesError,
    RegisteredAnnotation,
    _reject_unknown_chromosomes,
    annotation_dir,
    annotation_register_command,
    annotation_status,
    chromosome_check_summary,
    default_annotation,
    discard_merged_annotation,
    gene_list,
    gene_lists,
    list_annotations,
    list_broken_annotations,
    register_annotation,
    register_gtf,
    register_merged_gtf,
)
from genome.io.registration import AssemblyDir
from genome.io.utils import ChecksumMismatchError
from genome.metadata import AnnotationMetadata

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

    def test_a_plain_or_gzipped_gtf_is_copied_built_and_recorded(self, tmp_path: Path) -> None:
        plain = tmp_path / "ann.gtf"
        plain.write_text(_GTF)
        assembly = tmp_path / "asm"

        annotation = _register_by_path(assembly, plain, "WS298")

        assert annotation.gtf == annotation_dir(assembly, "WS298") / "WS298.gtf"
        assert annotation.gtf.read_text() == _GTF
        assert annotation.db.is_file()
        assert list(list_annotations(assembly)) == ["WS298"]

        record = read_record(annotation_dir(assembly, "WS298"))
        assert record is not None
        assert record.kind == "annotation"
        assert record.name == "WS298"
        assert sorted(record.files) == ["WS298.db", "WS298.gtf"]
        assert record.source_url == str(plain)

        # A gzipped source is decompressed on the way in — stored as a plain .gtf.
        gzipped = tmp_path / "ann.gtf.gz"
        with gzip.open(gzipped, "wt") as handle:
            handle.write(_GTF)
        gzip_assembly = tmp_path / "asm-gz"

        gzip_annotation = _register_by_path(gzip_assembly, gzipped, "WS298")

        assert gzip_annotation.gtf.suffix == ".gtf"
        assert gzip_annotation.gtf.read_text() == _GTF
        assert gzip_annotation.db.is_file()
        assert list(list_annotations(gzip_assembly)) == ["WS298"]

        with pytest.raises(FileNotFoundError, match="GTF file not found"):
            _register_by_path(tmp_path / "asm-missing", tmp_path / "nope.gtf", "X")

    def test_reregistering_is_a_no_op_unless_forced(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly = tmp_path / "asm"
        first = _register_by_path(assembly, src, "WS298")
        built_at = first.db.stat().st_mtime_ns

        # Silently: `filterwarnings = ["error"]` fails the test on any warning at all.
        second = _register_by_path(assembly, src, "WS298")
        assert second == first
        assert second.db.stat().st_mtime_ns == built_at  # not rebuilt

        forced = _register_by_path(assembly, src, "WS298", force=True)
        assert forced.db.is_file()
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
        # gffutils is a Python library, not an External tool resolved on PATH, so its
        # version is provenance in details and never a tool version.
        assert record.tool_versions == {}
        assert record.details["gffutils_version"] == gffutils.__version__
        assert record.details["provider"] == "UCSC"
        assert record.details["version"] == "ensGene.v101"

        # ...and a row that pins no digest at all still records whatever arrived.
        unpinned = _register_by_name(
            tmp_path, "tiny", "unpinned", progressbar=False, metadata=_row(name="unpinned")
        )
        unpinned_record = read_record(annotation_dir(tmp_path, "unpinned"))
        assert unpinned_record is not None
        assert unpinned_record.sha256 == _TINY_GTF_SHA256
        assert unpinned.db.is_file()

        # Against the shipped table, which lists exactly one annotation for sacCer3.
        with pytest.raises(ValueError, match="no annotation named 'nope'") as no_row:
            _register_by_name(tmp_path / "sacCer3", "sacCer3", "nope", progressbar=False)
        no_row_message = str(no_row.value)
        assert "ensgene_v101" in no_row_message
        assert "register-gtf" in no_row_message  # the way in for one no row lists
        assert "register_path" in no_row_message  # ...and the same from Python

    def test_the_database_it_builds_answers_queries_from_a_gzipped_or_plain_source(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        annotation = _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=_row())

        database = gffutils.FeatureDB(str(annotation.db))
        try:
            transcripts = list(database.features_of_type("transcript"))
            assert len(transcripts) == 18
            assert {feature.seqid for feature in transcripts} == {"chrI", "chrII", "chrIII"}
        finally:
            database.conn.close()

        # An uncompressed URL is placed exactly as it arrives, not run through gunzip.
        fake_fetch.serve("tiny.gtf")
        url = "https://mirror.example.invalid/annotations/tiny.gtf"
        plain = _register_by_name(
            tmp_path, "tiny", "plain", progressbar=False, metadata=_row(name="plain", url=url)
        )
        assert plain.gtf.read_text().startswith("chrII\tensGene.v101\ttranscript")

    def test_a_wrong_checksum_or_a_disk_disagreement_raises_naming_both_digests(
        self, tmp_path: Path
    ) -> None:
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

        # And once a good build has landed, a disk that no longer agrees with its own
        # record is the same kind of disagreement, raised the same way.
        row = _row(sha256=_TINY_GTF_SHA256)
        annotation = _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)
        annotation.db.write_bytes(b"truncated")
        with pytest.raises(RegistrationMismatchError, match="disagrees with its"):
            _register_by_name(tmp_path, "tiny", _NAME, progressbar=False, metadata=row)

    def test_reregistering_a_valid_one_is_a_silent_no_op_and_a_half_built_one_is_broken(
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

        # A gffutils build killed part-way: a database file, and no record.
        directory = annotation_dir(tmp_path / "half", _NAME)
        directory.mkdir(parents=True)
        (directory / f"{_NAME}.db").write_bytes(b"half a database")

        with pytest.raises(UnfinishedRegistrationError) as excinfo:
            _register_by_name(tmp_path / "half", "tiny", _NAME, progressbar=False, metadata=_row())
        assert f"genome register-annotation tiny {_NAME} --force" in str(excinfo.value)

        repaired = _register_by_name(
            tmp_path / "half", "tiny", _NAME, progressbar=False, force=True, metadata=_row()
        )
        assert read_record(directory) is not None
        assert repaired.db.stat().st_size > len(b"half a database")

    def test_force_keeps_a_matching_gtf_but_refetches_one_nothing_can_prove(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        pinned = _row(name="pinned", sha256=_TINY_GTF_SHA256)
        _register_by_name(tmp_path, "tiny", "pinned", progressbar=False, metadata=pinned)
        record_path(annotation_dir(tmp_path, "pinned")).unlink()

        _register_by_name(
            tmp_path, "tiny", "pinned", progressbar=False, force=True, metadata=pinned
        )
        assert len(fake_fetch.calls) == 1  # the GTF on disk proved itself; nothing refetched

        unpinned = _row(name="unpinned")
        _register_by_name(tmp_path, "tiny", "unpinned", progressbar=False, metadata=unpinned)
        _register_by_name(
            tmp_path, "tiny", "unpinned", progressbar=False, force=True, metadata=unpinned
        )
        assert len(fake_fetch.calls) == 3  # nothing on disk could prove itself; refetched


class TestListAnnotations:
    """Registered means a record that agrees with disk — never a database file."""

    def test_registered_means_a_record_agreeing_with_disk_never_just_a_database_file(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        annotation = _register_by_path(tmp_path, src, "finished")
        assert list(list_annotations(tmp_path)) == ["finished"]

        # A database file with no record beside it (a build killed part-way) is not
        # registered...
        halfway = annotation_dir(tmp_path, "halfway")
        halfway.mkdir(parents=True)
        (halfway / "halfway.db").write_bytes(b"half a database")
        assert list(list_annotations(tmp_path)) == ["finished"]

        # ...and neither is a record that no longer agrees with what is on disk.
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

    def test_discarding_returns_false_for_a_hand_registered_or_an_unregistered_name(
        self, tmp_path: Path
    ) -> None:
        # The name comes from a previous build's record, and a name is not ownership: only
        # a record showing a merge wrote it is.
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        _register_by_path(tmp_path, src, "a+b")
        assert discard_merged_annotation(tmp_path, "a+b") is False
        assert list(list_annotations(tmp_path)) == ["a+b"]

        # A name nothing is registered under at all is not an error either.
        assert discard_merged_annotation(tmp_path, "never-registered") is False


class TestListBrokenAnnotations:
    """The complement of ``list_annotations``: what is on disk and cannot be trusted.

    Every directory under ``gtf/`` is registered, broken, or not begun; this reports
    the middle one, which is otherwise invisible to anything that lists.
    """

    def test_broken_means_untrusted_files_never_a_fresh_or_finished_directory(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)

        _register_by_path(tmp_path, src, "no-record")
        record_path(annotation_dir(tmp_path, "no-record")).unlink()
        broken = list_broken_annotations(tmp_path, "tiny")
        assert list(broken) == ["no-record"]
        assert broken["no-record"].directory == annotation_dir(tmp_path, "no-record")
        assert "holds files but no .completion.json" in broken["no-record"].problem

        disagreeing = _register_by_path(tmp_path, src, "disagreeing")
        disagreeing.db.write_bytes(b"truncated")
        broken = list_broken_annotations(tmp_path, "tiny")
        assert "disagreeing.db" in broken["disagreeing"].problem

        # A finished annotation is never reported broken...
        _register_by_path(tmp_path / "finished", src, "mine")
        assert list_broken_annotations(tmp_path / "finished", "tiny") == {}

        # ...and neither is a fresh or empty one: ADR-0007 says an absent or empty
        # directory is a fresh registration, not a broken one. A run interrupted before
        # it downloaded anything must not be reported.
        annotation_dir(tmp_path, "fresh").mkdir(parents=True)
        work_dir(annotation_dir(tmp_path, "fresh")).mkdir()
        assert "fresh" not in list_broken_annotations(tmp_path, "tiny")

        # And an assembly with no annotations at all answers empty rather than raising.
        assert list_broken_annotations(tmp_path / "no-such-directory", "sacCer3") == {}

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

    def test_the_repair_it_names_depends_on_whether_the_table_or_a_path_can_be_shown(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)

        # A name the table offers is repaired by name.
        _register_by_path(tmp_path, src, "ensgene_v101")
        record_path(annotation_dir(tmp_path, "ensgene_v101")).unlink()
        offered = list_broken_annotations(tmp_path, "sacCer3")
        assert offered["ensgene_v101"].repair == (
            "genome register-annotation sacCer3 ensgene_v101 --force"
        )
        assert offered["ensgene_v101"].repair in offered["ensgene_v101"].problem

        # An unlisted one is repaired from the path its own record remembers.
        unlisted = _register_by_path(tmp_path, src, "mine")
        unlisted.db.write_bytes(b"truncated")
        broken = list_broken_annotations(tmp_path, "tiny")
        assert broken["mine"].repair == f"genome register-gtf tiny {src} mine --force"

        # No record survives to say which GTF it was built from, so there is no path to
        # print: the command is named with the one thing it still needs filled in,
        # rather than a path that would not run.
        _register_by_path(tmp_path, src, "unknowable")
        record_path(annotation_dir(tmp_path, "unknowable")).unlink()
        broken = list_broken_annotations(tmp_path, "tiny")
        assert broken["unknowable"].repair == "genome register-gtf tiny <path> unknowable --force"

        # A record survives, but the path it names is gone — same placeholder, for the
        # same reason: the path it remembers would not run either.
        gone = _register_by_path(tmp_path, src, "gone")
        gone.db.write_bytes(b"truncated")
        src.unlink()
        broken = list_broken_annotations(tmp_path, "tiny")
        assert str(src) not in broken["gone"].repair
        assert broken["gone"].repair == "genome register-gtf tiny <path> gone --force"


class TestDefaultAnnotation:
    """The one rule that decides a default, wherever the question is asked from."""

    def test_the_flag_the_sole_registration_or_an_explicit_choice_decides_it(self) -> None:
        assert default_annotation([_row()], []) == _NAME
        # Everyone in the lab reaches for the same one, whatever this machine happens
        # to hold...
        assert default_annotation([_row()], ["something_else"]) == _NAME
        # ...unless a caller explicitly chooses another, which wins over the flag.
        assert default_annotation([_row()], [_NAME], explicit="something_else") == "something_else"

        # With nothing flagged, the sole registered annotation wins, or there is none.
        assert default_annotation([_row(default=False)], ["only_one"]) == "only_one"
        assert default_annotation([], []) is None
        assert default_annotation([], ["one", "two"]) is None


class TestAnnotationStatus:
    """What an assembly's table offers, set against what is registered on this machine."""

    def test_a_fresh_machine_reports_the_shipped_table_serialized_and_untouched(
        self, fake_fetch: FakeFetch, liulab_data: Path
    ) -> None:
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

        # `--json` is this report rendered, so a row's fields and the payload's keys are
        # one spelling: a surface reads attributes and never names a key of its own.
        assert payload.as_json() == {
            "assembly": "sacCer3",
            "directory": str(liulab_data / "genome" / "sacCer3"),
            "default_annotation": "ensgene_v101",
            "annotations": [asdict(row) for row in payload.annotations],
        }

        # What the closing line of `genome annotations` needs: the default's own state,
        # so "not registered here" and "broken here" are told apart by the report itself.
        default = payload.default_row
        assert default is payload.annotations[0]
        assert default is not None
        assert not default.registered

        # No default decided, so no row is about one — the state a fresh unlisted
        # assembly is in, and not an error.
        nothing = annotation_status("tiny")
        assert nothing.default_annotation is None
        assert nothing.default_row is None

        # And asking creates nothing and fetches nothing.
        annotation_status("hg38")
        assert fake_fetch.calls == []
        assert not (liulab_data / "genome" / "hg38").exists()

    def test_a_registered_or_broken_annotation_is_reported_whether_offered_or_not(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        offered_dir, unlisted_dir = tmp_path / "offered", tmp_path / "unlisted"
        annotation = _register_by_path(offered_dir, src, "ensgene_v101")
        _register_by_path(unlisted_dir, src, "mine")

        offered = annotation_status("sacCer3", cache_dir=offered_dir)
        unlisted = annotation_status("tiny", cache_dir=unlisted_dir)

        offered_rows = offered.annotations
        assert [(r.name, r.offered, r.registered) for r in offered_rows] == [
            ("ensgene_v101", True, True)
        ]
        assert offered_rows[0].path == str(annotation.gtf)

        unlisted_rows = unlisted.annotations
        assert [(r.name, r.offered, r.registered) for r in unlisted_rows] == [("mine", False, True)]
        assert unlisted_rows[0].provider is None
        assert unlisted_rows[0].broken is False
        assert unlisted_rows[0].problem is None
        assert unlisted.default_annotation == "mine"  # nothing flagged, and it is alone

        # The bug this closes: half-registered and never-fetched looked identical here.
        # No row lists the second one and no record vouches for it either, so nothing
        # used to mention it at all — a broken annotation is reported as broken, and
        # never as simply absent, whether it is offered or not.
        record_path(annotation_dir(offered_dir, "ensgene_v101")).unlink()
        broken_unlisted = _register_by_path(unlisted_dir, src, "gone")
        broken_unlisted.db.write_bytes(b"truncated")

        broken_offered = annotation_status("sacCer3", cache_dir=offered_dir)
        broken_unlisted_payload = annotation_status("tiny", cache_dir=unlisted_dir)

        broken_offered_rows = broken_offered.annotations
        assert [(r.name, r.offered, r.registered, r.broken) for r in broken_offered_rows] == [
            ("ensgene_v101", True, False, True)
        ]
        assert broken_offered_rows[0].repair == (
            "genome register-annotation sacCer3 ensgene_v101 --force"
        )
        assert "holds files but no .completion.json" in str(broken_offered_rows[0].problem)
        assert broken_offered_rows[0].path is None

        unlisted_rows = [row for row in broken_unlisted_payload.annotations if row.name == "gone"]
        assert [(r.name, r.offered, r.registered, r.broken) for r in unlisted_rows] == [
            ("gone", False, False, True)
        ]
        assert unlisted_rows[0].repair == f"genome register-gtf tiny {src} gone --force"

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

    def test_the_four_states_are_settled_in_one_construction_bound_to_its_own_directory(
        self, tmp_path: Path
    ) -> None:
        self._registered(tmp_path, "ensgene_v101")
        self._registered(tmp_path, "damaged")
        record_path(annotation_dir(tmp_path, "damaged")).unlink()

        registry = AnnotationRegistry.locate("sacCer3", tmp_path)

        assert registry.registered == ["ensgene_v101"]
        assert [entry.name for entry in registry.broken] == ["damaged"]
        assert [record.name for record in registry.offered] == ["ensgene_v101"]
        assert registry.default == "ensgene_v101"

        # The assembly dir travels with the registry rather than being re-derived from
        # the data root at each question.
        elsewhere = tmp_path / "elsewhere"
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(elsewhere, source, "mine")

        assert AnnotationRegistry.locate("tiny", elsewhere).registered == ["mine"]
        assert AnnotationRegistry.locate("tiny", tmp_path / "unused").registered == []

        # Its status is exactly what the status report answers, and nothing is
        # created by asking.
        assert registry.status() == annotation_status("sacCer3", cache_dir=tmp_path)
        nothing = AnnotationRegistry.locate("sacCer3", tmp_path / "unused" / "genome" / "sacCer3")
        assert nothing.registered == []
        assert nothing.broken == []
        assert nothing.default == "ensgene_v101"
        assert not (tmp_path / "unused" / "genome").exists()

    def test_path_resolves_registered_names_and_explains_unregistered_or_broken_ones(
        self, tmp_path: Path
    ) -> None:
        self._registered(tmp_path, "mine")
        registry = AnnotationRegistry.locate("tiny", tmp_path)
        assert registry.path("mine") == annotation_dir(tmp_path, "mine") / "mine.gtf"

        offering = AnnotationRegistry.locate("sacCer3", tmp_path)
        with pytest.raises(AnnotationNotRegisteredError) as unregistered:
            offering.path("no_such_annotation")
        message = str(unregistered.value)
        assert "no_such_annotation" in message
        assert "ensgene_v101" in message  # what the table does offer

        # Not the command that would itself raise and demand --force: the one that works.
        self._registered(tmp_path, "ensgene_v101")
        record_path(annotation_dir(tmp_path, "ensgene_v101")).unlink()
        with pytest.raises(AnnotationNotRegisteredError) as broken:
            AnnotationRegistry.locate("sacCer3", tmp_path).path("ensgene_v101")
        assert "genome register-annotation sacCer3 ensgene_v101 --force" in str(broken.value)

    def test_registering_adopts_the_result_and_the_default_is_the_flag_unless_overridden(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        by_path = AnnotationRegistry.locate("tiny", tmp_path / "path")
        assert by_path.registered == []

        annotation = by_path.register_path(source, "mine")

        assert by_path.registered == ["mine"]
        assert by_path.path("mine") == annotation.gtf
        assert by_path.default == "mine"  # nothing flagged, and it is alone

        fake_fetch.serve("tiny.gtf.gz")
        by_name = AnnotationRegistry.locate("tiny", tmp_path / "name")

        fetched = by_name.register(_NAME, progressbar=False, metadata=_row())

        assert by_name.registered == [_NAME]
        assert by_name.path(_NAME) == fetched.gtf

        flagged = AnnotationRegistry.locate("sacCer3", tmp_path / "flagged")
        assert flagged.default == "ensgene_v101"  # the table's flag

        flagged.register_path(source, "mine")
        assert flagged.default == "ensgene_v101"  # never displaced by a registration

        explicit = AnnotationRegistry.locate("sacCer3", tmp_path / "explicit", default="mine")
        assert explicit.default == "mine"
        assert explicit.registered == []  # need not be registered to win

    def test_a_broken_directory_either_names_a_repair_command_or_is_repaired_by_force(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        source = self._registered(tmp_path, "mine")
        record_path(annotation_dir(tmp_path, "mine")).unlink()
        registry = AnnotationRegistry.locate("tiny", tmp_path)
        assert [entry.name for entry in registry.broken] == ["mine"]

        registry.register_path(source, "mine", force=True)
        assert registry.broken == []
        assert registry.registered == ["mine"]

        # A registry always knows its assembly, so the repair a half-built directory
        # names is a command a shell can run rather than a Python call with the
        # assembly left to guess at.
        directory = annotation_dir(tmp_path, "WS298")
        directory.mkdir(parents=True)
        (directory / "WS298.db").write_bytes(b"half a database")
        with pytest.raises(UnfinishedRegistrationError) as excinfo:
            registry.register_path(source, "WS298")
        assert f"genome register-gtf tiny {source} WS298 --force" in str(excinfo.value)

        # Chrom.sizes defaults to the assembly's own...
        chrom_dir = tmp_path / "chrom"
        _write_chrom_sizes(chrom_dir, "chrI", "chrII", "chrIII")
        with pytest.raises(ChromosomeMismatchError):
            AnnotationRegistry.locate("tiny", chrom_dir).register_path(
                data_dir / "ensembl_style.gtf", _NAME
            )

        # ...but what a Genome opened somewhere of its own passes — the chrom.sizes it
        # actually prepared, rather than the one the layout would name for its
        # assembly — is the one checked against.
        elsewhere = tmp_path / "elsewhere"
        sizes = _write_chrom_sizes(elsewhere, "chrI", "chrII", "chrIII", assembly="elsewhere")
        with pytest.raises(ChromosomeMismatchError):
            AnnotationRegistry(
                AssemblyDir.locate("tiny", elsewhere), chrom_sizes=sizes
            ).register_path(data_dir / "ensembl_style.gtf", _NAME)


class TestChromosomeNames:
    """A GTF's sequence names must be the assembly's, and are checked before the build.

    The mismatch case is the committed ``ensembl_style.gtf`` — ``tiny.gtf``'s own 85
    features with the ``chr`` prefix stripped — against a ``chrI``/``chrII``/``chrIII``
    assembly: the UCSC-versus-Ensembl case in real bytes.
    """

    #: How the fixture assembly spells its three sequences.
    _UCSC = ("chrI", "chrII", "chrIII")

    def test_a_mismatch_is_refused_naming_every_offender_and_costs_nothing_by_path_or_name(
        self, fake_fetch: FakeFetch, tmp_path: Path, data_dir: Path
    ) -> None:
        sizes = _write_chrom_sizes(tmp_path, *self._UCSC)

        with pytest.raises(ChromosomeMismatchError) as excinfo:
            _register_by_path(tmp_path, data_dir / "ensembl_style.gtf", _NAME, chrom_sizes=sizes)

        assert excinfo.value.missing == ("I", "II", "III")
        message = str(excinfo.value)
        assert "I, II, III" in message
        assert "chrI" in message  # what the assembly spells them as
        assert "check_chromosomes=False" in message
        # The check runs before the build, and before anything is placed: the
        # annotation directory is left exactly as it was found.
        assert not annotation_dir(tmp_path, _NAME).exists()
        assert list(tmp_path.rglob("*.db")) == []

        by_name_dir = tmp_path / "by-name"
        fake_fetch.serve("ensembl_style.gtf")
        _write_chrom_sizes(by_name_dir, *self._UCSC)
        row = _row(url="https://mirror.example.invalid/annotations/ensembl_style.gtf")

        with pytest.raises(ChromosomeMismatchError):
            _register_by_name(by_name_dir, "tiny", _NAME, progressbar=False, metadata=row)

        directory = annotation_dir(by_name_dir, _NAME)
        assert not (directory / f"{_NAME}.gtf").exists()  # never placed
        assert list(by_name_dir.rglob("*.db")) == []  # never paid for the build
        assert read_record(directory) is None

        # Running it again reports the same problem, not an interrupted registration.
        with pytest.raises(ChromosomeMismatchError):
            _register_by_name(by_name_dir, "tiny", _NAME, progressbar=False, metadata=row)

    def test_extra_sequences_comments_and_a_gzipped_source_are_handled_and_mismatches_summarized(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        # Strict one way only: the GTF names chrI alone, the assembly carries five.
        sizes = _write_chrom_sizes(tmp_path, *self._UCSC, "chrM", "scaffold_17")
        source = tmp_path / "one-chromosome.gtf"
        source.write_text(_GTF)

        annotation = _register_by_path(tmp_path, source, "WS298", chrom_sizes=sizes)

        assert annotation.db.is_file()
        record = read_record(annotation_dir(tmp_path, "WS298"))
        assert record is not None
        assert record.details["chromosomes_checked"] is True

        # A header comment is not a chromosome name either.
        commented_dir = tmp_path / "commented"
        commented_sizes = _write_chrom_sizes(commented_dir, *self._UCSC)
        commented = commented_dir / "commented.gtf"
        commented.write_text("##description: a header\n#!genome-build tiny\n" + _GTF)
        assert _register_by_path(
            commented_dir, commented, "WS298", chrom_sizes=commented_sizes
        ).db.is_file()

        # A gzipped source is checked without being unpacked first.
        gzip_dir = tmp_path / "gzip"
        gzipped = gzip_dir / "ensembl_style.gtf.gz"
        gzip_dir.mkdir()
        with gzip.open(gzipped, "wt") as handle:
            handle.write((data_dir / "ensembl_style.gtf").read_text())
        gzip_sizes = _write_chrom_sizes(gzip_dir, *self._UCSC)
        with pytest.raises(ChromosomeMismatchError):
            _register_by_path(gzip_dir, gzipped, _NAME, chrom_sizes=gzip_sizes)
        assert not annotation_dir(gzip_dir, _NAME).exists()

        # A wholesale mismatch lists ten names in the message and counts the rest.
        offenders = [f"scaffold_{n}" for n in range(25)]
        many_dir = tmp_path / "many"
        many_sizes = _write_chrom_sizes(many_dir, *self._UCSC)
        many = many_dir / "many.gtf"
        many.write_text(
            "".join(f'{name}\ttest\texon\t1\t100\t.\t+\t.\tgene_id "g1";\n' for name in offenders)
        )
        with pytest.raises(ChromosomeMismatchError) as excinfo:
            _register_by_path(many_dir, many, "WS299", chrom_sizes=many_sizes)
        assert len(excinfo.value.missing) == 25  # every one of them is on the exception
        assert "(and 15 more)" in str(excinfo.value)  # ten of them are in the message

    def test_the_override_registers_a_mismatch_by_path_or_name_and_the_record_says_why(
        self, fake_fetch: FakeFetch, tmp_path: Path, data_dir: Path
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

        by_name_dir = tmp_path / "by-name"
        fake_fetch.serve("ensembl_style.gtf")
        _write_chrom_sizes(by_name_dir, *self._UCSC)
        row = _row(url="https://mirror.example.invalid/annotations/ensembl_style.gtf")

        by_name = _register_by_name(
            by_name_dir, "tiny", _NAME, progressbar=False, metadata=row, check_chromosomes=False
        )

        assert by_name.db.is_file()
        by_name_record = read_record(annotation_dir(by_name_dir, _NAME))
        assert by_name_record is not None
        assert by_name_record.details["chromosomes_checked"] is False
        assert by_name_record.details["chromosomes_unchecked_because"] == "caller-override"

        matching_dir = tmp_path / "matching"
        fake_fetch.serve("tiny.gtf.gz")
        _write_chrom_sizes(matching_dir, *self._UCSC)
        _register_by_name(matching_dir, "tiny", _NAME, progressbar=False, metadata=_row())
        matching_record = read_record(annotation_dir(matching_dir, _NAME))
        assert matching_record is not None
        assert matching_record.details["chromosomes_checked"] is True
        # A check that ran and did not raise passed, so there is no reason beside it.
        assert matching_record.details["chromosomes_unchecked_because"] is None

        # An annotation registered before its assembly was prepared: no chrom.sizes
        # exists, so the names cannot be checked. The record says they were not, and
        # says it was for want of that file rather than because anyone asked to skip it.
        unchecked_dir = tmp_path / "unchecked"
        fake_fetch.serve("ensembl_style.gtf")
        row = _row(url="https://mirror.example.invalid/annotations/ensembl_style.gtf")
        annotation = _register_by_name(
            unchecked_dir, "tiny", _NAME, progressbar=False, metadata=row
        )
        assert annotation.db.is_file()
        unchecked_record = read_record(annotation_dir(unchecked_dir, _NAME))
        assert unchecked_record is not None
        assert unchecked_record.details["chromosomes_checked"] is False
        assert unchecked_record.details["chromosomes_unchecked_because"] == "no-chrom-sizes"


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

    def test_it_returns_and_serializes_the_record_plus_where_it_landed(
        self, fake_fetch: FakeFetch, tmp_path: Path, liulab_data: Path
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

        # The `--json` payload is the completion record under its own on-disk key
        # names, with the two facts a record does not hold about itself. A type wraps
        # those names; it never renames them, because lab directories are read by both.
        assert payload.as_json() == {
            **asdict(payload.record),
            "assembly": "tiny",
            "directory": str(directory),
        }
        assert list(payload.as_json())[-2:] == ["assembly", "directory"]

        # With no cache_dir given, it files under the assembly's own Data dir.
        default_payload = register_annotation(
            "tiny", "elsewhere", progressbar=False, metadata=_row(name="elsewhere")
        )
        assert default_payload.directory == liulab_data / "genome" / "tiny" / "gtf" / "elsewhere"

    def test_the_chromosome_check_is_read_off_the_record_from_both_ways_in(
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

        mismatch_dir = tmp_path / "mismatch"
        fake_fetch.serve("ensembl_style.gtf")
        _write_chrom_sizes(mismatch_dir, "chrI", "chrII", "chrIII")
        row = _row(url="https://mirror.example.invalid/annotations/ensembl_style.gtf")

        with pytest.raises(ChromosomeMismatchError):
            register_annotation(
                "tiny", _NAME, cache_dir=mismatch_dir, progressbar=False, metadata=row
            )

        override_payload = register_annotation(
            "tiny",
            _NAME,
            cache_dir=mismatch_dir,
            progressbar=False,
            metadata=row,
            check_chromosomes=False,
        )
        assert override_payload.record.details == {
            "provider": "UCSC",
            "version": "ensGene.v101",
            "gffutils_version": gffutils.__version__,
            "chromosomes_checked": False,
            "chromosomes_unchecked_because": "caller-override",
        }

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

    def test_it_returns_the_record_plus_where_it_landed_and_checks_chromosomes(
        self, tmp_path: Path, liulab_data: Path, data_dir: Path
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

        elsewhere = tmp_path / "elsewhere"
        overridden = register_gtf("tiny", source, "WS298", cache_dir=elsewhere)
        assert overridden.directory == annotation_dir(elsewhere, "WS298")

        # Naming the assembly is what says where its chrom.sizes is, so an
        # Ensembl-spelled GTF is refused rather than silently registered unchecked.
        mismatch_dir = tmp_path / "mismatch"
        _write_chrom_sizes(mismatch_dir, "chrI", "chrII", "chrIII")

        with pytest.raises(ChromosomeMismatchError):
            register_gtf("tiny", data_dir / "ensembl_style.gtf", _NAME, cache_dir=mismatch_dir)

        override_payload = register_gtf(
            "tiny",
            data_dir / "ensembl_style.gtf",
            _NAME,
            cache_dir=mismatch_dir,
            check_chromosomes=False,
        )
        assert override_payload.record.details == {
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

        with pytest.raises(FileNotFoundError, match="GTF file not found"):
            register_gtf("tiny", tmp_path / "nope.gtf", "WS298", cache_dir=tmp_path)

    def test_a_broken_directory_names_its_repair_and_force_applies_it(self, tmp_path: Path) -> None:
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
        tmp_path.mkdir(parents=True, exist_ok=True)
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(tmp_path, source, _CURATED, assembly=assembly)
        return AnnotationRegistry.locate(assembly, tmp_path)

    def test_a_plain_annotation_answers_with_one_unattributed_source_and_the_curated_ids(
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
        # The two sentences travel with the ids, because they are what says whether
        # these ids mean what the caller's metric needs.
        assert answer.sources[0].description.strip()
        assert answer.sources[0].source.strip()
        assert answer.gene_ids == list(_curated_ids(_CURATED, category))
        assert answer.gene_ids  # never an empty answer: a declared category has genes

    def test_the_two_kinds_of_absence_are_told_apart_and_each_names_whats_available(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path)

        with pytest.raises(GeneCategoryNotDeclaredError) as declared:
            registry.gene_list("no_such_category", _CURATED)
        message = str(declared.value)
        assert "no_such_category" in message
        for category in _declared(_CURATED):
            assert category in message

        # The fact #111 exists for: *no categories are declared* and *this category is
        # not declared* are different, and neither is an empty answer.
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(tmp_path, source, "mine")
        nothing_declared_registry = AnnotationRegistry.locate("tiny", tmp_path)
        with pytest.raises(NoGeneCategoriesError) as nothing_declared:
            nothing_declared_registry.gene_list("rRNA", "mine")
        nothing_message = str(nothing_declared.value)
        assert "mine" in nothing_message
        assert _CURATED in nothing_message  # …and which annotations do declare categories

        assert isinstance(declared.value, GeneCategoryNotDeclaredError)
        assert isinstance(nothing_declared.value, NoGeneCategoriesError)
        assert not isinstance(declared.value, NoGeneCategoriesError)
        assert not isinstance(nothing_declared.value, GeneCategoryNotDeclaredError)

    def test_naming_no_annotation_asks_the_default_and_an_unregistered_or_missing_one_names_it(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path / "default")
        category = _declared(_CURATED)[0]

        assert registry.default == _CURATED
        assert registry.gene_list(category).annotation == _CURATED

        # Callers pass names, never paths, and a name nothing registered is resolved by
        # the same `path` every other question goes through — so the message is the
        # same one.
        unregistered = AnnotationRegistry.locate(_CURATED_ASSEMBLY, tmp_path / "unregistered")
        with pytest.raises(AnnotationNotRegisteredError) as excinfo:
            unregistered.gene_list("rRNA", _CURATED)
        assert f"genome register-annotation {_CURATED_ASSEMBLY} {_CURATED}" in str(excinfo.value)

        no_default = AnnotationRegistry.locate("tiny", tmp_path / "no-default")
        with pytest.raises(ValueError, match="annotation") as no_default_excinfo:
            no_default.gene_list("rRNA")
        assert "default_gtf" in str(no_default_excinfo.value)

        # A name is unique only within its assembly, so a list found by name alone is not
        # yet known to be about this reference — and answering would hand back another
        # species' genes under this one's name.
        wrong_assembly = self._registry(tmp_path / "wrong", assembly="tiny")

        with pytest.raises(GeneListAssemblyMismatchError) as excinfo:
            wrong_assembly.gene_list(_declared(_CURATED)[0], _CURATED)

        message = str(excinfo.value)
        assert "tiny" in message
        assert _CURATED_ASSEMBLY in message

    def test_gene_lists_returns_every_declared_category_or_raises_if_none_are(
        self, tmp_path: Path
    ) -> None:
        answers = self._registry(tmp_path).gene_lists(_CURATED)
        assert [answer.category for answer in answers] == list(_declared(_CURATED))
        assert all(answer.gene_ids for answer in answers)

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
        tmp_path.mkdir(parents=True, exist_ok=True)
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_merged(tmp_path, name, source)
        return AnnotationRegistry.locate(_CHIMERA, tmp_path)

    def test_each_contributor_answers_for_its_own_component_and_ids_are_never_deduped(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path)
        shared = next(category for category in _declared(_WORM) if category in _declared(_FOOD))

        answer = registry.gene_list(shared, f"{_WORM}+{_FOOD}")

        assert [(source.component, source.annotation) for source in answer.sources] == [
            (_WORM_COMPONENT, _WORM),
            (_FOOD_COMPONENT, _FOOD),
        ]
        assert answer.annotation == f"{_WORM}+{_FOOD}"
        assert answer.assembly == _CHIMERA
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

    def test_a_category_no_contributor_declares_raises_and_gene_lists_is_their_union(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path)

        with pytest.raises(GeneCategoryNotDeclaredError):
            registry.gene_list("no_such_category", f"{_WORM}+{_FOOD}")

        union = list(dict.fromkeys([*_declared(_WORM), *_declared(_FOOD)]))
        assert [answer.category for answer in registry.gene_lists()] == union

        # The name of a merge is the +-join of what went in, but it cannot say which
        # component each half came from — and that is exactly what attribution needs, so
        # who contributed comes from the record and not from splitting the name.
        renamed = self._registry(tmp_path / "renamed", name="merged")
        shared = next(category for category in _declared(_WORM) if category in _declared(_FOOD))

        answer = renamed.gene_list(shared, "merged")

        assert [source.component for source in answer.sources] == [
            _WORM_COMPONENT,
            _FOOD_COMPONENT,
        ]


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
        tmp_path.mkdir(parents=True, exist_ok=True)
        source = tmp_path / f"{name}.gtf"
        source.write_text(gtf)
        _register_by_path(tmp_path, source, name)
        return AnnotationRegistry.locate("tiny", tmp_path)

    def test_a_stem_resolves_to_the_versioned_id_or_to_itself_when_unversioned(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path, _gtf_declaring(_ALONE))

        answer = registry.resolve_gene_ids([_ALONE_STEM], "mine")

        assert answer.resolved == {_ALONE_STEM: (_ALONE,)}
        assert answer.unresolved == ()
        assert (answer.assembly, answer.annotation) == ("tiny", "mine")

        # WormBase and SGD never versioned a gene id, and an Ensembl-shaped assumption
        # would leave both unresolvable. An id with no version is its own stem.
        unversioned_dir = tmp_path / "unversioned"
        unversioned_dir.mkdir()
        unversioned = self._registry(unversioned_dir, _GTF)
        unversioned_answer = unversioned.resolve_gene_ids(["g1", "g2"], "mine")
        assert unversioned_answer.resolved == {"g1": ("g1",)}
        assert unversioned_answer.unresolved == ("g2",)

    def test_an_unresolvable_stem_rides_back_and_order_and_repeats_are_handled(
        self, tmp_path: Path
    ) -> None:
        # The `gencode_v50lift37` case: nine stems name two genes each, eight of them a
        # pseudoautosomal pair. Answering with the first would hand back the X copy of a Y
        # gene without ever saying a choice had been made.
        par_stem_registry = self._registry(
            tmp_path / "par-stem", _gtf_declaring(_PAR_Y, _PAR_X, _ALONE)
        )
        par_stem_answer = par_stem_registry.resolve_gene_ids([_PAR_STEM], "mine")
        assert par_stem_answer.resolved[_PAR_STEM] == (_PAR_X, _PAR_Y)

        registry = self._registry(tmp_path, _gtf_declaring(_ALONE))

        answer = registry.resolve_gene_ids([_ALONE_STEM, _ABSENT_STEM], "mine")
        assert answer.unresolved == (_ABSENT_STEM,)
        assert _ABSENT_STEM not in answer.resolved
        # …and what did resolve is still there, so an absence costs the caller nothing else.
        assert answer.gene_ids == [_ALONE]

        empty_answer = registry.resolve_gene_ids([], "mine")
        assert (dict(empty_answer.resolved), empty_answer.unresolved) == ({}, ())

        # A caller passing a few thousand at once reads its own list against the answer.
        par_dir = tmp_path / "par"
        par_dir.mkdir()
        par_registry = self._registry(par_dir, _gtf_declaring(_ALONE, _PAR_X))
        par_answer = par_registry.resolve_gene_ids(
            [_PAR_STEM, _ABSENT_STEM, _ALONE_STEM, _PAR_STEM], "mine"
        )
        assert list(par_answer.resolved) == [_PAR_STEM, _ALONE_STEM]
        assert par_answer.gene_ids == [_PAR_X, _ALONE]
        assert par_answer.unresolved == (_ABSENT_STEM,)

    def test_an_unregistered_name_or_no_name_is_resolved_by_the_same_lookup(
        self, tmp_path: Path
    ) -> None:
        # An exon-level GTF registers as exons alone, since reconstructing the features
        # above them is off by default. Every stem would come back unresolved, which reads
        # as *this annotation has none of your genes* and is a different fact entirely.
        bare_registry = self._registry(tmp_path / "bare", _BARE_GTF)
        with pytest.raises(NoGeneFeaturesError) as no_genes:
            bare_registry.resolve_gene_ids(["g1"], "mine")
        no_genes_message = str(no_genes.value)
        assert "mine" in no_genes_message
        assert "--infer-genes" in no_genes_message  # the argument that rebuilds it with genes

        # Resolved through the same lookup every other question goes through, so the
        # message and its repair are the ones that surface already has.
        registry = AnnotationRegistry.locate(_CURATED_ASSEMBLY, tmp_path)
        with pytest.raises(AnnotationNotRegisteredError) as excinfo:
            registry.resolve_gene_ids([_ALONE_STEM], _CURATED)
        assert f"genome register-annotation {_CURATED_ASSEMBLY} {_CURATED}" in str(excinfo.value)

        default_dir = tmp_path / "default"
        default_dir.mkdir()
        default_registry = self._registry(default_dir, _gtf_declaring(_ALONE))
        assert default_registry.default == "mine"
        assert default_registry.resolve_gene_ids([_ALONE_STEM]).annotation == "mine"

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


class TestAddressedByAssembly:
    """``gene_list`` and ``gene_lists`` addressed by assembly name rather than opened.

    A registry for the length of the call, exactly as ``annotation_status`` is, so there
    is no second code path to keep in step.
    """

    def test_it_answers_for_an_assembly_named_rather_than_opened_or_raises_if_nothing_is(
        self, tmp_path: Path
    ) -> None:
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(tmp_path, source, _CURATED, assembly=_CURATED_ASSEMBLY)
        category = _declared(_CURATED)[0]

        answer = gene_list(_CURATED_ASSEMBLY, category, cache_dir=tmp_path)

        assert answer.annotation == _CURATED
        assert [entry.category for entry in gene_lists(_CURATED_ASSEMBLY, cache_dir=tmp_path)] == (
            list(_declared(_CURATED))
        )

        with pytest.raises(AnnotationNotRegisteredError):
            gene_list(_CURATED_ASSEMBLY, "rRNA", annotation=_CURATED, cache_dir=tmp_path / "empty")


# ---------------------------------------------------------------------------------------
# What registering an annotation answers with, built by hand from a record's own fields
# ---------------------------------------------------------------------------------------

# Nothing below registers anything: these are the values registration *returns*, so they
# are built from the fields a record carries. The one exception writes a real annotation,
# because reading back a record an older version wrote is a claim about a file on disk and
# not about a dataclass.


def _completion(kind: str, name: str, **details: object) -> CompletionRecord:
    """A completion record with everything filled in, so ``as_json`` has every key."""
    return CompletionRecord(
        kind=kind,
        name=name,
        files={f"{name}.db": 34, f"{name}.gtf": 12},
        source_url="https://example.org/x.gz",
        sha256="1a2b3c",
        tool_versions={"samtools": "1.21"},
        package_version="2026.8.0",
        completed_at="2026-08-12T09:00:00+00:00",
        details=dict(details),
    )


def _state_row(**overrides: object) -> AnnotationStatusRow:
    """An offered-but-not-registered status row, with any field overridden by keyword."""
    fields: dict[str, object] = {
        "name": "gencode_v50",
        "offered": True,
        "registered": False,
        "broken": False,
        "default": True,
        "provider": "GENCODE",
        "version": "v50",
        "url": "https://example.org/gencode_v50.gtf.gz",
        "sha256": None,
        "path": None,
        "problem": None,
        "repair": None,
    }
    fields.update(overrides)
    return AnnotationStatusRow(**fields)  # type: ignore[arg-type]


def _state(default: str | None, *rows: AnnotationStatusRow) -> AnnotationStatus:
    """A status report over ``rows``, naming ``default`` as the default annotation."""
    return AnnotationStatus(
        assembly="hg38",
        directory=Path("/data/genome/hg38"),
        default_annotation=default,
        annotations=rows,
    )


class TestTheStateARowIsIn:
    """``AnnotationStatusRow.state`` — the four-way state, in the words a surface prints.

    It was the CLI's to derive, which meant the precedence between ``broken`` and
    ``registered`` was written twice: once as the row's invariant, once as an ``if``
    ordering in a renderer. It is the row's, and the renderer chooses column widths.
    """

    def test_the_four_states_and_brokens_precedence_over_both(self) -> None:
        assert _state_row(registered=True).state == "registered"
        assert _state_row().state == "offered, not registered"
        assert _state_row(offered=False, registered=True, default=False).state == (
            "registered, not offered"
        )
        # Broken beats registered rather than reading as neither: a directory nothing
        # vouches for is not registered, so reporting the absence of a registration would
        # be true and useless. What needs acting on is that it is broken.
        assert _state_row(broken=True).state == "broken"
        assert _state_row(offered=False, broken=True).state == "broken"  # and beats unlisted


class TestTheDefaultAnnotationLine:
    """``AnnotationStatus.default_summary`` — the closing line, and what to do about it.

    Four answers, and the two that name a command take it off an interface: the broken
    one from the row's own ``repair``, the absent one from
    :func:`annotation_register_command`. Neither is concatenated here or in a surface.
    """

    def test_the_four_summary_lines(self) -> None:
        assert _state(None).default_summary == "default: (none)"
        assert (
            _state("gencode_v50", _state_row(registered=True)).default_summary
            == "default: gencode_v50"
        )

        absent = _state("gencode_v50", _state_row())
        assert absent.default_summary == (
            "default: gencode_v50 — not registered here; register it with "
            "`genome register-annotation hg38 gencode_v50`"
        )
        # The command it names is the one the package spells once, not a copy of it.
        assert annotation_register_command("hg38", "gencode_v50") in absent.default_summary

        repair = "genome register-annotation hg38 gencode_v50 --force"
        broken = _state("gencode_v50", _state_row(broken=True, repair=repair))
        assert broken.default_summary == (
            f"default: gencode_v50 — broken here; repair it with `{repair}`"
        )

        # Named by the table, and no row here is about it: there is no row to read
        # `registered` or `broken` off, and the answer is the one that registers it.
        fresh = _state("gencode_v50")
        assert fresh.default_row is None
        assert "not registered here" in fresh.default_summary


class TestTheJsonKeysAndTheirOrder:
    """``as_json`` — every ``--json`` surface here, pinned key for key and in order.

    ``--json`` is what a script parses, so a key renamed, dropped or reordered is a break
    whether or not anything in this suite notices. These assert the whole list rather than
    a key inside it, which is the only form that fails on an addition.
    """

    def test_a_registered_annotation_is_a_record_plus_what_a_record_does_not_hold(self) -> None:
        # The same shape a registered assembly serializes in, deliberately: a record plus
        # the two facts a record does not hold about itself. test_download pins that half.
        annotation = RegisteredAnnotation(
            assembly="hg38",
            directory=Path("/data/genome/hg38/gtf/gencode_v50"),
            record=_completion("annotation", "gencode_v50", chromosomes_checked=True),
        )

        assert list(annotation.as_json()) == [
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
        assert annotation.as_json()["directory"] == "/data/genome/hg38/gtf/gencode_v50"

    def test_a_status_row_and_the_status_around_it_serialize_with_no_derived_keys(self) -> None:
        assert list(_state_row().as_json()) == [
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
        # Writing `state` out would be a second spelling of the precedence for a parser to
        # disagree with, and the three fields it comes from are all here already.
        assert "state" not in _state_row(broken=True).as_json()

        status = _state(
            "gencode_v50", _state_row(), _state_row(name="mine", offered=False, registered=True)
        )
        assert list(status.as_json()) == [
            "assembly",
            "directory",
            "default_annotation",
            "annotations",
        ]
        assert [row["name"] for row in status.as_json()["annotations"]] == ["gencode_v50", "mine"]
        assert status.as_json()["directory"] == "/data/genome/hg38"

    def test_the_whole_report_survives_json_dumps_unchanged(self) -> None:
        # The last step the CLI takes, taken here: nothing in a report is a type json
        # cannot write, which is what makes the `--json` path total.
        status = _state("gencode_v50", _state_row(broken=True, repair="x", problem="y"))

        assert json.loads(json.dumps(status.as_json())) == status.as_json()


class TestReadingBackWhatWasChecked:
    """``chromosome_check_summary`` — one sentence per state, and never the wrong one.

    The states differ in what a reader should do about them, which is why they are told
    apart at all: an annotation registered before its assembly is waiting for the
    assembly, and one whose check the caller stood down is waiting for nothing.
    """

    _ADVICE = "register the assembly first"

    def test_each_check_state_reads_as_its_own_sentence_and_a_registration_reads_off_its_record(
        self,
    ) -> None:
        # Silence is not how a pass is reported: a surface printing nothing about the
        # check reads exactly like one printing that it passed.
        checked = chromosome_check_summary(
            {"chromosomes_checked": True, "chromosomes_unchecked_because": None}
        )
        assert "chromosomes checked" in checked
        assert self._ADVICE not in checked

        no_sizes = chromosome_check_summary(
            {
                "chromosomes_checked": False,
                "chromosomes_unchecked_because": UNCHECKED_NO_CHROM_SIZES,
            }
        )
        assert "chromosomes not checked" in no_sizes
        assert self._ADVICE in no_sizes

        # An override is never told to register the assembly: it may well be registered,
        # and the caller turned the check off on purpose.
        overridden = chromosome_check_summary(
            {
                "chromosomes_checked": False,
                "chromosomes_unchecked_because": UNCHECKED_CALLER_OVERRIDE,
            }
        )
        assert "stood down" in overridden
        assert self._ADVICE not in overridden

        # A registration answers off its own record's details — the surface never spells
        # the two `details` keys itself.
        registered = RegisteredAnnotation(
            assembly="hg38",
            directory=Path("/data/genome/hg38/gtf/gencode_v50"),
            record=_completion(
                "annotation",
                "gencode_v50",
                chromosomes_checked=False,
                chromosomes_unchecked_because=UNCHECKED_CALLER_OVERRIDE,
            ),
        )
        assert registered.chromosome_check == chromosome_check_summary(registered.record.details)
        assert "stood down" in registered.chromosome_check

    def test_every_known_state_is_distinct_and_an_unknown_reason_reads_as_unknown(self) -> None:
        summaries = {
            chromosome_check_summary(details)
            for details in (
                {"chromosomes_checked": True, "chromosomes_unchecked_because": None},
                {
                    "chromosomes_checked": False,
                    "chromosomes_unchecked_because": UNCHECKED_NO_CHROM_SIZES,
                },
                {
                    "chromosomes_checked": False,
                    "chromosomes_unchecked_because": UNCHECKED_CALLER_OVERRIDE,
                },
                {"chromosomes_checked": False},
            )
        }
        assert len(summaries) == 4

        # Forward as well as backward: a record from a later version claiming some third
        # reason is one this version cannot report, which is the same as not knowing.
        future = chromosome_check_summary(
            {"chromosomes_checked": False, "chromosomes_unchecked_because": "some-later-reason"}
        )
        assert "does not say why" in future

    def test_a_record_written_before_the_reason_existed_reads_as_unknown(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        # The real back-compatibility case, on a record that is on disk: an older version
        # wrote the bare bool, and which of the two reasons it stood for is not knowable.
        # It must read as neither, and reading it must not raise.
        _register_by_path(tmp_path, data_dir / "tiny.gtf", _NAME)
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


class TestARecordIsCarriedWholeRatherThanCopiedOut:
    """The properties that exist so a surface never re-reads a directory."""

    def test_a_registered_annotation_answers_off_the_record_it_carries(self) -> None:
        annotation = RegisteredAnnotation(
            assembly="hg38",
            directory=Path("/data/genome/hg38/gtf/gencode_v50"),
            record=_completion("annotation", "gencode_v50"),
        )

        assert annotation.name == "gencode_v50"
        assert annotation.source_url == "https://example.org/x.gz"
        assert annotation.sha256 == "1a2b3c"
        assert annotation.file_names == ["gencode_v50.db", "gencode_v50.gtf"]
        first = annotation.file_names
        first.append("intruder")
        assert annotation.file_names == ["gencode_v50.db", "gencode_v50.gtf"]  # fresh each call


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


def test_registering_an_annotation_imports_nothing_from_the_tf_context() -> None:
    # The same guard for the other edge this module carried, and it cost more than it
    # looked: `genome/tf/__init__.py` imports the motif link table, which imports both
    # gene-keyed halves and the whole motif tree down to the scan, its worker pool and its
    # Parquet sink. Two import lines here meant registering an annotation loaded all of it.
    # A prefix rather than a fixed set, since every module under `genome.tf` is out of
    # bounds and a new one must not be able to arrive unnoticed.
    #
    # What it costs to keep: which species selects a shipped table, and what a row of one
    # says, belong to whoever ships the table. This module answers *which gene ids does
    # this stem name here* — `resolve_gene_ids`, which knows nothing about what it is
    # handed a list of — and a topic that wants an annotation's own gene ids crosses that
    # from its own directory, as `genome.tf.gene`, `genome.tf.cofactor` and
    # `genome.homology.annotation` all do.
    reached = {name for name in _module_level_imports(gtf_module) if name.startswith("genome.tf")}

    assert reached == set()
