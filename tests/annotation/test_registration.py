"""Tests for genome.annotation.registration — placing a GTF, and the record that ends it.

Both ways in, the chromosome check that guards them, the merged annotation a chimera build
derives, and what registering answers with. The **Annotation database** has no test module
of its own: it is exercised through the real builds these registrations perform, and
through the helper that reads a built database's feature types back.

What a broken directory tells the caller to run is asserted here and over in
``test_registry``, which owns the four-way state it is read out of.
"""

from __future__ import annotations

import ast
import gzip
import json
from dataclasses import asdict
from pathlib import Path

import gffutils
import pytest
from hypothesis import given
from hypothesis import strategies as st

import genome
from genome.annotation import curated as curated_module
from genome.annotation import database as database_module
from genome.annotation import metadata as metadata_module
from genome.annotation import registration as registration_module
from genome.annotation import registry as registry_module
from genome.annotation import stems as stems_module
from genome.annotation.registration import (
    UNCHECKED_CALLER_OVERRIDE,
    UNCHECKED_NO_CHROM_SIZES,
    ChromosomeMismatchError,
    MergeSource,
    RegisteredAnnotation,
    _reject_unknown_chromosomes,
    annotation_dir,
    chromosome_check_summary,
    discard_merged_annotation,
    register_annotation,
    register_gtf,
    register_merged_gtf,
)
from genome.annotation.registry import list_annotations
from genome.store.checksum import ChecksumMismatchError
from genome.store.completion import (
    CompletionRecord,
    RegistrationMismatchError,
    UnfinishedRegistrationError,
    read_record,
    record_path,
    work_dir,
)

from ..assembly.test_source import _module_level_imports
from ..conftest import FakeFetch
from .conftest import (
    _BARE_GTF,
    _GTF,
    _NAME,
    _PINNED_URL,
    _TINY_GTF_SHA256,
    _register_by_name,
    _register_by_path,
    _row,
    _write_chrom_sizes,
)


class TestRegisterByPath:
    """``AnnotationRegistry.register_path`` — the way in for a GTF the table does not list.

    Placing, decompressing and rebuilding, asserted on the paths the registry answers
    with. What a broken directory tells the caller to run is asserted once, over in
    ``test_registry``, where the assembly name that composes it lives.
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
        assert f"genome annotation register tiny {_NAME} --force" in str(excinfo.value)

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


def _feature_types(database_path: Path) -> list[str]:
    """The kinds of feature a built database holds, with the connection closed behind us."""
    database = gffutils.FeatureDB(str(database_path))
    try:
        return sorted(database.featuretypes())
    finally:
        database.conn.close()


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
        assert f"genome annotation register-gtf tiny {source} WS298 --force" in str(excinfo.value)

        payload = register_gtf("tiny", source, "WS298", cache_dir=tmp_path, force=True)
        assert payload.name == "WS298"
        assert list(list_annotations(tmp_path)) == ["WS298"]


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


class TestTheJsonKeysAndTheirOrder:
    """``as_json`` — what registering one annotation serializes as, key for key.

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
# The edges this package must not grow back
# ---------------------------------------------------------------------------------------

#: Every module of the annotation package. Each guard below is a claim about the package
#: rather than about one file, and checking one module would leave the others as the place
#: an edge could arrive unnoticed.
_ANNOTATION_MODULES = (
    registration_module,
    registry_module,
    stems_module,
    database_module,
    curated_module,
    metadata_module,
)


def _files_importing(package_root: Path, library: str) -> list[str]:
    """Every file under ``package_root`` whose source imports ``library``, at any depth.

    Read rather than observed, and *not* restricted to module-level statements: a deferred
    import is still an import of the library, and hiding one inside a function is exactly
    how a dependency that is meant to sit behind one adapter gets a second entrance.
    """
    found: list[str] = []
    for path in sorted(package_root.rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text())):
            names = []
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module is not None:
                names = [node.module]
            if any(name == library or name.startswith(f"{library}.") for name in names):
                found.append(path.relative_to(package_root.parent).as_posix())
                break
    return found


def test_registering_an_annotation_imports_nothing_that_downloads_an_assembly() -> None:
    # The cycle, asserted closed. `assembly.download` imports `assembly.source` at the top
    # of the file, `assembly.chimera_build` imports this package's placement half, and that
    # half used to import the downloader back — once for the annotations subdirectory name,
    # which `assembly.registration` defines and the downloader merely re-exports, and once
    # for the package's one fetch step, which `store.fetch` now holds. Each was a line, and
    # each grows back the moment somebody reaches for a name that happens to be importable
    # from the downloader. Were the edge to return, asking what a chimera is made of would
    # drag the whole annotation build stack in behind it, and the downloader's
    # module-level import of the resolution would go back behind a deferred one.
    forbidden = {
        "genome.assembly.download",
        "genome.assembly.chimera_build",
        "genome.assembly.genome",
    }

    for module in _ANNOTATION_MODULES:
        assert _module_level_imports(module) & forbidden == set()


def test_entering_this_package_first_does_not_run_the_open_a_genome_stack() -> None:
    # The half of the cycle above that lives in the other package's `__init__`, and it
    # fails loudly rather than subtly: placement imports `genome.assembly.registration` for
    # the **Assembly dir**, which runs `genome/assembly/__init__.py` first. Were that file
    # to re-export `Genome` or the chimera build — the two assembly modules that import
    # *this* package — then `import genome.annotation` would come back round to a package
    # it is halfway through importing, and raise. `Genome` is exported from `genome`
    # itself, which is where a caller holds it anyway.
    package_root = Path(genome.__file__).parent
    init = package_root / "assembly" / "__init__.py"
    reached = {
        node.module
        for node in ast.parse(init.read_text()).body
        if isinstance(node, ast.ImportFrom) and node.module is not None
    }

    assert "genome.assembly.genome" not in reached
    assert "genome.assembly.chimera_build" not in reached


def test_the_annotation_fetch_is_the_packages_one_fetch_step() -> None:
    # The positive half of the guard above: the edge is gone because the fetch moved to a
    # module of its own, not because this package started spelling a download itself.
    # Placement is where it lands, since placement is the only half of this package that
    # fetches anything at all.
    assert "genome.store.fetch" in _module_level_imports(registration_module)
    for module in _ANNOTATION_MODULES:
        if module is not registration_module:
            assert "genome.store.fetch" not in _module_level_imports(module)


def test_registering_an_annotation_imports_nothing_from_the_tf_context() -> None:
    # The same guard for the other edge this package carried, and it cost more than it
    # looked: `genome/tf/__init__.py` imports the motif link table, which imports both
    # gene-keyed halves and the whole motif tree down to the scan, its worker pool and its
    # Parquet sink. Two import lines here meant registering an annotation loaded all of it.
    # A prefix rather than a fixed set, since every module under `genome.tf` is out of
    # bounds and a new one must not be able to arrive unnoticed.
    #
    # What it costs to keep: which species selects a shipped table, and what a row of one
    # says, belong to whoever ships the table. This package answers *which gene ids does
    # this stem name here* — `resolve_gene_ids`, which knows nothing about what it is
    # handed a list of — and a topic that wants an annotation's own gene ids crosses that
    # from its own directory, as `genome.tf.gene`, `genome.tf.cofactor` and
    # `genome.homology.annotation` all do.
    for module in _ANNOTATION_MODULES:
        reached = {name for name in _module_level_imports(module) if name.startswith("genome.tf")}
        assert reached == set()


def test_the_gffutils_dependency_has_exactly_one_entrance() -> None:
    # What makes the Annotation database an adapter rather than something embedded: two
    # calls into the library are all this package makes — the build, and the read the stem
    # pass walks — and both live in `annotation.database`. A second importer anywhere in
    # `src/` and the API surface starts to say "gffutils", and what backs an Annotation
    # database stops being one file's business.
    package_root = Path(genome.__file__).parent

    assert _files_importing(package_root, "gffutils") == ["genome/annotation/database.py"]
