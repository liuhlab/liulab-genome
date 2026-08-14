"""Tests for genome.io.completion — the record a finished build writes.

Everything here works on plain directories of small files: the record is about presence,
sizes and provenance, so nothing needs a genome, a network or a native tool. The
round-trip through write and read is checked as a property over generated records.
"""

from __future__ import annotations

import json
import tempfile
from datetime import datetime
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome import __version__
from genome.io.completion import (
    RECORD_NAME,
    WORK_DIR_NAME,
    CompletionRecord,
    RegistrationMismatchError,
    UnfinishedRegistrationError,
    build_record,
    check_registration,
    clear_work_dir,
    disagreements,
    read_record,
    record_path,
    tool_versions,
    work_dir,
    write_record,
)

# JSON-representable leaves, so a generated record survives a round trip through JSON.
_LEAVES = st.none() | st.booleans() | st.integers() | st.text()

_RECORDS = st.builds(
    CompletionRecord,
    kind=st.text(),
    name=st.text(),
    files=st.dictionaries(st.text(), st.integers(min_value=0)),
    source_url=st.none() | st.text(),
    sha256=st.none() | st.text(),
    tool_versions=st.dictionaries(st.text(), st.text()),
    package_version=st.text(),
    completed_at=st.text(),
    details=st.dictionaries(st.text(), _LEAVES | st.lists(_LEAVES)),
)


def _build(directory: Path, *names: str, size: int = 3) -> list[Path]:
    """Create ``names`` under ``directory`` (parents included) and return their paths."""
    paths = []
    for name in names:
        path = directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x" * size)
        paths.append(path)
    return paths


# --- the record round-trips --------------------------------------------------


@given(record=_RECORDS)
def test_a_record_survives_write_and_read_unchanged(record: CompletionRecord) -> None:
    # A fresh directory per example: @given may not lean on a function-scoped fixture.
    with tempfile.TemporaryDirectory() as name:
        directory = Path(name) / "roundtrip"

        write_record(directory, record)

        assert read_record(directory) == record


def test_the_record_is_written_under_one_well_known_name(tmp_path: Path) -> None:
    (fasta,) = _build(tmp_path, "tiny.fa")

    written = write_record(
        tmp_path, build_record(tmp_path, kind="genome", name="tiny", files=[fasta])
    )

    assert written == record_path(tmp_path) == tmp_path / RECORD_NAME
    assert written.is_file()


def test_writing_leaves_no_temporary_file_behind(tmp_path: Path) -> None:
    (fasta,) = _build(tmp_path, "tiny.fa")
    record = build_record(tmp_path, kind="genome", name="tiny", files=[fasta])

    write_record(tmp_path, record)
    write_record(tmp_path, record)  # a second run replaces the record in place

    assert sorted(p.name for p in tmp_path.iterdir()) == [RECORD_NAME, "tiny.fa"]


def test_a_reader_never_sees_a_half_written_record(tmp_path: Path) -> None:
    # The record is renamed over its destination, so whatever is at that path is
    # always a complete document — including while a rewrite is in flight.
    (fasta,) = _build(tmp_path, "tiny.fa")
    first = build_record(tmp_path, kind="genome", name="tiny", files=[fasta], source_url="a://one")
    write_record(tmp_path, first)
    second = build_record(tmp_path, kind="genome", name="tiny", files=[fasta], source_url="a://two")

    write_record(tmp_path, second)

    payload = json.loads((tmp_path / RECORD_NAME).read_text())
    assert payload["source_url"] == "a://two"


def test_no_record_reads_as_unfinished(tmp_path: Path) -> None:
    assert read_record(tmp_path) is None


@pytest.mark.parametrize("content", ["", "not json at all", "[]", '{"kind": "genome"}'])
def test_an_unusable_record_reads_as_unfinished(tmp_path: Path, content: str) -> None:
    # Truncated, malformed or the wrong shape: all of them mean "not finished" here.
    # Which of them should raise instead is the repairing caller's decision, not this one's.
    (tmp_path / RECORD_NAME).write_text(content)

    assert read_record(tmp_path) is None


# --- what a record carries ---------------------------------------------------


def test_a_record_carries_the_provenance_of_the_build(tmp_path: Path) -> None:
    fasta, twobit = _build(tmp_path, "tiny.fa", "tiny.2bit", size=7)

    record = build_record(
        tmp_path,
        kind="genome",
        name="tiny",
        files=[fasta, twobit],
        source_url="https://example.org/tiny.fa.gz",
        sha256="1a2b3c",
        details={"note": "kind-specific"},
    )

    assert record.kind == "genome"
    assert record.name == "tiny"
    assert record.files == {"tiny.fa": 7, "tiny.2bit": 7}
    assert record.source_url == "https://example.org/tiny.fa.gz"
    assert record.sha256 == "1a2b3c"
    assert record.package_version == __version__
    assert record.details == {"note": "kind-specific"}
    # ISO-8601 in UTC, so the moment is unambiguous months later.
    assert datetime.fromisoformat(record.completed_at).utcoffset() is not None


def test_claimed_paths_are_relative_so_the_directory_stays_movable(tmp_path: Path) -> None:
    built = tmp_path / "before"
    fasta, index = _build(built, "tiny.fa", "index/star/SA")
    write_record(built, build_record(built, kind="genome", name="tiny", files=[fasta, index]))

    moved = tmp_path / "after"
    built.rename(moved)

    record = read_record(moved)
    assert record is not None
    assert sorted(record.files) == ["index/star/SA", "tiny.fa"]
    assert disagreements(moved, record) == []


def test_recording_a_file_outside_the_directory_is_refused(tmp_path: Path) -> None:
    (stray,) = _build(tmp_path, "stray.fa")
    directory = tmp_path / "build"
    directory.mkdir()

    with pytest.raises(ValueError, match="only files inside its own directory"):
        build_record(directory, kind="genome", name="tiny", files=[stray])


def test_recording_a_file_that_is_not_there_yet_is_refused(tmp_path: Path) -> None:
    # The record is written last; claiming a file that does not exist is a build bug.
    with pytest.raises(FileNotFoundError, match="written last"):
        build_record(tmp_path, kind="genome", name="tiny", files=[tmp_path / "tiny.fa"])


def test_a_tool_that_cannot_be_run_is_left_out_rather_than_raising(tmp_path: Path) -> None:
    (fasta,) = _build(tmp_path, "tiny.fa")

    record = build_record(
        tmp_path, kind="genome", name="tiny", files=[fasta], tools=["definitelyNotInstalled"]
    )

    assert record.tool_versions == {}


def test_tool_versions_reports_a_tool_that_answers() -> None:
    versions = tool_versions(["python", "definitelyNotInstalled"])

    assert "definitelyNotInstalled" not in versions
    assert versions["python"].startswith("Python")


def test_a_tool_that_will_not_identify_itself_is_left_out_too(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # Installed, and rejecting `--version` the way several UCSC binaries do. An absent
    # key means *unknown*, and it must mean that for both reasons a version can be
    # unknown — otherwise a record would carry an empty string as if it were a fact.
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    declines = bin_dir / "faToTwoBit"
    declines.write_text("#!/bin/sh\necho '--version is not a valid option' >&2\nexit 255\n")
    declines.chmod(0o755)
    monkeypatch.setenv("PATH", str(bin_dir))

    assert tool_versions(["faToTwoBit"]) == {}


# --- holding a directory to its record ---------------------------------------


def test_a_directory_that_matches_its_record_has_nothing_to_report(tmp_path: Path) -> None:
    fasta, twobit = _build(tmp_path, "tiny.fa", "tiny.2bit")
    record = build_record(tmp_path, kind="genome", name="tiny", files=[fasta, twobit])

    assert disagreements(tmp_path, record) == []


def test_a_deleted_file_is_reported_by_name(tmp_path: Path) -> None:
    fasta, twobit = _build(tmp_path, "tiny.fa", "tiny.2bit")
    record = build_record(tmp_path, kind="genome", name="tiny", files=[fasta, twobit])

    twobit.unlink()

    (bad,) = disagreements(tmp_path, record)
    assert bad.path == "tiny.2bit"
    assert bad.actual is None
    assert str(bad) == "tiny.2bit: recorded 3 bytes, missing"


def test_a_truncated_file_is_reported_with_both_sizes(tmp_path: Path) -> None:
    (fasta,) = _build(tmp_path, "tiny.fa", size=100)
    record = build_record(tmp_path, kind="genome", name="tiny", files=[fasta])

    fasta.write_text("x")

    (bad,) = disagreements(tmp_path, record)
    assert (bad.path, bad.expected, bad.actual) == ("tiny.fa", 100, 1)
    assert str(bad) == "tiny.fa: recorded 100 bytes, found 1"


def test_every_offender_is_reported_not_just_the_first(tmp_path: Path) -> None:
    fasta, twobit, sizes = _build(tmp_path, "tiny.fa", "tiny.2bit", "tiny.chrom.sizes")
    record = build_record(tmp_path, kind="genome", name="tiny", files=[fasta, twobit, sizes])

    twobit.unlink()
    sizes.write_text("longer than before")

    assert [bad.path for bad in disagreements(tmp_path, record)] == [
        "tiny.2bit",
        "tiny.chrom.sizes",
    ]


def test_a_same_size_edit_goes_unnoticed_because_contents_are_never_read(tmp_path: Path) -> None:
    # The check is presence and size only, deliberately: reading contents would make
    # reopening a prepared human genome a multi-gigabyte read rather than milliseconds.
    (fasta,) = _build(tmp_path, "tiny.fa", size=100)
    record = build_record(tmp_path, kind="genome", name="tiny", files=[fasta])

    fasta.write_text("y" * 100)

    assert disagreements(tmp_path, record) == []


# --- finished, fresh, or broken ----------------------------------------------

_REPAIR = "genome register tiny --force"


def test_a_finished_build_answers_with_its_record(tmp_path: Path) -> None:
    (fasta,) = _build(tmp_path, "tiny.fa")
    written = build_record(tmp_path, kind="genome", name="tiny", files=[fasta])
    write_record(tmp_path, written)

    assert check_registration(tmp_path, repair=_REPAIR) == written


def test_an_absent_directory_is_fresh_rather_than_broken(tmp_path: Path) -> None:
    assert check_registration(tmp_path / "never-built", repair=_REPAIR) is None


def test_an_empty_directory_is_fresh_rather_than_broken(tmp_path: Path) -> None:
    assert check_registration(tmp_path, repair=_REPAIR) is None


def test_a_directory_holding_only_a_download_is_still_fresh(tmp_path: Path) -> None:
    # The working area is working state, not a claimed output, so an interrupted
    # download does not make a directory that was never registered look broken.
    _build(work_dir(tmp_path), "tiny.fa.gz")

    assert check_registration(tmp_path, repair=_REPAIR) is None


def test_files_with_no_record_raise_and_name_the_repair(tmp_path: Path) -> None:
    _build(tmp_path, "tiny.fa", "tiny.2bit")

    with pytest.raises(UnfinishedRegistrationError) as excinfo:
        check_registration(tmp_path, repair=_REPAIR)

    message = str(excinfo.value)
    assert "tiny.fa" in message
    assert _REPAIR in message


def test_a_record_that_disagrees_raises_naming_which_file_and_how(tmp_path: Path) -> None:
    fasta, twobit = _build(tmp_path, "tiny.fa", "tiny.2bit", size=100)
    write_record(
        tmp_path, build_record(tmp_path, kind="genome", name="tiny", files=[fasta, twobit])
    )
    twobit.write_text("x")  # truncated behind our back

    with pytest.raises(RegistrationMismatchError) as excinfo:
        check_registration(tmp_path, repair=_REPAIR)

    message = str(excinfo.value)
    assert "tiny.2bit: recorded 100 bytes, found 1" in message
    assert "tiny.fa:" not in message  # the file that still agrees is not accused
    assert _REPAIR in message


# --- the working area --------------------------------------------------------


def test_the_working_area_is_a_hidden_directory_inside_the_build(tmp_path: Path) -> None:
    assert work_dir(tmp_path) == tmp_path / WORK_DIR_NAME
    assert WORK_DIR_NAME.startswith(".")


def test_clearing_the_working_area_removes_everything_in_it(tmp_path: Path) -> None:
    _build(work_dir(tmp_path), "tiny.fa.gz")
    (kept,) = _build(tmp_path, "tiny.fa")

    clear_work_dir(tmp_path)

    assert not work_dir(tmp_path).exists()
    assert kept.is_file()


def test_clearing_a_working_area_that_was_never_made_is_fine(tmp_path: Path) -> None:
    clear_work_dir(tmp_path)  # no raise


def test_an_annotation_registered_first_does_not_make_its_assembly_look_broken(
    tmp_path: Path,
) -> None:
    # An Assembly dir hosts the gtf/ and index/ subtrees other contexts own, and each
    # carries its own record. Registering an annotation before its assembly is a
    # documented flow, so the assembly must still read as a fresh registration rather
    # than as a run that was interrupted.
    (tmp_path / "gtf" / "gencode_v50").mkdir(parents=True)
    (tmp_path / "index" / "star_gencode_v50").mkdir(parents=True)

    assert check_registration(tmp_path, repair="...", ignore={"gtf", "index"}) is None


def test_a_foreign_subtree_does_not_hide_a_real_interrupted_run(tmp_path: Path) -> None:
    # Ignoring those subtrees must not ignore the assembly's own leftovers beside them.
    (tmp_path / "gtf" / "gencode_v50").mkdir(parents=True)
    (tmp_path / "hg38.fa").write_text(">chr1\nACGT\n")

    with pytest.raises(UnfinishedRegistrationError, match=r"hg38\.fa"):
        check_registration(tmp_path, repair="genome register hg38 --force", ignore={"gtf", "index"})
