"""Tests for genome.store.completion — the record a finished build writes.

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
from genome.store.completion import (
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

from ..conftest import StubBinary

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


def test_writing_is_atomic_under_one_well_known_name_and_leaves_no_temp_file(
    tmp_path: Path,
) -> None:
    (fasta,) = _build(tmp_path, "tiny.fa")
    first = build_record(tmp_path, kind="genome", name="tiny", files=[fasta], source_url="a://one")

    written = write_record(tmp_path, first)

    assert written == record_path(tmp_path) == tmp_path / RECORD_NAME
    assert written.is_file()

    # A second run replaces the record in place, leaving no other file behind.
    second = build_record(tmp_path, kind="genome", name="tiny", files=[fasta], source_url="a://two")
    write_record(tmp_path, second)
    assert sorted(p.name for p in tmp_path.iterdir()) == [RECORD_NAME, "tiny.fa"]

    # The record is renamed over its destination, so whatever is at that path is
    # always a complete document — including while a rewrite is in flight.
    payload = json.loads((tmp_path / RECORD_NAME).read_text())
    assert payload["source_url"] == "a://two"


def test_no_record_or_an_unusable_one_reads_as_unfinished(tmp_path: Path) -> None:
    assert read_record(tmp_path) is None

    # Truncated, malformed or the wrong shape: all of them mean "not finished" here.
    # Which of them should raise instead is the repairing caller's decision, not this
    # one's.
    for content in ("", "not json at all", "[]", '{"kind": "genome"}'):
        (tmp_path / RECORD_NAME).write_text(content)
        assert read_record(tmp_path) is None, f"content {content!r} should read as unfinished"


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


def test_recording_a_file_outside_or_not_yet_in_the_directory_is_refused(tmp_path: Path) -> None:
    (stray,) = _build(tmp_path, "stray.fa")
    directory = tmp_path / "build"
    directory.mkdir()

    with pytest.raises(ValueError, match="only files inside its own directory"):
        build_record(directory, kind="genome", name="tiny", files=[stray])

    # The record is written last; claiming a file that does not exist is a build bug.
    with pytest.raises(FileNotFoundError, match="written last"):
        build_record(tmp_path, kind="genome", name="tiny", files=[tmp_path / "tiny.fa"])


def test_tool_versions_reports_what_answers_and_omits_what_cannot(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, stub_binary: StubBinary
) -> None:
    versions = tool_versions(["python", "definitelyNotInstalled"])
    assert "definitelyNotInstalled" not in versions
    assert versions["python"].startswith("Python")

    (fasta,) = _build(tmp_path, "tiny.fa")
    record = build_record(
        tmp_path, kind="genome", name="tiny", files=[fasta], tools=["definitelyNotInstalled"]
    )
    assert record.tool_versions == {}

    # Installed, and rejecting `--version` the way several UCSC binaries do. An absent
    # key means *unknown*, and it must mean that for both reasons a version can be
    # unknown — otherwise a record would carry an empty string as if it were a fact.
    bin_dir = tmp_path / "bin"
    stub_binary(bin_dir, "faToTwoBit", "echo '--version is not a valid option' >&2\nexit 255")
    monkeypatch.setenv("PATH", str(bin_dir))
    assert tool_versions(["faToTwoBit"]) == {}


# --- holding a directory to its record ---------------------------------------


def test_a_matching_directory_or_a_same_size_edit_has_nothing_to_report(tmp_path: Path) -> None:
    (fasta,) = _build(tmp_path, "tiny.fa", size=100)
    record = build_record(tmp_path, kind="genome", name="tiny", files=[fasta])

    assert disagreements(tmp_path, record) == []

    # The check is presence and size only, deliberately: reading contents would make
    # reopening a prepared human genome a multi-gigabyte read rather than milliseconds.
    fasta.write_text("y" * 100)
    assert disagreements(tmp_path, record) == []


def test_a_deleted_or_a_truncated_file_is_reported_with_the_right_message(tmp_path: Path) -> None:
    fasta, twobit = _build(tmp_path, "tiny.fa", "tiny.2bit")
    record = build_record(tmp_path, kind="genome", name="tiny", files=[fasta, twobit])
    twobit.unlink()

    (deleted,) = disagreements(tmp_path, record)
    assert deleted.path == "tiny.2bit"
    assert deleted.actual is None
    assert str(deleted) == "tiny.2bit: recorded 3 bytes, missing"

    (fasta_alone,) = _build(tmp_path, "tiny.fa", size=100)
    truncated_record = build_record(tmp_path, kind="genome", name="tiny", files=[fasta_alone])
    fasta_alone.write_text("x")

    (truncated,) = disagreements(tmp_path, truncated_record)
    assert (truncated.path, truncated.expected, truncated.actual) == ("tiny.fa", 100, 1)
    assert str(truncated) == "tiny.fa: recorded 100 bytes, found 1"


def test_every_offender_is_reported_not_just_the_first(tmp_path: Path) -> None:
    fasta, twobit, sizes = _build(tmp_path, "tiny.fa", "tiny.2bit", "tiny.chrom.sizes")
    record = build_record(tmp_path, kind="genome", name="tiny", files=[fasta, twobit, sizes])

    twobit.unlink()
    sizes.write_text("longer than before")

    assert [bad.path for bad in disagreements(tmp_path, record)] == [
        "tiny.2bit",
        "tiny.chrom.sizes",
    ]


# --- finished, fresh, or broken ----------------------------------------------

_REPAIR = "genome assembly register tiny --force"


def test_a_finished_build_answers_with_its_record(tmp_path: Path) -> None:
    (fasta,) = _build(tmp_path, "tiny.fa")
    written = build_record(tmp_path, kind="genome", name="tiny", files=[fasta])
    write_record(tmp_path, written)

    assert check_registration(tmp_path, repair=_REPAIR) == written


def test_an_absent_or_empty_directory_or_one_holding_only_a_download_is_fresh(
    tmp_path: Path,
) -> None:
    assert check_registration(tmp_path / "never-built", repair=_REPAIR) is None
    assert check_registration(tmp_path, repair=_REPAIR) is None

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


def test_the_working_area_is_hidden_and_clearing_it_removes_only_its_own_files(
    tmp_path: Path,
) -> None:
    assert work_dir(tmp_path) == tmp_path / WORK_DIR_NAME
    assert WORK_DIR_NAME.startswith(".")

    clear_work_dir(tmp_path)  # never made yet — no raise

    _build(work_dir(tmp_path), "tiny.fa.gz")
    (kept,) = _build(tmp_path, "tiny.fa")

    clear_work_dir(tmp_path)

    assert not work_dir(tmp_path).exists()
    assert kept.is_file()


def test_an_ignored_subtree_hides_its_own_leftovers_but_not_a_real_interrupted_run(
    tmp_path: Path,
) -> None:
    # An Assembly dir hosts the gtf/ and index/ subtrees other contexts own, and each
    # carries its own record. Registering an annotation before its assembly is a
    # documented flow, so the assembly must still read as a fresh registration rather
    # than as a run that was interrupted.
    (tmp_path / "gtf" / "gencode_v50").mkdir(parents=True)
    (tmp_path / "index" / "star_gencode_v50").mkdir(parents=True)
    assert check_registration(tmp_path, repair="...", ignore={"gtf", "index"}) is None

    # But ignoring those subtrees must not ignore the assembly's own leftovers beside
    # them.
    (tmp_path / "hg38.fa").write_text(">chr1\nACGT\n")
    with pytest.raises(UnfinishedRegistrationError, match=r"hg38\.fa"):
        check_registration(
            tmp_path, repair="genome assembly register hg38 --force", ignore={"gtf", "index"}
        )
