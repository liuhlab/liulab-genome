"""Tests for :func:`genome.assembly.assembly_status` — the table set against this disk.

Two questions in one answer, and the second of them is the one nothing in the package
could ask before: *which assemblies are prepared on this machine*. So the cases here are
mostly about what the assembly tree is allowed to hold — a registration, a directory
nothing vouches for, a name no shipped row lists, and entries that are not assemblies at
all.

The shipped table answers throughout rather than one stood up for the test: reporting it
is half of what this call is for.
"""

from __future__ import annotations

from pathlib import Path

from genome.assembly import (
    AssemblyStatus,
    AssemblyStatusRow,
    assembly_status,
    assembly_table,
    is_prepared,
)
from genome.store.completion import build_record, write_record

#: Every name the shipped table lists, in table order — what a fresh machine reports.
_OFFERED = tuple(row.assembly_name for row in assembly_table())


def _register(root: Path, assembly: str) -> Path:
    """Write a record in ``root/assembly`` claiming that assembly finished there."""
    directory = root / assembly
    directory.mkdir(parents=True, exist_ok=True)
    fasta = directory / f"{assembly}.fa"
    fasta.write_text(">chrI\nACGT\n")
    write_record(directory, build_record(directory, kind="genome", name=assembly, files=[fasta]))
    return directory


def _rows(status: AssemblyStatus) -> dict[str, AssemblyStatusRow]:
    """The report's rows keyed by name, for a test asking after one of them."""
    return {row.assembly_name: row for row in status.assemblies}


class TestNothingPreparedHere:
    """A fresh install: the table is the whole answer, and nothing is read off the disk."""

    def test_every_offered_assembly_is_listed_and_nothing_is_created_to_answer(
        self, liulab_data: Path
    ) -> None:
        status = assembly_status()

        assert status.directory == liulab_data / "genome"
        assert tuple(row.assembly_name for row in status.assemblies) == _OFFERED
        assert all(row.offered for row in status.assemblies)
        assert not any(row.registered or row.present for row in status.assemblies)
        assert {row.state for row in status.assemblies} == {"offered, not registered"}
        assert status.registered == ()
        # The question is answered without preparing anything — the tree is not even there.
        assert not status.directory.exists()

    def test_the_summary_names_the_command_that_registers_one(self) -> None:
        status = assembly_status()

        assert "genome assembly register <name>" in status.summary
        # Nothing is registered, so there is no re-check to point at and no directory
        # nothing vouches for to explain.
        assert "verify" not in status.summary
        assert status.unregistered_note is None


class TestWhatThisMachineHolds:
    """The half nothing in the package could ask before: what is prepared here."""

    def test_a_registered_assembly_reads_as_registered_and_carries_its_directory(
        self, liulab_data: Path
    ) -> None:
        directory = _register(liulab_data / "genome", "sacCer3")

        status = assembly_status()
        row = _rows(status)["sacCer3"]

        assert (row.offered, row.registered, row.present) == (True, True, True)
        assert row.state == "registered"
        assert row.directory == str(directory)
        assert status.registered == ("sacCer3",)
        assert is_prepared("sacCer3")
        # A registration is what the summary points a doubter at re-checking.
        assert "genome assembly verify <name>" in status.summary

    def test_a_name_no_row_lists_is_listed_after_the_offered_ones_with_the_columns_blank(
        self, liulab_data: Path
    ) -> None:
        _register(liulab_data / "genome", "danRer11")

        status = assembly_status()
        rows = _rows(status)

        assert tuple(rows)[: len(_OFFERED)] == _OFFERED
        assert tuple(rows)[-1] == "danRer11"
        unlisted = rows["danRer11"]
        assert (unlisted.offered, unlisted.registered) == (False, True)
        assert unlisted.state == "registered, not offered"
        assert (unlisted.species, unlisted.ncbi_name, unlisted.source_url, unlisted.sha256) == (
            None,
            None,
            None,
            None,
        )

    def test_unlisted_names_are_reported_in_sorted_order(self, liulab_data: Path) -> None:
        for name in ("zebrafish", "danRer11", "myref"):
            _register(liulab_data / "genome", name)

        listed = [row.assembly_name for row in assembly_status().assemblies]

        assert listed[len(_OFFERED) :] == ["danRer11", "myref", "zebrafish"]


class TestADirectoryNothingVouchesFor:
    """The third case enumeration meets and asking after one name never did.

    A directory with no record is neither a registration nor absent. Reporting it as
    absent is a lie to someone looking at a full disk; reporting it as an assembly is a
    lie about what is trustworthy. It is reported as being here, unregistered — the same
    record-alone rule, applied to a name that was found rather than asked for.
    """

    def test_an_unlisted_directory_with_no_record_is_here_but_not_registered(
        self, liulab_data: Path
    ) -> None:
        directory = liulab_data / "genome" / "half_built"
        directory.mkdir(parents=True)
        (directory / "half_built.fa").write_text(">chrI\nACGT\n")

        row = _rows(assembly_status())["half_built"]

        assert (row.offered, row.registered, row.present) == (False, False, True)
        assert row.state == "here, not registered"
        assert row.directory == str(directory)

    def test_an_offered_assembly_whose_directory_is_here_is_not_reported_as_absent(
        self, liulab_data: Path
    ) -> None:
        # The state that would otherwise read `offered, not registered` — the same words
        # as an assembly nobody has ever fetched, on a directory that is right there.
        (liulab_data / "genome" / "hg38" / "gtf" / "gencode_v50").mkdir(parents=True)

        row = _rows(assembly_status())["hg38"]

        assert (row.offered, row.registered, row.present) == (True, False, True)
        assert row.state == "here, not registered"
        assert assembly_status().registered == ()

    def test_a_lost_record_reads_as_here_rather_than_as_never_registered(
        self, liulab_data: Path
    ) -> None:
        directory = _register(liulab_data / "genome", "sacCer3")
        (directory / ".completion.json").unlink()

        row = _rows(assembly_status())["sacCer3"]

        assert row.state == "here, not registered"
        assert not is_prepared("sacCer3")


class TestWhatTheTreeIsAllowedToHold:
    """The rule for telling an assembly from anything else filed in the same directory."""

    def test_files_and_hidden_entries_are_not_assemblies(self, liulab_data: Path) -> None:
        root = liulab_data / "genome"
        _register(root, "sacCer3")
        (root / "notes.txt").write_text("mine\n")
        (root / ".work").mkdir()

        listed = {row.assembly_name for row in assembly_status().assemblies}

        assert "notes.txt" not in listed
        assert ".work" not in listed
        assert "sacCer3" in listed

    def test_an_absent_tree_is_answered_as_an_empty_one(self, tmp_path: Path) -> None:
        status = assembly_status(root=tmp_path / "definitely-not-a-data-root")

        assert status.registered == ()
        assert tuple(row.assembly_name for row in status.assemblies) == _OFFERED


class TestTheRowIsOneShape:
    """One row shape whatever state a row is in, and ``state`` derived rather than stored."""

    def test_the_json_payload_carries_every_field_and_not_the_derived_state(
        self, liulab_data: Path
    ) -> None:
        _register(liulab_data / "genome", "sacCer3")

        status = assembly_status()
        payload = status.as_json()

        assert payload["directory"] == str(liulab_data / "genome")
        assert [row["assembly_name"] for row in payload["assemblies"]] == [
            row.assembly_name for row in status.assemblies
        ]
        row = next(row for row in payload["assemblies"] if row["assembly_name"] == "sacCer3")
        assert row == {
            "assembly_name": "sacCer3",
            "offered": True,
            "registered": True,
            "present": True,
            "directory": str(liulab_data / "genome" / "sacCer3"),
            "species": "Saccharomyces cerevisiae",
            "ucsc_name": "sacCer3",
            "ncbi_name": "R64-1-1",
            "source_url": "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz",
            "sha256": "6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3",
        }
        assert "state" not in row

    def test_an_absent_assembly_carries_no_directory(self) -> None:
        row = _rows(assembly_status())["hg38"]

        assert row.present is False
        assert row.directory is None
        # The table's columns are there whether or not this machine holds the assembly.
        assert row.species == "Homo sapiens"
        assert row.sha256 is not None


class TestTheRootOverride:
    """``root`` points the whole report at another tree, as ``cache_dir`` points at one."""

    def test_an_explicit_root_is_read_instead_of_the_layout(
        self, liulab_data: Path, tmp_path: Path
    ) -> None:
        elsewhere = tmp_path / "elsewhere"
        _register(elsewhere, "mm39")
        _register(liulab_data / "genome", "hg38")

        status = assembly_status(root=elsewhere)

        assert status.directory == elsewhere
        assert status.registered == ("mm39",)
        assert _rows(status)["hg38"].registered is False
