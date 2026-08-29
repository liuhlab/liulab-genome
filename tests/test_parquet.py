"""Tests for genome.tf.motif.parquet — the sink a scan too large to hold streams to.

One claim carries this module and everything here serves it: **what comes off the disk is
the table that would have come back in memory** — the same rows, the same column order, the
same compact dtypes down to the category order, and the same provenance. Anything less and
a saved scan and a live one could not be compared, which is the whole reason a scan is
allowed to go to disk at all.

Two things make that non-obvious and both are asserted rather than assumed. A categorical
column's index width follows its cardinality, so batches written one at a time would
otherwise disagree on the schema. And ``frame.attrs`` does not survive pandas' own Parquet
round trip, so the provenance is carried in the file's key-value metadata and put back by
:func:`~genome.tf.motif.parquet.read_hits` — :func:`pandas.read_parquet` alone gives the
rows and drops the meaning.

**Nothing here reprs a full hit table.** pandas overflows casting a ``float16`` column for
display, and this suite turns warnings into errors; the score column is compared as numbers
and never printed.

The unit lane, unmarked: nothing here needs a binary.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import pytest

from genome.tf.motif import (
    HIT_COLUMNS,
    HIT_DTYPES,
    HIT_PROVENANCE,
    MotifSet,
    hit_count,
    parse_transfac,
    provenance_of,
    read_hits,
)
from genome.tf.motif.parquet import HIT_PROVENANCE_KEY, write_hits
from genome.tf.motif.scan import empty_hits

from .test_scan import FIXTURE, MOTIF_FIXTURE, TOO_SHORT, read_records, rows, word_motif


@pytest.fixture
def motifs(data_dir: Path) -> MotifSet:
    """The ten committed JASPAR records as a plain **Motif set**."""
    return MotifSet(parse_transfac((data_dir / MOTIF_FIXTURE).read_text(encoding="utf-8")))


@pytest.fixture
def planted(data_dir: Path) -> Path:
    """The committed FASTA with motifs planted in it."""
    return data_dir / FIXTURE


@pytest.fixture
def planted_records(planted: Path) -> dict[str, str]:
    """The committed FASTA read into ``{name: bases}``, case as committed."""
    return read_records(planted.read_text(encoding="utf-8"))


@pytest.fixture
def in_memory(motifs: MotifSet, planted: Path) -> pd.DataFrame:
    """The committed FASTA scanned into a table, which the file has to match."""
    return motifs.scan_fasta(planted)


@pytest.fixture
def written(motifs: MotifSet, planted: Path, tmp_path: Path) -> Path:
    """The same scan, streamed to Parquet instead."""
    return motifs.scan_fasta(planted, output=tmp_path / "hits.parquet")


class TestAPathIsAnswerEnough:
    def test_passing_an_output_path_returns_the_path(self, written: Path, tmp_path: Path) -> None:
        assert written == tmp_path / "hits.parquet"
        assert written.is_file()

    def test_a_string_path_works_as_well_as_a_path(
        self, motifs: MotifSet, planted: Path, tmp_path: Path
    ) -> None:
        found = motifs.scan_fasta(planted, output=str(tmp_path / "hits.parquet"))
        assert isinstance(found, Path)

    def test_the_parent_directory_is_created(
        self, motifs: MotifSet, planted: Path, tmp_path: Path
    ) -> None:
        found = motifs.scan_fasta(planted, output=tmp_path / "runs" / "today" / "hits.parquet")
        assert found.is_file()

    def test_a_directory_is_refused_by_name(self, tmp_path: Path) -> None:
        with pytest.raises(IsADirectoryError, match="Name the file itself"):
            write_hits([empty_hits()], tmp_path, {})

    def test_all_three_entry_points_take_one(
        self, motifs: MotifSet, planted: Path, planted_records: dict[str, str], tmp_path: Path
    ) -> None:
        one = motifs.scan(planted_records["plantedI"], "plantedI", output=tmp_path / "a.parquet")
        many = motifs.scan_sequences(planted_records, output=tmp_path / "b.parquet")
        whole = motifs.scan_fasta(planted, output=tmp_path / "c.parquet")
        assert [path.name for path in (one, many, whole)] == ["a.parquet", "b.parquet", "c.parquet"]
        assert rows(read_hits(many)) == rows(read_hits(whole))

    def test_reading_a_file_that_is_not_there_says_how_to_make_one(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="Scan with output="):
            read_hits(tmp_path / "never-written.parquet")


class TestTheFileIsTheTable:
    def test_the_rows_are_the_rows_the_scan_would_have_returned(
        self, written: Path, in_memory: pd.DataFrame
    ) -> None:
        assert rows(read_hits(written)) == rows(in_memory)

    def test_the_frame_is_equal_dtypes_included(
        self, written: Path, in_memory: pd.DataFrame
    ) -> None:
        # The strongest form of the claim: pandas' own comparison, which checks the column
        # order, every dtype, and a category column's categories and their order.
        pd.testing.assert_frame_equal(read_hits(written), in_memory)

    def test_the_dtypes_are_the_compact_ones(self, written: Path) -> None:
        found = read_hits(written)
        assert {name: str(dtype) for name, dtype in found.dtypes.items()} == dict(HIT_DTYPES)
        assert list(found.columns) == list(HIT_COLUMNS)

    def test_the_score_survives_as_float16_and_not_as_something_wider(
        self, written: Path, in_memory: pd.DataFrame
    ) -> None:
        # 19 bytes a row against about 100 is the reason the dtypes are the contract; a
        # round trip that quietly widened the score column would cost most of that.
        found = read_hits(written)
        assert str(found["score"].dtype) == "float16"
        assert list(found["score"].to_numpy()) == list(in_memory["score"].to_numpy())

    def test_pandas_alone_gives_the_rows_and_drops_the_meaning(self, written: Path) -> None:
        # Why there is a reader here at all.
        assert pd.read_parquet(written).attrs == {}
        assert read_hits(written).attrs != {}

    def test_a_scan_that_found_nothing_still_writes_a_file(
        self, motifs: MotifSet, tmp_path: Path
    ) -> None:
        found = read_hits(motifs.scan("N" * 300, output=tmp_path / "none.parquet"))
        assert len(found) == 0
        assert {name: str(dtype) for name, dtype in found.dtypes.items()} == dict(HIT_DTYPES)
        assert found.attrs["motifs_skipped"] == TOO_SHORT

    def test_the_index_is_a_fresh_range(self, written: Path) -> None:
        found = read_hits(written)
        assert list(found.index) == list(range(len(found)))


class TestTheCategoriesComeBackInOrder:
    def test_categories_a_later_batch_introduced_are_still_sorted(self, tmp_path: Path) -> None:
        # Arrow keeps the categories in the order the row groups introduced them; the
        # in-memory table sorts them, and the two have to agree. Written out of order on
        # purpose: 'chrZ' arrives first and must still sort second.
        motifs = MotifSet([word_motif("GATTACAG")])
        peaks = {"chrZ": "TTTTGATTACAGTTTT", "chrA": "AAAAGATTACAGAAAA"}
        in_memory = motifs.scan_sequences(peaks)
        written = motifs.scan_sequences(peaks, output=tmp_path / "hits.parquet")
        assert list(in_memory["sequence_name"].cat.categories) == ["chrA", "chrZ"]
        pd.testing.assert_frame_equal(read_hits(written), in_memory)

    def test_more_categories_than_a_narrow_index_could_hold(self, tmp_path: Path) -> None:
        # A dictionary index widens with cardinality, so two batches would disagree on the
        # schema unless the writer pins one width. Three hundred names needs more than int8.
        motifs = MotifSet([word_motif("GATTACAG")])
        peaks = {f"peak{index:04d}": "TTTTGATTACAGTTTT" for index in range(300)}
        in_memory = motifs.scan_sequences(peaks)
        written = motifs.scan_sequences(peaks, output=tmp_path / "hits.parquet")
        assert len(in_memory) == 300
        pd.testing.assert_frame_equal(read_hits(written), in_memory)


class TestTheProvenanceTravels:
    def test_every_provenance_key_arrives(self, written: Path) -> None:
        assert set(read_hits(written).attrs) == set(HIT_PROVENANCE)

    def test_it_is_the_provenance_the_in_memory_table_carries(
        self, written: Path, in_memory: pd.DataFrame
    ) -> None:
        assert read_hits(written).attrs == in_memory.attrs

    def test_the_tuples_come_back_as_tuples_and_not_as_lists(self, written: Path) -> None:
        # JSON has one sequence type; a background that read back as a list would not
        # compare equal to the one a live scan records, and the two could not be reconciled.
        found = read_hits(written)
        assert isinstance(found.attrs["background"], tuple)
        assert isinstance(found.attrs["motifs_scanned"], tuple)
        assert isinstance(found.attrs["motifs_skipped"], tuple)

    def test_a_derived_background_survives_the_round_trip(
        self, motifs: MotifSet, planted_records: dict[str, str], tmp_path: Path
    ) -> None:
        peaks = {f"peak{index}": planted_records["plantedI"] for index in range(20)}
        in_memory = motifs.scan_sequences(peaks)
        written = motifs.scan_sequences(peaks, output=tmp_path / "hits.parquet")
        assert in_memory.attrs["background"] != (0.25, 0.25, 0.25, 0.25)
        assert read_hits(written).attrs["background"] == in_memory.attrs["background"]

    def test_a_release_still_names_itself_from_disk(
        self, motifs: MotifSet, planted: Path, tmp_path: Path
    ) -> None:
        found = read_hits(motifs.scan_fasta(planted, output=tmp_path / "hits.parquet"))
        assert (found.attrs["release"], found.attrs["tax_group"]) == (None, None)

    def test_it_lives_under_its_own_key_in_the_file(self, written: Path) -> None:
        assert HIT_PROVENANCE_KEY in (pq.read_schema(written).metadata or {})

    def test_a_parquet_written_by_something_else_simply_carries_none(
        self, in_memory: pd.DataFrame, tmp_path: Path
    ) -> None:
        path = tmp_path / "foreign.parquet"
        in_memory.to_parquet(path)
        assert read_hits(path).attrs == {}


class TestNothingIsMaterialised:
    def test_the_batches_are_written_as_they_arrive(self, tmp_path: Path) -> None:
        # The sink drains the same per-sequence iterator the collector does. Draining it
        # lazily is what keeps a 550-million-row scan off the heap, so the writer must not
        # be handed a list — assert it consumed a generator one batch at a time.
        drained: list[int] = []

        def batches() -> Iterator[pd.DataFrame]:
            for index in range(3):
                assert len(drained) == index  # nothing was pulled ahead of being asked for
                drained.append(index)
                yield empty_hits()

        write_hits(batches(), tmp_path / "hits.parquet", {})
        assert drained == [0, 1, 2]

    def test_a_scan_never_builds_the_table_it_is_streaming(
        self, motifs: MotifSet, planted: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # The collector is what concatenates every batch into one frame; streaming must not
        # reach it, or "the whole result is never materialised" would not be true.
        from genome.tf.motif import scan as scan_mod

        def refuse(*args: object, **kwargs: object) -> pd.DataFrame:
            raise AssertionError("the whole table was collected")

        monkeypatch.setattr(scan_mod, "_collect", refuse)
        assert motifs.scan_fasta(planted, output=tmp_path / "hits.parquet").is_file()


class TestNoGuardAndNoRefusal:
    def test_a_large_result_is_written_without_complaint(self, tmp_path: Path) -> None:
        # A genome-scale scan is the caller's decision. Nothing here counts the rows and
        # decides they are too many — the only limit is the disk.
        motifs = MotifSet([word_motif("GATTACAG")])
        peaks = {f"peak{index:05d}": "TTTTGATTACAGTTTT" for index in range(2000)}
        found = read_hits(motifs.scan_sequences(peaks, output=tmp_path / "hits.parquet"))
        assert len(found) == len(motifs.scan_sequences(peaks)) == 2000


class TestTheFooterAnswersWithoutTheRows:
    """What a written scan was, and how much of it there is, read off the footer alone.

    The case a written table exists for is the case where reading it back is fatal — 550
    million rows is hg38 against a full vertebrate release — so the two facts a finished
    scan is summarised from must cost the same there as on an empty file.
    """

    def test_the_provenance_is_the_one_read_hits_puts_on_the_frame(self, written: Path) -> None:
        assert provenance_of(written) == read_hits(written).attrs

    def test_the_count_is_the_number_of_rows_the_table_holds(self, written: Path) -> None:
        assert hit_count(written) == len(read_hits(written)) > 0

    def test_a_scan_that_found_nothing_counts_zero_rather_than_raising(
        self, motifs: MotifSet, tmp_path: Path
    ) -> None:
        found = motifs.scan("N" * 400, output=tmp_path / "none.parquet")
        assert hit_count(found) == 0
        assert set(provenance_of(found)) == set(HIT_PROVENANCE)

    def test_neither_reads_a_row(self, written: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        # The claim itself: pandas' reader is the only thing here that materialises rows,
        # so making it raise proves neither of these went near one.
        def refuse(*args: object, **kwargs: object) -> pd.DataFrame:
            raise AssertionError("the rows were read back")

        monkeypatch.setattr(pd, "read_parquet", refuse)
        assert hit_count(written) > 0
        assert provenance_of(written)["threshold"] == 1e-4

    def test_a_parquet_written_by_something_else_carries_no_provenance(
        self, in_memory: pd.DataFrame, tmp_path: Path
    ) -> None:
        path = tmp_path / "foreign.parquet"
        in_memory.to_parquet(path)
        assert provenance_of(path) == {}
        assert hit_count(path) == len(in_memory)

    @pytest.mark.parametrize("reader", [provenance_of, hit_count])
    def test_a_file_that_is_not_there_says_how_to_make_one(
        self, reader: object, tmp_path: Path
    ) -> None:
        with pytest.raises(FileNotFoundError, match="Scan with output="):
            reader(tmp_path / "never-written.parquet")  # type: ignore[operator]
