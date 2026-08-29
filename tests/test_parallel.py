"""Tests for genome.tf.motif.parallel — scanning across processes, and the shard arithmetic.

**One property carries this file: serial and parallel produce the identical table.** Not
the same hits as a set — the identical frame, row for row, with the same dtypes and the same
provenance, so that choosing two workers is a choice about wall time and about nothing else.
Everything else here exists because that property can fail quietly: a hit lying across a
shard boundary reported twice, or not at all; a shard's positions reported in its own frame
instead of the sequence's; the shards of one sequence handed back interleaved rather than in
the order a serial scan would have emitted them.

:func:`~genome.tf.motif.parallel.plan_shards` is the arithmetic on its own, and is tested
property-based with no process anywhere near it: the owned regions must **partition** the
sequence, and every piece must run far enough past its own region to score a hit starting at
its last owned position.

The end-to-end tests use **two workers on the 600-base fixture**, which is what keeps them
cheap enough to sit in the unit lane. The default shard length is five megabases, so the
fixture would never be cut; the tests that are about cutting set it small and assert that it
really was cut, rather than passing because nothing happened.

The unit lane, unmarked: a process pool is not a binary this package ships.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.tf.motif import MotifSet, parse_transfac, read_hits
from genome.tf.motif import parallel as parallel_mod
from genome.tf.motif.parallel import plan_shards

from .test_scan import FIXTURE, MOTIF_FIXTURE, PLANTED, read_records, rows, sites, word_motif


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


def shard_lengths(monkeypatch: pytest.MonkeyPatch, length: int) -> None:
    """Cut sequences into ``length``-base pieces, so a 600-base record is really sharded."""
    monkeypatch.setattr(parallel_mod, "_SHARD_LENGTH", length)


# ---------------------------------------------------------------------------
# The shard arithmetic, on its own
# ---------------------------------------------------------------------------


class TestPlanShards:
    def test_a_sequence_that_fits_is_one_shard(self) -> None:
        assert plan_shards(100, overlap=14, shard_length=1000) == [(0, 100, 100)]

    def test_a_longer_one_is_cut_with_an_overlap(self) -> None:
        assert plan_shards(100, overlap=14, shard_length=40) == [
            (0, 54, 40),
            (40, 94, 40),
            (80, 100, 20),
        ]

    def test_a_shard_shorter_than_the_overlap_is_raised_to_clear_it(self) -> None:
        # Otherwise a shard would be almost entirely overlap and the work would double.
        plan = plan_shards(100, overlap=30, shard_length=4)
        assert [owned for _offset, _stop, owned in plan][:2] == [31, 31]

    @given(
        length=st.integers(min_value=0, max_value=5000),
        overlap=st.integers(min_value=0, max_value=60),
        shard_length=st.integers(min_value=1, max_value=500),
    )
    def test_the_owned_regions_partition_the_sequence(
        self, length: int, overlap: int, shard_length: int
    ) -> None:
        plan = plan_shards(length, overlap, shard_length)
        covered = 0
        for offset, _stop, owned in plan:
            assert offset == covered  # contiguous, in order, no gap and no repeat
            covered += owned
        assert covered == length

    @given(
        length=st.integers(min_value=1, max_value=5000),
        overlap=st.integers(min_value=0, max_value=60),
        shard_length=st.integers(min_value=1, max_value=500),
    )
    def test_every_piece_reaches_past_its_own_region_by_the_overlap(
        self, length: int, overlap: int, shard_length: int
    ) -> None:
        # What makes a hit starting at the last owned position scorable over every base it
        # covers: the piece runs to the end of that hit, or to the end of the sequence.
        for offset, stop, owned in plan_shards(length, overlap, shard_length):
            assert stop == min(offset + owned + overlap, length)
            assert stop > offset

    @given(
        length=st.integers(min_value=1, max_value=2000),
        motif_length=st.integers(min_value=7, max_value=40),
        shard_length=st.integers(min_value=1, max_value=300),
        start=st.integers(min_value=0, max_value=1999),
    )
    def test_a_hit_anywhere_is_owned_by_exactly_one_shard_that_can_see_it(
        self, length: int, motif_length: int, shard_length: int, start: int
    ) -> None:
        # The whole correctness claim, stated over every possible hit: exactly one shard
        # keeps it, and that shard's piece holds every base of it.
        if start + motif_length > length:
            return
        plan = plan_shards(length, motif_length - 1, shard_length)
        owners = [stop for offset, stop, owned in plan if offset <= start < offset + owned]
        assert len(owners) == 1
        assert start + motif_length <= owners[0]


# ---------------------------------------------------------------------------
# Serial and parallel, which must be the same table
# ---------------------------------------------------------------------------


class TestSerialAndParallelAgree:
    def test_two_workers_and_one_agree_when_sequences_are_split_between_them(
        self, motifs: MotifSet, planted: Path
    ) -> None:
        serial = motifs.scan_fasta(planted)
        shared = motifs.scan_fasta(planted, workers=2)
        pd.testing.assert_frame_equal(shared, serial)
        assert shared.attrs == serial.attrs

    @pytest.mark.parametrize("shard_length", [100, 110])
    def test_two_workers_and_one_agree_when_a_record_is_cut_up(
        self,
        motifs: MotifSet,
        planted: Path,
        monkeypatch: pytest.MonkeyPatch,
        shard_length: int,
    ) -> None:
        # 100 puts each planted site exactly on a shard boundary; 110 puts one across one.
        serial = motifs.scan_fasta(planted)
        shard_lengths(monkeypatch, shard_length)
        longest = max(len(motif) for motif in motifs)
        assert len(plan_shards(600, longest - 1, shard_length)) > 1  # it really was cut
        shared = motifs.scan_fasta(planted, workers=2)
        pd.testing.assert_frame_equal(shared, serial)

    def test_a_site_on_a_shard_boundary_is_reported_exactly_once(
        self, motifs: MotifSet, planted: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        shard_lengths(monkeypatch, 110)  # the CTCF site at [100, 115) spans this boundary
        shared = motifs.scan_fasta(planted, workers=2)
        for record, start, end, strand, motif_id, _word in PLANTED:
            found = [row for row in rows(shared) if row[:1] == (motif_id,)]
            wanted = (motif_id, record, start, end, strand)
            assert [row for row in found if row[2:6] == wanted[1:]] != []
            assert sum(1 for row in found if row[2:6] == wanted[1:]) == 1

    def test_the_provenance_is_the_same_provenance(self, motifs: MotifSet, planted: Path) -> None:
        assert motifs.scan_fasta(planted, workers=2).attrs == motifs.scan_fasta(planted).attrs

    def test_a_parallel_scan_streams_to_parquet_the_same_way(
        self, motifs: MotifSet, planted: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # The two features compose: the sink drains the parallel source without knowing it.
        serial = motifs.scan_fasta(planted)
        shard_lengths(monkeypatch, 100)
        written = motifs.scan_fasta(planted, workers=2, output=tmp_path / "hits.parquet")
        pd.testing.assert_frame_equal(read_hits(written), serial)
        assert read_hits(written).attrs == serial.attrs


# ---------------------------------------------------------------------------
# One worker starts none
# ---------------------------------------------------------------------------


class TestOneWorkerIsSerial:
    def test_the_default_scan_starts_no_process(
        self, motifs: MotifSet, planted: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # Importing this package and scanning must never start a pool: under spawn, a pool
        # re-imports the caller's script, and an unguarded one would re-execute itself.
        def refuse(*args: object, **kwargs: object) -> object:
            raise AssertionError("a process pool was started")

        monkeypatch.setattr(parallel_mod, "ProcessPoolExecutor", refuse)
        assert len(motifs.scan_fasta(planted)) > 0
        assert len(motifs.scan_fasta(planted, workers=1)) > 0

    def test_zero_workers_is_refused_before_anything_is_read(self, motifs: MotifSet) -> None:
        with pytest.raises(ValueError, match="at least 1"):
            motifs.scan("ACGT" * 20, workers=0)

    def test_a_resolved_count_is_used(
        self, motifs: MotifSet, planted: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # None means "work it out"; inside a one-CPU allocation that is one, and one is
        # serial — so this asserts the resolution is really wired into the scan.
        def refuse(*args: object, **kwargs: object) -> object:
            raise AssertionError("a process pool was started")

        monkeypatch.setattr(parallel_mod, "ProcessPoolExecutor", refuse)
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
        assert rows(motifs.scan_fasta(planted, workers=None)) == rows(motifs.scan_fasta(planted))


# ---------------------------------------------------------------------------
# The batches the parallel source hands back
# ---------------------------------------------------------------------------


class TestTheBatchesAreStillOnePerSequence:
    def test_a_sequence_with_no_hit_contributes_nothing_and_breaks_nothing(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        shard_lengths(monkeypatch, 20)
        motifs = MotifSet([word_motif("GATTACAG")])
        peaks = {"quiet": "N" * 200, "loud": "TTTT" + "GATTACAG" + "TTTT"}
        shared = motifs.scan_sequences(peaks, workers=2)
        assert rows(shared) == rows(motifs.scan_sequences(peaks))
        assert set(shared["sequence_name"]) == {"loud"}

    def test_the_sequences_come_back_in_the_order_they_arrived(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        shard_lengths(monkeypatch, 20)
        motifs = MotifSet([word_motif("GATTACAG")])
        peaks = {f"peak{index}": "TTTT" + "GATTACAG" + "TTTT" * 20 for index in range(6)}
        shared = motifs.scan_sequences(peaks, workers=2)
        assert list(dict.fromkeys(shared["sequence_name"])) == list(peaks)

    def test_an_empty_motif_set_scans_in_parallel_to_an_empty_table(self, planted: Path) -> None:
        found = MotifSet([]).scan_fasta(planted, workers=2)
        assert len(found) == 0
        assert found.attrs["motifs_scanned"] == ()

    def test_sharding_finds_the_same_sites_a_whole_record_does(
        self, motifs: MotifSet, planted_records: dict[str, str], monkeypatch: pytest.MonkeyPatch
    ) -> None:
        whole = sites(motifs.scan_sequences(planted_records))
        shard_lengths(monkeypatch, 37)  # deliberately awkward against every motif length
        assert sites(motifs.scan_sequences(planted_records, workers=2)) == whole
