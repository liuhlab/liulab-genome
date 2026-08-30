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
really was cut, rather than passing because nothing happened. A real worker pool costs real
process-spawn time, so the agreement claims below share as few pool spawns as the distinct
claims allow: one test asserts several agreement properties from a single spawn rather than
paying that cost once per property.

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

from .test_scan import FIXTURE, MOTIF_FIXTURE, read_records, rows, word_motif

# Every test here spawns its own worker processes. Under `--dist=loadgroup` that pins
# them to ONE xdist worker, so they run one at a time rather than eight of them forking
# pools at once — which oversubscribed the box and made the lane's wall bimodal, 13.5 s
# or 16.1 s depending purely on how they happened to be scheduled.
pytestmark = pytest.mark.xdist_group("spawns_parallel")


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
    def test_plan_shards_fits_cuts_with_overlap_and_widens_a_too_short_shard(self) -> None:
        assert plan_shards(100, overlap=14, shard_length=1000) == [(0, 100, 100)]
        assert plan_shards(100, overlap=14, shard_length=40) == [
            (0, 54, 40),
            (40, 94, 40),
            (80, 100, 20),
        ]
        # A shard shorter than the overlap is raised to clear it, or it would be almost
        # entirely overlap and the work would double.
        plan = plan_shards(100, overlap=30, shard_length=4)
        assert [owned for _offset, _stop, owned in plan][:2] == [31, 31]

    @given(
        length=st.integers(min_value=0, max_value=5000),
        overlap=st.integers(min_value=0, max_value=60),
        shard_length=st.integers(min_value=1, max_value=500),
    )
    def test_the_owned_regions_partition_the_sequence_and_each_piece_clears_the_overlap(
        self, length: int, overlap: int, shard_length: int
    ) -> None:
        # What makes a hit starting at the last owned position scorable over every base it
        # covers: the piece runs to the end of that hit, or to the end of the sequence.
        plan = plan_shards(length, overlap, shard_length)
        covered = 0
        for offset, stop, owned in plan:
            assert offset == covered  # contiguous, in order, no gap and no repeat
            assert stop == min(offset + owned + overlap, length)
            covered += owned
        assert covered == length

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
    def test_two_workers_agree_with_serial_rows_and_provenance(
        self, motifs: MotifSet, planted: Path
    ) -> None:
        serial = motifs.scan_fasta(planted)
        shared = motifs.scan_fasta(planted, workers=2)
        pd.testing.assert_frame_equal(shared, serial)
        assert shared.attrs == serial.attrs

    @pytest.mark.parametrize("shard_length", [100, 110, 37])
    def test_two_workers_agree_when_a_record_is_cut_across_a_shard_boundary(
        self,
        motifs: MotifSet,
        planted: Path,
        monkeypatch: pytest.MonkeyPatch,
        shard_length: int,
    ) -> None:
        # 100 puts each planted site exactly on a shard boundary; 110 puts one across one;
        # 37 is deliberately awkward against every motif length, cutting the record into
        # many shards at once. `assert_frame_equal` is row-for-row, so it already proves a
        # boundary site is reported exactly once and that sharding finds the same sites a
        # whole record does — a duplicate or a drop would change the row count or order.
        serial = motifs.scan_fasta(planted)
        shard_lengths(monkeypatch, shard_length)
        longest = max(len(motif) for motif in motifs)
        assert len(plan_shards(600, longest - 1, shard_length)) > 1  # it really was cut
        shared = motifs.scan_fasta(planted, workers=2)
        pd.testing.assert_frame_equal(shared, serial)

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
    def test_no_process_starts_for_one_worker_none_or_a_resolved_count_of_one(
        self, motifs: MotifSet, planted: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # Importing this package and scanning must never start a pool: under spawn, a pool
        # re-imports the caller's script, and an unguarded one would re-execute itself.
        def refuse(*args: object, **kwargs: object) -> object:
            raise AssertionError("a process pool was started")

        monkeypatch.setattr(parallel_mod, "ProcessPoolExecutor", refuse)
        assert len(motifs.scan_fasta(planted)) > 0
        assert len(motifs.scan_fasta(planted, workers=1)) > 0

        # None means "work it out"; inside a one-CPU allocation that is one, and one is
        # serial — so this asserts the resolution is really wired into the scan.
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
        assert rows(motifs.scan_fasta(planted, workers=None)) == rows(motifs.scan_fasta(planted))

    def test_zero_workers_is_refused_before_anything_is_read(self, motifs: MotifSet) -> None:
        with pytest.raises(ValueError, match="at least 1"):
            motifs.scan("ACGT" * 20, workers=0)


# ---------------------------------------------------------------------------
# The batches the parallel source hands back
# ---------------------------------------------------------------------------


class TestTheBatchesAreStillOnePerSequence:
    def test_a_silent_sequence_contributes_nothing_and_order_survives_sharding(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        shard_lengths(monkeypatch, 20)
        motifs = MotifSet([word_motif("GATTACAG")])
        peaks = {"quiet": "N" * 200} | {
            f"peak{index}": "TTTT" + "GATTACAG" + "TTTT" * 20 for index in range(6)
        }
        serial = motifs.scan_sequences(peaks)
        shared = motifs.scan_sequences(peaks, workers=2)
        assert rows(shared) == rows(serial)
        assert set(shared["sequence_name"]) == set(peaks) - {"quiet"}
        assert list(dict.fromkeys(shared["sequence_name"])) == [
            name for name in peaks if name != "quiet"
        ]

    def test_an_empty_motif_set_scans_in_parallel_to_an_empty_table(self, planted: Path) -> None:
        found = MotifSet([]).scan_fasta(planted, workers=2)
        assert len(found) == 0
        assert found.attrs["motifs_scanned"] == ()
