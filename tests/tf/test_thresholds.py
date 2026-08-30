"""Tests for genome.tf.motif.thresholds — the disk cache of per-motif cutoffs.

Converting a **Threshold** into one score per matrix is the engine's one slow step and a
pure function of ``(matrices, background, p)``, so what matters here is that the cache is
keyed on exactly that triple and on nothing else: the same triple must not recompute, a
different one must not collide, and two backgrounds that the 0.001 grid makes equal must
land on one entry — the whole reason the grid exists.

The engine is called for real, as everywhere else in this feature. What is counted is how
often: :func:`MOODS.tools.threshold_from_p` is wrapped rather than replaced, so a hit is
observed as an engine call that did not happen rather than as an internal that was reached.

The **Data dir** is a temporary one in every test, from the autouse fixture in
``conftest.py``, so the cache these tests fill is their own.

The unit lane, unmarked: the scan engine is a library, not a binary.
"""

from __future__ import annotations

from pathlib import Path

import MOODS.tools
import numpy as np
import pytest

from genome.tf.motif.motif import Motif
from genome.tf.motif.thresholds import cutoffs_for, threshold_cache_dir

#: The background every call here shares unless it is the thing being varied.
UNIFORM = (0.25, 0.25, 0.25, 0.25)


def word_matrix(bases: str) -> list[list[float]]:
    """The log-odds matrix, **in nats**, of a motif fixed on every base of ``bases``."""
    counts = np.zeros((4, len(bases)))
    for column, base in enumerate(bases):
        counts["ACGT".index(base), column] = 100.0
    return (Motif("MA9999.1", "Testin", counts).log_odds() * np.log(2)).tolist()


@pytest.fixture
def matrices() -> list[list[list[float]]]:
    """Two matrices, as a scan's doubled forward-and-reverse list would hand them over."""
    return [word_matrix("GATTACAG"), word_matrix("CTGTAATC")]


@pytest.fixture
def engine_calls(monkeypatch: pytest.MonkeyPatch) -> list[float]:
    """Record every real conversion, without standing in for it — the engine still runs."""
    calls: list[float] = []
    real = MOODS.tools.threshold_from_p

    def counted(matrix: object, background: object, p: float, *args: object) -> float:
        calls.append(p)
        return real(matrix, background, p, *args)

    monkeypatch.setattr(MOODS.tools, "threshold_from_p", counted)
    return calls


def entries() -> list[Path]:
    """Every cache entry currently on disk."""
    directory = threshold_cache_dir()
    return sorted(directory.glob("*.json")) if directory.is_dir() else []


class TestWhereTheCacheLives:
    def test_it_lives_under_the_motif_subtree_and_only_appears_on_first_use(
        self, liulab_data: Path, matrices: list[list[list[float]]]
    ) -> None:
        # A motif belongs to no assembly, so neither does anything derived from one.
        assert threshold_cache_dir() == liulab_data / "motif" / "thresholds"
        assert not threshold_cache_dir().exists()  # asking where it is creates nothing
        cutoffs_for(matrices, UNIFORM, 1e-4)
        assert len(entries()) == 1


class TestTheCacheIsHitOnARepeat:
    def test_a_repeat_does_not_reach_the_engine_and_yields_one_cutoff_per_matrix(
        self, matrices: list[list[list[float]]], engine_calls: list[float]
    ) -> None:
        first = cutoffs_for(matrices, UNIFORM, 1e-4)
        assert len(engine_calls) == len(matrices) == len(first)
        second = cutoffs_for(matrices, UNIFORM, 1e-4)
        assert len(engine_calls) == len(matrices)  # nothing was converted the second time
        assert second == first
        # One cutoff per matrix, in the matrices' own order — the index split that
        # recovers Strand depends on this alignment.
        assert len(cutoffs_for(matrices[:1], UNIFORM, 1e-4)) == 1

    def test_the_answer_off_the_disk_is_the_answer_the_engine_gave(
        self, matrices: list[list[list[float]]]
    ) -> None:
        computed = cutoffs_for(matrices, UNIFORM, 1e-4)
        for path in entries():
            path.unlink()
        assert cutoffs_for(matrices, UNIFORM, 1e-4) == computed


class TestWhatTheKeyIs:
    def test_each_part_of_the_key_makes_a_different_entry(
        self, matrices: list[list[list[float]]], engine_calls: list[float]
    ) -> None:
        # Each call below changes exactly one axis away from this same baseline, so a key
        # that silently dropped one axis (background, say) would collide with the
        # baseline here rather than growing the entry count — which is what each assert
        # below would catch.
        baseline = cutoffs_for(matrices, UNIFORM, 1e-4)
        assert len(entries()) == 1

        strict = cutoffs_for(matrices, UNIFORM, 1e-6)  # p differs; background, matrices don't
        assert len(entries()) == 2
        assert all(a > b for a, b in zip(strict, baseline, strict=True))

        cutoffs_for(matrices, (0.3, 0.2, 0.2, 0.3), 1e-4)  # background differs; p, matrices don't
        assert len(entries()) == 3

        cutoffs_for([word_matrix("GATTACAGTC")], UNIFORM, 1e-4)  # matrices differ; p, bg don't
        assert len(entries()) == 4

    def test_two_backgrounds_the_grid_makes_equal_share_one_entry(
        self, matrices: list[list[list[float]]], engine_calls: list[float]
    ) -> None:
        # Two peak sets from one genome, whose compositions agree to three decimals.
        cutoffs_for(matrices, (0.3, 0.2, 0.2, 0.3), 1e-4)
        cutoffs_for(matrices, (0.3, 0.2, 0.2, 0.3), 1e-4)
        assert len(entries()) == 1
        assert len(engine_calls) == len(matrices)


class TestACacheIsNeverADependency:
    def test_a_bad_entry_on_disk_is_a_miss_and_not_an_error(
        self, matrices: list[list[list[float]]]
    ) -> None:
        computed = cutoffs_for(matrices, UNIFORM, 1e-4)
        (entry,) = entries()
        # Garbage, a truncated (one cutoff where two matrices were asked about, which would
        # silently mis-call a strand), and an older format all have to read as a miss.
        for payload in (
            "{not json at all",
            '{"version": 1, "cutoffs": [0.0]}',
            '{"version": 0, "cutoffs": [1.0, 1.0]}',
        ):
            entry.write_text(payload, encoding="utf-8")
            assert cutoffs_for(matrices, UNIFORM, 1e-4) == computed

    def test_a_cache_that_cannot_be_written_makes_the_scan_slow_and_not_broken(
        self, liulab_data: Path, matrices: list[list[list[float]]], engine_calls: list[float]
    ) -> None:
        # A read-only data root is an ordinary cluster situation; it must not fail a scan.
        (liulab_data / "motif").write_text("not a directory", encoding="utf-8")
        first = cutoffs_for(matrices, UNIFORM, 1e-4)
        second = cutoffs_for(matrices, UNIFORM, 1e-4)
        assert first == second
        assert len(engine_calls) == 2 * len(matrices)  # every call paid for itself
