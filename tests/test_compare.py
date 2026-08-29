"""Tests for genome.tf.motif.compare — asking what a motif looks like.

No network: the ten motifs come from the committed transfac fixture, parsed in
process. **tomtom is called for real**, exactly as the scan engine is — it is a declared
dependency, deterministic, and fast enough on ten motifs that faking it would only test
the adapter against our own assumptions about the engine.

The split between the two kinds of test here is deliberate. The *engine* tests hand real
motifs to tomtom and assert what comes back. The *pure* tests hand a hand-built
:class:`xarray.Dataset` to :class:`MotifComparison` and assert how it is ranked and
flattened, which is arithmetic on a labelled array and owes nothing to the engine — so
that is where hypothesis goes, and no property test pays for a tomtom call.

The unit lane, unmarked.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from hypothesis import given
from hypothesis import strategies as st

from genome.tf.motif import (
    Motif,
    MotifComparison,
    MotifSet,
    RaggedComparisonError,
    parse_transfac,
)
from genome.tf.motif.compare import COMPARISON_VARIABLES, FRAME_COLUMNS, _neg_log10

#: The committed fixture, cut from JASPAR 2024's `all` union file. Ten real motifs,
#: three of them named CTCF — one 15 columns long and two that embed it — which is what
#: makes it a real comparison rather than ten unrelated matrices.
FIXTURE = "tiny_jaspar_transfac.txt"

#: The nested CTCF pair. Both embed `MA0139.2`, and that is the point of naming them: a
#: long motif's strongest target by p-value can be a shorter motif nested inside it.
NESTED_CTCFS = ("MA1929.2", "MA1930.2")


@pytest.fixture
def fixture_motifs(data_dir: Path) -> tuple[Motif, ...]:
    """The committed transfac fixture, parsed — ten motifs, no file left open."""
    return parse_transfac((data_dir / FIXTURE).read_text(encoding="utf-8"))


@pytest.fixture
def targets(fixture_motifs: tuple[Motif, ...]) -> MotifSet:
    """The ten fixture motifs as the set every comparison here runs against."""
    return MotifSet(fixture_motifs)


# ---------------------------------------------------------------------------
# Hand-built comparisons, so the pure ranking and flattening can be exercised
# without the engine. The schema is spelled out here rather than imported so
# that the test is an independent statement of the contract.
# ---------------------------------------------------------------------------


def wide(
    neg_log10_p: list[list[float]],
    *,
    score: list[list[float]] | None = None,
    queries: list[str] | None = None,
    targets_: list[str] | None = None,
) -> MotifComparison:
    """A complete comparison over ``query x target``, from the p-values alone."""
    values = np.asarray(neg_log10_p, dtype=np.float16)
    n_queries, n_targets = values.shape
    query_ids = queries or [f"Q{index}" for index in range(n_queries)]
    target_ids = targets_ or [f"T{index}" for index in range(n_targets)]
    scores = np.asarray(score if score is not None else values, dtype=np.float32)
    dims = ("query", "target")
    return MotifComparison(
        xr.Dataset(
            {
                "neg_log10_p": (dims, values),
                "score": (dims, scores),
                "offset": (dims, np.zeros(values.shape, dtype=np.int16)),
                "overlap": (dims, np.ones(values.shape, dtype=np.int16)),
                "strand": (dims, np.full(values.shape, "+", dtype="<U1")),
            },
            coords={"query": query_ids, "target": target_ids},
            attrs={"targets_compared": n_targets},
        )
    )


def ragged(
    neg_log10_p: list[list[float]], target_ids: list[list[str]], *, targets_compared: int = 10
) -> MotifComparison:
    """A limited comparison over ``query x rank``, its target axis per query."""
    values = np.asarray(neg_log10_p, dtype=np.float16)
    n_queries, n_ranks = values.shape
    dims = ("query", "rank")
    return MotifComparison(
        xr.Dataset(
            {
                "target": (dims, np.asarray(target_ids)),
                "neg_log10_p": (dims, values),
                "score": (dims, values.astype(np.float32)),
                "offset": (dims, np.zeros(values.shape, dtype=np.int16)),
                "overlap": (dims, np.ones(values.shape, dtype=np.int16)),
                "strand": (dims, np.full(values.shape, "+", dtype="<U1")),
            },
            coords={
                "query": [f"Q{index}" for index in range(n_queries)],
                "rank": np.arange(n_ranks),
            },
            attrs={"targets_compared": targets_compared},
        )
    )


# ---------------------------------------------------------------------------
# Strategies — over p-value grids, never over motifs: no property test here
# calls the engine.
# ---------------------------------------------------------------------------

_p_value = st.floats(min_value=1e-300, max_value=1.0, allow_nan=False, allow_infinity=False)

_stored = st.floats(min_value=0.0, max_value=300.0, allow_nan=False, allow_infinity=False)


@st.composite
def grids(draw: st.DrawFn, *, max_queries: int = 4, max_targets: int = 5) -> list[list[float]]:
    """A ``query x target`` grid of stored negative log10 p-values."""
    n_queries = draw(st.integers(min_value=1, max_value=max_queries))
    n_targets = draw(st.integers(min_value=1, max_value=max_targets))
    return [
        draw(st.lists(_stored, min_size=n_targets, max_size=n_targets)) for _ in range(n_queries)
    ]


# ---------------------------------------------------------------------------
# What compare() accepts, and what it refuses
# ---------------------------------------------------------------------------


class TestWhatIsAccepted:
    def test_one_motif(self, targets: MotifSet) -> None:
        comparison = targets.compare(targets["MA0139.2"])
        assert comparison.query_ids == ("MA0139.2",)
        assert comparison.target_ids == targets.motif_ids

    def test_several_motifs(self, targets: MotifSet) -> None:
        comparison = targets.compare([targets["MA0139.2"], targets["SP1"]])
        assert comparison.query_ids == ("MA0139.2", "MA0079.5")

    def test_a_whole_motif_set(self, targets: MotifSet) -> None:
        comparison = targets.compare(targets)
        assert comparison.query_ids == targets.motif_ids
        assert comparison.data.sizes == {"query": 10, "target": 10}

    def test_a_de_novo_motif_that_no_release_published(self, targets: MotifSet) -> None:
        # The use case: a matrix out of a model, carrying none of the six annotations.
        de_novo = Motif("pattern_0", "", np.hstack([np.eye(4) * 20 + 1] * 2))
        assert targets.compare(de_novo).query_ids == ("pattern_0",)

    def test_no_queries_is_refused(self, targets: MotifSet) -> None:
        with pytest.raises(ValueError, match="no motifs to compare"):
            targets.compare([])

    def test_comparing_against_an_empty_set_is_refused(self, targets: MotifSet) -> None:
        with pytest.raises(ValueError, match="holds no motifs"):
            MotifSet([]).compare(targets["MA0139.2"])

    def test_two_queries_sharing_an_id_are_refused(self, targets: MotifSet) -> None:
        one = targets["MA0139.2"]
        with pytest.raises(ValueError, match="share the motif id"):
            targets.compare([one, one])

    def test_a_top_below_one_is_refused(self, targets: MotifSet) -> None:
        with pytest.raises(ValueError, match="top must be >= 1"):
            targets.compare(targets, top=0)

    def test_a_top_above_the_target_count_is_refused(self, targets: MotifSet) -> None:
        # Guards the engine: asking tomtom for more neighbours than there are targets
        # raises a SystemError out of numba, which names nothing a caller can act on.
        with pytest.raises(ValueError, match="only 10 targets"):
            targets.compare(targets["MA0139.2"], top=11)


# ---------------------------------------------------------------------------
# The anchor: a motif compared against itself
# ---------------------------------------------------------------------------


class TestAMotifAgainstItself:
    def test_a_motif_is_its_own_best_target(self, targets: MotifSet) -> None:
        best = targets.compare(targets["MA0139.2"]).to_frame()
        assert len(best) == 1
        assert best.loc[0, "target"] == "MA0139.2"

    def test_the_self_pair_is_a_gapless_forward_alignment_of_the_whole_motif(
        self, targets: MotifSet
    ) -> None:
        data = targets.compare(targets).data
        for motif in targets:
            pair = data.sel(query=motif.motif_id, target=motif.motif_id)
            assert int(pair["offset"]) == 0
            assert int(pair["overlap"]) == len(motif)
            assert str(pair["strand"].item()) == "+"

    def test_the_self_pair_is_the_strongest_score_in_its_row(self, targets: MotifSet) -> None:
        # True for all ten, unlike the p-value below: nothing aligns to a motif better
        # than the motif itself, column for column.
        data = targets.compare(targets).data
        strongest = data["score"].argmax(dim="target")
        assert list(data["target"][strongest].to_numpy()) == list(targets.motif_ids)

    def test_a_nested_motif_can_beat_the_self_pair_on_p_value(self, targets: MotifSet) -> None:
        # Not a defect and worth pinning: TOMTOM's p-value rewards a short dense
        # alignment, so the 15-column CTCF that both of these embed is ranked above
        # them by p even though they align to themselves perfectly. Documented on
        # MotifSet.compare, and the reason the anchor above ranks on score.
        frame = targets.compare(targets).to_frame()
        best = dict(zip(frame["query"], frame["target"], strict=True))
        assert [best[motif_id] for motif_id in NESTED_CTCFS] == ["MA0139.2", "MA0139.2"]
        unnested = set(targets.motif_ids) - set(NESTED_CTCFS)
        assert all(best[motif_id] == motif_id for motif_id in unnested)


# ---------------------------------------------------------------------------
# The labelled array, in both of its shapes
# ---------------------------------------------------------------------------


class TestTheLabelledArray:
    def test_the_complete_shape_is_query_by_target(self, targets: MotifSet) -> None:
        comparison = targets.compare(targets["MA0139.2"])
        assert comparison.is_ragged is False
        assert comparison.top is None
        assert comparison.data.sizes == {"query": 1, "target": 10}
        assert comparison.targets_compared == 10

    def test_the_limited_shape_is_query_by_rank(self, targets: MotifSet) -> None:
        comparison = targets.compare(targets, top=3)
        assert comparison.is_ragged is True
        assert comparison.top == 3
        assert comparison.data.sizes == {"query": 10, "rank": 3}
        assert comparison.targets_compared == 10

    def test_it_is_indexable_by_motif_id(self, targets: MotifSet) -> None:
        data = targets.compare(targets).data
        pair = data.sel(query="MA0139.2", target="MA1929.2")
        assert float(pair["neg_log10_p"]) > 0.0
        # And the same cell, reached the other way round, is a different comparison:
        # TOMTOM's p-value is not symmetric.
        mirrored = data.sel(query="MA1929.2", target="MA0139.2")
        assert float(mirrored["neg_log10_p"]) != float(pair["neg_log10_p"])

    def test_a_limited_array_is_still_indexable_by_query_id(self, targets: MotifSet) -> None:
        row = targets.compare(targets, top=2).data.sel(query="MA1930.2")
        assert list(row["target"].to_numpy()) == ["MA0139.2", "MA1930.2"]

    def test_a_limited_array_names_a_different_target_per_query(self, targets: MotifSet) -> None:
        # This is what ragged means: rank 0 is a different motif for each query, so
        # there is no shared target axis to index.
        best = targets.compare(targets, top=1).data["target"].to_numpy().ravel()
        assert len(set(best)) > 1

    def test_the_variables_and_their_dtypes_are_the_contract(self, targets: MotifSet) -> None:
        assert COMPARISON_VARIABLES == ("neg_log10_p", "score", "offset", "overlap", "strand")
        data = targets.compare(targets).data
        assert set(data.data_vars) == set(COMPARISON_VARIABLES)
        assert data["neg_log10_p"].dtype == np.float16
        assert data["score"].dtype == np.float32
        assert data["offset"].dtype == np.int16
        assert data["overlap"].dtype == np.int16
        assert data["strand"].dtype.kind == "U"

    def test_a_limited_array_carries_its_targets_as_a_variable(self, targets: MotifSet) -> None:
        data = targets.compare(targets, top=2).data
        assert set(data.data_vars) == {*COMPARISON_VARIABLES, "target"}
        assert data["target"].dims == ("query", "rank")

    def test_the_target_axis_of_a_limited_result_cannot_be_had(self, targets: MotifSet) -> None:
        comparison = targets.compare(targets, top=2)
        with pytest.raises(RaggedComparisonError, match="recompute"):
            _ = comparison.target_ids

    def test_a_dataset_that_is_neither_shape_is_refused(self) -> None:
        with pytest.raises(ValueError, match="query"):
            MotifComparison(xr.Dataset({"neg_log10_p": (("a", "b"), np.zeros((1, 1)))}))

    def test_repr_says_which_side_is_which(self, targets: MotifSet) -> None:
        assert repr(targets.compare(targets["MA0139.2"], top=2)) == (
            "MotifComparison(queries=1, targets=10, top=2)"
        )


# ---------------------------------------------------------------------------
# Negative log10 p, and why raw p is not stored
# ---------------------------------------------------------------------------


class TestNegativeLog10P:
    def test_a_p_value_half_precision_would_flush_to_zero_survives(self, targets: MotifSet) -> None:
        # The self-comparison of a real motif lands around p = 2e-14. Stored raw, in the
        # half precision this array uses, that is indistinguishable from zero; stored as
        # its negative log10 it is an ordinary small number.
        stored = float(
            targets.compare(targets["MA0139.2"])
            .data["neg_log10_p"]
            .sel(query="MA0139.2", target="MA0139.2")
        )
        raw = 10.0**-stored
        assert raw < np.finfo(np.float16).tiny
        assert np.float16(raw) == 0.0
        assert stored == pytest.approx(13.7, abs=0.5)

    def test_a_p_value_of_one_stores_as_zero(self) -> None:
        assert _neg_log10(np.array([1.0]))[0] == 0.0

    def test_an_underflowed_p_value_stores_as_infinite_and_warns_about_nothing(self) -> None:
        # The suite turns warnings into errors, so log10(0) must never be evaluated.
        assert np.isinf(_neg_log10(np.array([0.0]))[0])

    @given(p=_p_value)
    def test_it_is_never_negative(self, p: float) -> None:
        assert _neg_log10(np.array([p]))[0] >= 0.0

    @given(p=_p_value, q=_p_value)
    def test_the_smaller_p_value_stores_larger(self, p: float, q: float) -> None:
        smaller, larger = min(p, q), max(p, q)
        assert _neg_log10(np.array([smaller]))[0] >= _neg_log10(np.array([larger]))[0]

    @given(exponent=st.integers(min_value=0, max_value=300))
    def test_a_power_of_ten_round_trips_exactly(self, exponent: int) -> None:
        assert _neg_log10(np.array([10.0**-exponent]))[0] == pytest.approx(exponent, abs=1e-9)


# ---------------------------------------------------------------------------
# Flattening to one row per pair
# ---------------------------------------------------------------------------


class TestTheFlatFrame:
    def test_it_defaults_to_one_row_per_query(self, targets: MotifSet) -> None:
        frame = targets.compare(targets).to_frame()
        assert len(frame) == len(targets)
        assert list(frame["query"]) == list(targets.motif_ids)
        assert set(frame["rank"]) == {0}

    def test_the_columns_are_the_contract(self, targets: MotifSet) -> None:
        assert FRAME_COLUMNS == (
            "query",
            "target",
            "rank",
            "neg_log10_p",
            "score",
            "offset",
            "overlap",
            "strand",
        )
        assert tuple(targets.compare(targets).to_frame().columns) == FRAME_COLUMNS

    def test_a_larger_limit_gives_that_many_rows_per_query(self, targets: MotifSet) -> None:
        frame = targets.compare(targets).to_frame(top=3)
        assert len(frame) == 30
        assert list(frame["rank"][:3]) == [0, 1, 2]

    def test_no_limit_gives_every_pair(self, targets: MotifSet) -> None:
        assert len(targets.compare(targets).to_frame(top=None)) == 100

    def test_rows_are_ordered_best_first_within_each_query(self, targets: MotifSet) -> None:
        frame = targets.compare(targets).to_frame(top=None)
        for _, group in frame.groupby("query", sort=False):
            assert list(group["neg_log10_p"]) == sorted(group["neg_log10_p"], reverse=True)
            assert list(group["rank"]) == list(range(len(targets)))

    def test_the_frame_says_what_the_array_says(self, targets: MotifSet) -> None:
        comparison = targets.compare(targets)
        frame = comparison.to_frame(top=None)
        for row in frame.sample(12, random_state=0).to_dict("records"):
            cell = comparison.data.sel(query=row["query"], target=row["target"])
            assert float(cell["neg_log10_p"]) == float(row["neg_log10_p"])
            assert int(cell["offset"]) == row["offset"]
            assert str(cell["strand"].item()) == row["strand"]

    def test_a_limit_beyond_the_target_count_is_simply_every_pair(self, targets: MotifSet) -> None:
        # Nothing is missing from a complete comparison, so asking for more pairs than
        # exist is answered with the ones that do.
        assert len(targets.compare(targets).to_frame(top=99)) == 100

    def test_it_flattens_a_limited_comparison_too(self, targets: MotifSet) -> None:
        frame = targets.compare(targets, top=3).to_frame()
        assert len(frame) == 10
        assert tuple(frame.columns) == FRAME_COLUMNS

    def test_no_limit_on_a_limited_comparison_gives_what_was_kept(self, targets: MotifSet) -> None:
        assert len(targets.compare(targets, top=3).to_frame(top=None)) == 30

    def test_widening_a_limited_comparison_is_refused(self, targets: MotifSet) -> None:
        comparison = targets.compare(targets, top=3)
        with pytest.raises(RaggedComparisonError, match="top=5"):
            comparison.to_frame(top=5)

    def test_the_two_shapes_agree_on_the_best_target(self, targets: MotifSet) -> None:
        # The whole point of the faster path: it must answer the common question the
        # same way the complete comparison does.
        complete = targets.compare(targets).to_frame()
        limited = targets.compare(targets, top=1).to_frame()
        assert list(complete["target"]) == list(limited["target"])

    def test_a_limit_asking_for_every_target_is_not_a_widening(self, targets: MotifSet) -> None:
        # Ragged, but nothing is missing from it, so it answers rather than refusing.
        assert len(targets.compare(targets, top=10).to_frame(top=10)) == 100


# ---------------------------------------------------------------------------
# Ranking and flattening as arithmetic, over generated grids
# ---------------------------------------------------------------------------


class TestFlatteningProperties:
    def test_the_best_target_is_the_largest_stored_value(self) -> None:
        frame = wide([[1.0, 9.0, 4.0], [7.0, 2.0, 3.0]]).to_frame()
        assert list(frame["target"]) == ["T1", "T0"]

    def test_ties_are_broken_by_score_then_by_target_order(self) -> None:
        tied = wide([[5.0, 5.0, 5.0]], score=[[1.0, 3.0, 3.0]]).to_frame(top=None)
        assert list(tied["target"]) == ["T1", "T2", "T0"]

    @given(grid=grids())
    def test_every_pair_is_flattened_exactly_once(self, grid: list[list[float]]) -> None:
        comparison = wide(grid)
        frame = comparison.to_frame(top=None)
        assert len(frame) == len(grid) * len(grid[0])
        assert set(zip(frame["query"], frame["target"], strict=True)) == {
            (query, target) for query in comparison.query_ids for target in comparison.target_ids
        }

    @given(grid=grids(), top=st.integers(min_value=1, max_value=5))
    def test_the_row_count_is_the_limit_times_the_queries(
        self, grid: list[list[float]], top: int
    ) -> None:
        frame = wide(grid).to_frame(top=top)
        assert len(frame) == len(grid) * min(top, len(grid[0]))

    @given(grid=grids())
    def test_each_query_is_ranked_best_first(self, grid: list[list[float]]) -> None:
        frame = wide(grid).to_frame(top=None)
        for _, group in frame.groupby("query", sort=False):
            values = list(group["neg_log10_p"])
            assert values == sorted(values, reverse=True)

    @given(grid=grids())
    def test_the_default_row_is_the_row_the_full_frame_leads_with(
        self, grid: list[list[float]]
    ) -> None:
        comparison = wide(grid)
        leading = comparison.to_frame(top=None).groupby("query", sort=False).head(1)
        assert list(comparison.to_frame()["target"]) == list(leading["target"])


class TestRaggedFlatteningProperties:
    def test_a_hand_built_limited_comparison_keeps_its_own_order(self) -> None:
        comparison = ragged([[9.0, 2.0], [8.0, 7.0]], [["T3", "T1"], ["T0", "T3"]])
        frame = comparison.to_frame(top=None)
        assert list(frame["target"]) == ["T3", "T1", "T0", "T3"]
        assert list(frame["rank"]) == [0, 1, 0, 1]

    def test_it_refuses_to_widen_and_names_the_call_that_would(self) -> None:
        comparison = ragged([[9.0, 2.0]], [["T3", "T1"]], targets_compared=10)
        with pytest.raises(RaggedComparisonError) as raised:
            comparison.to_frame(top=4)
        message = str(raised.value)
        assert "top=4" in message
        assert "2" in message
        assert "compare" in message

    def test_a_limited_comparison_holding_every_target_widens_fine(self) -> None:
        comparison = ragged([[9.0, 2.0]], [["T1", "T0"]], targets_compared=2)
        assert len(comparison.to_frame(top=7)) == 2

    @given(grid=grids(max_targets=3))
    def test_every_kept_pair_is_flattened_exactly_once(self, grid: list[list[float]]) -> None:
        # Ordered as the engine returned them, which is best first per query.
        ordered = [sorted(row, reverse=True) for row in grid]
        names = [[f"T{index}" for index in range(len(grid[0]))] for _ in grid]
        frame = ragged(ordered, names, targets_compared=len(grid[0])).to_frame(top=None)
        assert len(frame) == len(grid) * len(grid[0])
        # Compared at the array's own precision: the frame is the array flattened, so a
        # value that half precision rounded is rounded in both.
        assert list(frame["neg_log10_p"]) == [np.float16(value) for row in ordered for value in row]


# ---------------------------------------------------------------------------
# The comparison as a whole answers the question it exists for
# ---------------------------------------------------------------------------


class TestNamingADeNovoPattern:
    def test_a_pattern_spelling_a_known_consensus_is_named_by_it(
        self, targets: MotifSet, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        # The use case, end to end: a matrix out of a model, compared against a set, and
        # named by what it looks like. Built from CTCF's own counts under another id, so
        # what comes back is checkable by eye rather than by whatever tomtom prefers.
        ctcf = targets["MA0139.2"]
        pattern = Motif("pattern_0", "", ctcf.counts)
        named = targets.compare(pattern).to_frame()
        assert named.loc[0, "query"] == "pattern_0"
        assert named.loc[0, "target"] == "MA0139.2"
        assert named.loc[0, "overlap"] == len(ctcf)

    def test_a_filtered_set_is_a_smaller_target_axis(self, targets: MotifSet) -> None:
        vertebrates = targets.filter(tax_group="vertebrates")
        comparison = vertebrates.compare(targets["MA0261.1"])
        assert comparison.target_ids == vertebrates.motif_ids
        assert comparison.targets_compared == len(vertebrates)

    def test_the_frame_is_a_dataframe_a_caller_can_join_on(self, targets: MotifSet) -> None:
        frame = targets.compare(targets).to_frame()
        assert isinstance(frame, pd.DataFrame)
        assert frame.index.tolist() == list(range(len(targets)))
