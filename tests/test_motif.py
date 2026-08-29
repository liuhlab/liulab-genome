"""Tests for genome.tf.motif.motif — the Motif type, with nothing else loaded.

No file, no download, no network, and no scan engine: a motif is built from a count
matrix here and asked everything about itself. Hand-checkable cases for the answers a
reader can verify by eye — a consensus, a spacer that survives trimming — and hypothesis
for the invariants that must hold over every matrix: the counts to probabilities to
information chain, the bits bounds, and what trimming may and may not do to a length and
an offset.

The unit lane, unmarked. Plotting draws under a non-interactive backend, so no test opens
a window, and every figure a test makes is closed again — the suite turns warnings into
errors, and matplotlib warns once twenty are left open.
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

matplotlib.use("Agg")  # headless: no test may open a window

import matplotlib.pyplot as plt  # must follow the backend choice above
from matplotlib.axes import Axes

from genome.seq import DNA
from genome.tf.motif import BASES, MIN_MOTIF_LENGTH, Motif

# ---------------------------------------------------------------------------
# Hand-built columns. Two kinds is all most of these tests need: one that says
# everything (2 bits) and one that says nothing (0 bits).
# ---------------------------------------------------------------------------


def fixed(base: str, count: float = 20.0) -> np.ndarray:
    """One column with every observation on ``base`` — 2 bits."""
    column = np.zeros((4, 1))
    column[BASES.index(base), 0] = count
    return column


def flat(count: float = 5.0, positions: int = 1) -> np.ndarray:
    """``positions`` columns with the four bases equally observed — 0 bits each."""
    return np.full((4, positions), count)


def word(bases: str, count: float = 20.0) -> np.ndarray:
    """A count matrix spelling ``bases``, every position fixed."""
    return np.hstack([fixed(base, count) for base in bases])


def motif(counts: np.ndarray, **kwargs: object) -> Motif:
    """A motif over ``counts`` with a colourless id and name."""
    return Motif("MA9999.1", "Testin", counts, **kwargs)  # type: ignore[arg-type]


#: One column saying a little — about 0.43 bits, so it survives the default trim
#: threshold of 0.25 and not a threshold of 1.
WEAK = np.array([[6.0], [2.0], [1.0], [1.0]])


#: A dimer: two fixed half-sites either side of a spacer that says nothing, wrapped in
#: flanks that say nothing either. Trimming must take the flanks and leave the spacer.
DIMER_FLANK = 3
DIMER_HALF = "GATTACAG"
DIMER_SPACER = 4
DIMER = np.hstack(
    [
        flat(positions=DIMER_FLANK),
        word(DIMER_HALF),
        flat(positions=DIMER_SPACER),
        word(DIMER_HALF),
        flat(positions=DIMER_FLANK),
    ]
)


# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

_count = st.floats(
    min_value=0.0,
    max_value=1000.0,
    allow_nan=False,
    allow_infinity=False,
    allow_subnormal=False,
)


def _observed(column: list[float]) -> list[float]:
    """Give a column that observed nothing one observation of each base."""
    return column if sum(column) > 0 else [1.0, 1.0, 1.0, 1.0]


@st.composite
def count_matrices(draw: st.DrawFn, *, min_length: int = 1, max_length: int = 14) -> np.ndarray:
    """A valid 4 x L count matrix: finite, non-negative, every column observed."""
    length = draw(st.integers(min_value=min_length, max_value=max_length))
    columns = draw(
        st.lists(
            st.lists(_count, min_size=4, max_size=4).map(_observed),
            min_size=length,
            max_size=length,
        )
    )
    return np.array(columns, dtype=float).T


@st.composite
def motifs(draw: st.DrawFn, *, min_length: int = 1, max_length: int = 14) -> Motif:
    """A motif over a valid count matrix, with an arbitrary offset already on it."""
    counts = draw(count_matrices(min_length=min_length, max_length=max_length))
    offset = draw(st.integers(min_value=0, max_value=20))
    return motif(counts, offset=offset)


# ---------------------------------------------------------------------------
# Construction, and what a bad matrix says
# ---------------------------------------------------------------------------


class TestConstruction:
    def test_fields_are_kept_as_given(self) -> None:
        m = Motif(
            "MA0139.2",
            "CTCF",
            word("ACGTACG"),
            offset=4,
            tax_group="vertebrates",
            tf_class="C2H2 zinc finger factors",
            tf_family="More than 3 adjacent zinc fingers",
            uniprot_ids=("P49711",),
            pubmed_ids=("17512414", "27924024"),
            data_type="ChIP-seq",
        )
        assert (m.motif_id, m.motif_name, m.offset) == ("MA0139.2", "CTCF", 4)
        assert (m.tax_group, m.data_type) == ("vertebrates", "ChIP-seq")
        assert m.uniprot_ids == ("P49711",)
        assert m.pubmed_ids == ("17512414", "27924024")

    def test_offset_defaults_to_zero_and_annotations_to_nothing_stated(self) -> None:
        m = motif(word("ACGT"))
        assert m.offset == 0
        assert (m.tax_group, m.tf_class, m.tf_family, m.data_type) == ("", "", "", "")
        assert (m.uniprot_ids, m.pubmed_ids) == ((), ())

    def test_counts_become_float_even_from_integers(self) -> None:
        # The source carries fractional counts, so the matrix is float, never int.
        m = motif(np.array([[9, 1], [1, 1], [0, 7], [0, 1]]))
        assert m.counts.dtype == np.float64

    def test_fractional_counts_are_kept(self) -> None:
        m = motif(np.array([[9.5, 1.0], [0.5, 1.0], [0.0, 7.5], [0.0, 0.5]]))
        assert m.counts[0, 0] == 9.5

    def test_counts_are_a_copy_the_caller_cannot_reach(self) -> None:
        given_counts = word("ACGT")
        m = motif(given_counts)
        given_counts[0, 0] = 999.0
        assert m.counts[0, 0] == 20.0

    def test_counts_are_read_only(self) -> None:
        m = motif(word("ACGT"))
        assert m.counts.flags.writeable is False
        with pytest.raises(ValueError, match="read-only"):
            m.counts[0, 0] = 1.0

    def test_frozen(self) -> None:
        m = motif(word("ACGT"))
        with pytest.raises(AttributeError):
            m.motif_id = "MA0001.1"  # type: ignore[misc]

    def test_empty_id_raises_naming_what_an_id_is_for(self) -> None:
        with pytest.raises(ValueError, match="addressed by its id"):
            Motif("", "CTCF", word("ACGT"))

    def test_negative_offset_raises(self) -> None:
        with pytest.raises(ValueError, match="offset must be >= 0"):
            motif(word("ACGT"), offset=-1)

    @pytest.mark.parametrize(
        "bad",
        [
            np.ones((3, 5)),
            np.ones((5, 5)),
            np.ones((4, 0)),
            np.ones(4),
            np.ones((4, 5, 2)),
        ],
    )
    def test_wrong_shape_raises_and_names_the_transpose(self, bad: np.ndarray) -> None:
        with pytest.raises(ValueError, match="4 x L"):
            motif(bad)

    def test_positions_as_rows_is_caught_by_the_shape_check(self) -> None:
        # The one mistake this layout invites, since logomaker wants the transpose.
        with pytest.raises(ValueError, match="Transpose it"):
            motif(word("ACGTACG").T)

    def test_negative_counts_raise_and_name_the_position(self) -> None:
        counts = word("ACGT").copy()
        counts[0, 2] = -1.0
        with pytest.raises(ValueError, match=r"negative counts at position\(s\) 2"):
            motif(counts)

    def test_non_finite_counts_raise_and_name_the_position(self) -> None:
        counts = word("ACGT").copy()
        counts[0, 1] = np.nan
        with pytest.raises(ValueError, match=r"non-finite values at position\(s\) 1"):
            motif(counts)

    def test_a_position_with_no_observations_raises(self) -> None:
        counts = np.hstack([word("AC"), np.zeros((4, 1))])
        with pytest.raises(ValueError, match=r"no observations at position\(s\) 2"):
            motif(counts)

    def test_a_bare_string_of_ids_is_refused(self) -> None:
        # "P49711" would otherwise be stored letter by letter.
        with pytest.raises(ValueError, match="iterable of ids"):
            motif(word("ACGT"), uniprot_ids="P49711")

    def test_ids_are_frozen_into_a_tuple(self) -> None:
        m = motif(word("ACGT"), pubmed_ids=["1", "2"])
        assert m.pubmed_ids == ("1", "2")


# ---------------------------------------------------------------------------
# Identity: length, equality, hashing, repr
# ---------------------------------------------------------------------------


class TestIdentity:
    def test_len_and_length_agree_and_count_positions(self) -> None:
        m = motif(word("ACGTACG"))
        assert len(m) == m.length == 7

    def test_equal_motifs_compare_equal_and_hash_alike(self) -> None:
        a = motif(word("ACGT"), tax_group="vertebrates")
        b = motif(word("ACGT"), tax_group="vertebrates")
        assert a == b
        assert hash(a) == hash(b)
        assert len({a, b}) == 1

    def test_a_different_matrix_is_a_different_motif(self) -> None:
        assert motif(word("ACGT")) != motif(word("ACGA"))

    def test_a_different_annotation_is_a_different_motif(self) -> None:
        assert motif(word("ACGT"), tax_group="plants") != motif(word("ACGT"))

    def test_a_different_offset_is_a_different_motif(self) -> None:
        assert motif(word("ACGT"), offset=1) != motif(word("ACGT"))

    def test_comparison_with_a_non_motif_is_false_rather_than_an_error(self) -> None:
        assert motif(word("ACGT")) != 5
        assert motif(word("ACGT")) != "MA9999.1"

    def test_usable_as_a_dict_key(self) -> None:
        m = motif(word("ACGT"))
        assert {m: "kept"}[motif(word("ACGT"))] == "kept"

    def test_repr_names_the_identity_and_not_the_matrix(self) -> None:
        assert repr(motif(word("ACGT"), offset=2)) == (
            "Motif(motif_id='MA9999.1', motif_name='Testin', length=4, offset=2)"
        )


# ---------------------------------------------------------------------------
# Derived from the counts: probabilities, log-odds, information, consensus
# ---------------------------------------------------------------------------


class TestProbabilities:
    def test_columns_are_normalised(self) -> None:
        m = motif(np.array([[3.0], [1.0], [0.0], [0.0]]))
        assert m.probabilities.ravel().tolist() == [0.75, 0.25, 0.0, 0.0]

    def test_a_fresh_array_every_time_so_nothing_derived_is_stored(self) -> None:
        m = motif(word("ACGT"))
        first = m.probabilities
        first[0, 0] = 999.0
        assert m.probabilities[0, 0] != 999.0


class TestLogOdds:
    def test_uniform_background_scores_a_fixed_base_at_two_bits(self) -> None:
        m = motif(fixed("A", 100.0))
        assert m.log_odds()[0, 0] == pytest.approx(2.0, abs=1e-3)

    def test_a_gc_poor_background_raises_the_score_of_a_fixed_a(self) -> None:
        m = motif(fixed("A", 100.0))
        assert m.log_odds([0.4, 0.1, 0.1, 0.4])[0, 0] == pytest.approx(1.32, abs=1e-2)

    def test_the_background_is_an_argument_and_never_a_field(self) -> None:
        # One motif, two backgrounds, still one motif — which is the whole reason the
        # background is not stored on it.
        m = motif(word("ACGTACG"))
        uniform = m.log_odds()
        skewed = m.log_odds([0.4, 0.1, 0.1, 0.4])
        assert not np.allclose(uniform, skewed)
        assert m == motif(word("ACGTACG"))
        assert not hasattr(m, "background")

    def test_shape_matches_the_count_matrix(self) -> None:
        m = motif(word("ACGTACG"))
        assert m.log_odds().shape == m.counts.shape

    def test_a_pseudocount_keeps_an_unobserved_base_finite(self) -> None:
        assert np.isfinite(motif(fixed("A")).log_odds()).all()

    @pytest.mark.parametrize(
        ("background", "message"),
        [
            ([0.25, 0.25, 0.25], "4 frequencies"),
            ([0.25, 0.25, 0.25, 0.25, 0.0], "4 frequencies"),
            ([0.5, 0.5, 0.0, 0.0], "must all be > 0"),
            ([0.5, 0.5, 0.5, 0.5], "must sum to 1"),
        ],
    )
    def test_a_background_that_is_not_four_frequencies_raises(
        self, background: list[float], message: str
    ) -> None:
        with pytest.raises(ValueError, match=message):
            motif(word("ACGT")).log_odds(background)

    def test_a_pseudocount_of_zero_raises(self) -> None:
        with pytest.raises(ValueError, match="pseudocount must be > 0"):
            motif(word("ACGT")).log_odds(pseudocount=0.0)


class TestInformationContent:
    def test_a_fixed_position_says_two_bits(self) -> None:
        assert motif(fixed("A")).information_content.tolist() == [2.0]

    def test_a_uniform_position_says_nothing(self) -> None:
        assert motif(flat()).information_content.tolist() == [0.0]

    def test_one_value_per_position(self) -> None:
        m = motif(np.hstack([word("GATTACA"), flat(positions=3)]))
        assert m.information_content.shape == (10,)
        assert m.information_content.tolist() == [2.0] * 7 + [0.0] * 3

    def test_a_half_and_half_position_says_one_bit(self) -> None:
        m = motif(np.array([[5.0], [5.0], [0.0], [0.0]]))
        assert m.information_content[0] == pytest.approx(1.0)


class TestConsensus:
    def test_returns_a_typed_dna_and_never_a_str(self) -> None:
        result = motif(word("GATTACA")).consensus
        assert isinstance(result, DNA)
        assert result == DNA("GATTACA")

    def test_reads_the_most_observed_base_at_each_position(self) -> None:
        counts = np.array(
            [
                [9.0, 1.0, 0.0],  # A
                [1.0, 1.0, 2.0],  # C
                [0.0, 7.0, 3.0],  # G
                [0.0, 1.0, 5.0],  # T
            ]
        )
        assert motif(counts).consensus == DNA("AGT")

    def test_a_tie_goes_to_the_first_base_in_order(self) -> None:
        assert motif(flat(positions=3)).consensus == DNA("AAA")

    def test_an_uninformative_position_still_gets_a_letter(self) -> None:
        # No IUPAC code is ever produced: this package's alphabet is ACGT.
        m = motif(np.hstack([word("GA"), flat()]))
        assert m.consensus == DNA("GAA")
        assert DNA.outside_alphabet(m.consensus) == []


# ---------------------------------------------------------------------------
# Trimming
# ---------------------------------------------------------------------------


class TestTrim:
    def test_the_minimum_length_is_seven(self) -> None:
        # A 6-mer's best possible match has p = 2.44e-4 and cannot reach the default
        # 1e-4 threshold, so trimming must never produce one.
        assert MIN_MOTIF_LENGTH == 7

    def test_flanks_go_and_the_spacer_stays(self) -> None:
        trimmed = motif(DIMER).trim()
        assert len(trimmed) == 2 * len(DIMER_HALF) + DIMER_SPACER
        spacer = trimmed.information_content[len(DIMER_HALF) : len(DIMER_HALF) + DIMER_SPACER]
        assert spacer.tolist() == [0.0] * DIMER_SPACER
        assert trimmed.consensus == DNA(f"{DIMER_HALF}AAAA{DIMER_HALF}")

    def test_a_trimmed_motif_keeps_its_id_its_name_and_its_annotation(self) -> None:
        full = Motif(
            "MA0139.2",
            "CTCF",
            DIMER,
            tax_group="vertebrates",
            uniprot_ids=("P49711",),
            data_type="ChIP-seq",
        )
        trimmed = full.trim()
        assert (trimmed.motif_id, trimmed.motif_name) == ("MA0139.2", "CTCF")
        assert (trimmed.tax_group, trimmed.uniprot_ids, trimmed.data_type) == (
            "vertebrates",
            ("P49711",),
            "ChIP-seq",
        )

    def test_the_offset_maps_a_position_back_into_the_full_frame(self) -> None:
        full = motif(DIMER)
        trimmed = full.trim()
        assert trimmed.offset == DIMER_FLANK
        for position in range(len(trimmed)):
            assert (
                trimmed.information_content[position]
                == (full.information_content[position + trimmed.offset])
            )

    def test_trimming_a_trimmed_motif_composes(self) -> None:
        full = motif(
            np.hstack([flat(positions=2), WEAK, word("GATTACAG"), WEAK, flat(positions=2)])
        )
        once = full.trim(0.25)  # takes the two silent columns off each end
        twice = once.trim(1.0)  # takes the near-silent one off each end as well
        assert (len(once), once.offset) == (10, 2)
        assert (len(twice), twice.offset) == (8, 3)
        assert np.array_equal(twice.counts, full.counts[:, twice.offset : twice.offset + 8])

    def test_nothing_to_trim_returns_the_same_motif(self) -> None:
        m = motif(word("GATTACA"))
        assert m.trim() is m

    def test_a_motif_already_below_the_minimum_is_returned_untouched(self) -> None:
        m = motif(flat(positions=5))
        assert m.trim() is m

    def test_it_stops_at_the_minimum_rather_than_going_under(self) -> None:
        # Every position says nothing, so the threshold alone would take all of them.
        trimmed = motif(flat(positions=20)).trim()
        assert len(trimmed) == MIN_MOTIF_LENGTH

    def test_a_single_informative_position_still_yields_a_scannable_motif(self) -> None:
        m = motif(np.hstack([flat(positions=6), fixed("G"), flat(positions=6)]))
        trimmed = m.trim()
        assert len(trimmed) == MIN_MOTIF_LENGTH
        assert trimmed.information_content.max() == 2.0

    def test_the_maximum_length_is_honoured(self) -> None:
        trimmed = motif(word("GATTACAGATTACAGATTAC")).trim(max_length=10)
        assert len(trimmed) == 10

    def test_the_maximum_drops_the_less_informative_end(self) -> None:
        m = motif(np.hstack([np.repeat(WEAK, 5, axis=1), word("CAGATTAC")]))
        trimmed = m.trim(max_length=8)
        assert trimmed.consensus == DNA("CAGATTAC")
        assert trimmed.offset == 5

    def test_the_maximum_keeps_the_front_when_the_two_ends_say_the_same(self) -> None:
        trimmed = motif(word("GATTACAGATTACAGATTAC")).trim(max_length=8)
        assert (trimmed.consensus, trimmed.offset) == (DNA("GATTACAG"), 0)

    def test_a_maximum_below_the_minimum_raises(self) -> None:
        with pytest.raises(ValueError, match="below min_length"):
            motif(word("GATTACAG")).trim(max_length=3)

    def test_a_minimum_below_one_raises(self) -> None:
        with pytest.raises(ValueError, match="min_length must be >= 1"):
            motif(word("GATTACAG")).trim(min_length=0)

    def test_a_higher_threshold_trims_more(self) -> None:
        m = motif(np.hstack([WEAK, word("GATTACAG"), WEAK]))
        assert len(m.trim(0.25)) == 10
        assert len(m.trim(1.0)) == 8


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


class TestPlot:
    def test_returns_an_axes_and_makes_a_figure_when_given_none(self) -> None:
        axes = motif(word("GATTACA")).plot()
        try:
            assert isinstance(axes, Axes)
            assert axes.figure is not None
        finally:
            plt.close(axes.figure)  # type: ignore[arg-type]

    def test_draws_into_the_axes_it_is_handed(self) -> None:
        figure, axes = plt.subplots()
        try:
            returned = motif(word("GATTACA")).plot(ax=axes)
            assert returned is axes
            assert len(axes.patches) > 0
        finally:
            plt.close(figure)

    def test_a_caller_s_other_axes_is_left_alone(self) -> None:
        figure, (left, right) = plt.subplots(1, 2)
        try:
            motif(word("GATTACA")).plot(ax=left)
            assert len(left.patches) > 0
            assert len(right.patches) == 0
        finally:
            plt.close(figure)

    def test_the_drawn_extent_matches_the_motif_length(self) -> None:
        for length in (7, 12):
            figure, axes = plt.subplots()
            try:
                motif(word("GATTACAGATTAC"[:length])).plot(ax=axes)
                low, high = axes.get_xlim()
                assert high - low == pytest.approx(length)
            finally:
                plt.close(figure)

    def test_the_y_axis_is_information_content_in_bits(self) -> None:
        figure, axes = plt.subplots()
        try:
            motif(word("GATTACA")).plot(ax=axes)
            assert axes.get_ylabel() == "Information content (bits)"
            assert axes.get_ylim() == (0.0, 2.0)
        finally:
            plt.close(figure)

    def test_the_height_drawn_is_the_quantity_trimming_thresholds_on(self) -> None:
        # Column by column, the stacked glyphs are as tall as the information content.
        m = motif(np.hstack([word("GAT"), flat(), np.array([[5.0], [5.0], [0.0], [0.0]])]))
        figure, axes = plt.subplots()
        try:
            m.plot(ax=axes)
            drawn = np.zeros(len(m))
            for patch in axes.patches:
                bounds = patch.get_extents().transformed(axes.transData.inverted())
                drawn[round(bounds.x0 + bounds.width / 2)] += bounds.height
            assert drawn == pytest.approx(m.information_content, abs=1e-6)
        finally:
            plt.close(figure)

    def test_logomaker_keywords_are_passed_through(self) -> None:
        figure, axes = plt.subplots()
        try:
            motif(word("GATTACA")).plot(ax=axes, color_scheme="classic")
            assert len(axes.patches) > 0
        finally:
            plt.close(figure)


# ---------------------------------------------------------------------------
# Properties — the invariants that must hold over every matrix
# ---------------------------------------------------------------------------


class TestPropertiesCountsToInformation:
    @given(counts=count_matrices())
    def test_every_column_of_probabilities_sums_to_one(self, counts: np.ndarray) -> None:
        probabilities = motif(counts).probabilities
        assert probabilities.shape == counts.shape
        assert probabilities.sum(axis=0) == pytest.approx(np.ones(counts.shape[1]))
        assert (probabilities >= 0.0).all()

    @given(counts=count_matrices())
    def test_probabilities_stay_proportional_to_the_counts(self, counts: np.ndarray) -> None:
        m = motif(counts)
        assert m.probabilities * m.counts.sum(axis=0) == pytest.approx(m.counts, rel=1e-9)

    @given(counts=count_matrices(), factor=st.sampled_from([0.5, 2.0, 8.0]))
    def test_scaling_every_count_changes_nothing_derived(
        self, counts: np.ndarray, factor: float
    ) -> None:
        # Counts are observations: twice as many of them is the same motif.
        original, scaled = motif(counts), motif(counts * factor)
        assert scaled.probabilities == pytest.approx(original.probabilities, abs=1e-12)
        assert scaled.information_content == pytest.approx(original.information_content, abs=1e-9)
        assert scaled.consensus == original.consensus

    @given(counts=count_matrices())
    def test_information_is_in_bits_between_zero_and_two(self, counts: np.ndarray) -> None:
        information = motif(counts).information_content
        assert information.shape == (counts.shape[1],)
        assert (information >= 0.0).all()
        assert (information <= 2.0).all()

    @given(count=st.floats(min_value=1e-6, max_value=1e6, allow_nan=False))
    def test_a_uniform_column_is_exactly_zero_bits(self, count: float) -> None:
        assert motif(np.full((4, 1), count)).information_content.tolist() == [0.0]

    @given(base=st.sampled_from(BASES), count=st.floats(min_value=1e-6, max_value=1e6))
    def test_a_fixed_column_is_exactly_two_bits(self, base: str, count: float) -> None:
        assert motif(fixed(base, count)).information_content.tolist() == [2.0]

    @given(counts=count_matrices(), order=st.permutations(range(4)))
    def test_information_does_not_depend_on_which_base_it_is(
        self, counts: np.ndarray, order: list[int]
    ) -> None:
        permuted = motif(counts[list(order), :]).information_content
        assert permuted == pytest.approx(motif(counts).information_content, abs=1e-12)

    @given(counts=count_matrices())
    def test_consensus_is_one_base_per_position_and_each_is_a_column_maximum(
        self, counts: np.ndarray
    ) -> None:
        m = motif(counts)
        consensus = m.consensus
        assert isinstance(consensus, DNA)
        assert len(consensus) == len(m)
        for position, base in enumerate(consensus):
            assert m.counts[BASES.index(base), position] == m.counts[:, position].max()

    @given(counts=count_matrices())
    def test_log_odds_is_finite_and_shaped_like_the_matrix(self, counts: np.ndarray) -> None:
        scores = motif(counts).log_odds()
        assert scores.shape == counts.shape
        assert np.isfinite(scores).all()


class TestPropertiesTrim:
    @given(m=motifs(), threshold=st.floats(min_value=0.0, max_value=2.0))
    def test_never_longer_than_the_input(self, m: Motif, threshold: float) -> None:
        assert len(m.trim(threshold)) <= len(m)

    @given(m=motifs(), threshold=st.floats(min_value=0.0, max_value=2.0))
    def test_never_below_the_minimum_it_could_reach(self, m: Motif, threshold: float) -> None:
        # A motif already shorter than the minimum is handed back as it is, so the floor
        # for any one call is whichever of the two is smaller.
        assert len(m.trim(threshold)) >= min(MIN_MOTIF_LENGTH, len(m))

    @given(m=motifs(), threshold=st.floats(min_value=0.0, max_value=2.0))
    def test_the_offset_is_the_slice_that_was_actually_taken(
        self, m: Motif, threshold: float
    ) -> None:
        trimmed = m.trim(threshold)
        taken = trimmed.offset - m.offset
        assert taken >= 0
        assert taken + len(trimmed) <= len(m)
        assert np.array_equal(trimmed.counts, m.counts[:, taken : taken + len(trimmed)])

    @given(m=motifs(), threshold=st.floats(min_value=0.0, max_value=2.0))
    def test_the_identity_and_the_annotation_survive(self, m: Motif, threshold: float) -> None:
        trimmed = m.trim(threshold)
        assert (trimmed.motif_id, trimmed.motif_name) == (m.motif_id, m.motif_name)
        assert trimmed.tax_group == m.tax_group

    @given(
        m=motifs(min_length=1, max_length=20),
        threshold=st.floats(min_value=0.0, max_value=2.0),
        max_length=st.integers(min_value=MIN_MOTIF_LENGTH, max_value=20),
    )
    def test_a_maximum_length_is_never_exceeded(
        self, m: Motif, threshold: float, max_length: int
    ) -> None:
        assert len(m.trim(threshold, max_length=max_length)) <= max_length

    @given(m=motifs(), threshold=st.floats(min_value=0.0, max_value=2.0))
    def test_trimming_twice_composes_into_one_slice_of_the_full_frame(
        self, m: Motif, threshold: float
    ) -> None:
        twice = m.trim(threshold).trim(threshold)
        taken = twice.offset - m.offset
        assert np.array_equal(twice.counts, m.counts[:, taken : taken + len(twice)])
