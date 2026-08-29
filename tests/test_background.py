"""Tests for genome.tf.motif.background — where the **Background** comes from.

The **Background** decides the answer more than any other scan parameter, which is why it
is derived rather than assumed and recorded rather than remembered. Two things are pinned
here that nothing else can pin: that the floor is what switches derivation on, and that
**deciding does not eat the input** — the records read while deciding come back, in order,
or a scan would silently skip its own first chromosome.

Quantisation onto the 0.001 grid is property-based, because "still four positive
frequencies summing to 1" has to hold for every background rather than for the three in a
table.

The unit lane, unmarked: nothing here needs a binary.
"""

from __future__ import annotations

from collections.abc import Iterator

import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.tf.motif.background import (
    _SAMPLE_CAP,
    BACKGROUND_FLOOR,
    BACKGROUND_MODES,
    UNIFORM_BACKGROUND,
    derive_background,
    quantise_background,
    resolve_background,
)

#: Four frequencies above zero summing to 1 — what a caller may pass as a background.
frequencies = st.lists(st.floats(min_value=0.05, max_value=10.0), min_size=4, max_size=4).map(
    lambda values: [value / sum(values) for value in values]
)


def records(bases: str, count: int) -> list[tuple[str, str]]:
    """``count`` named records, each holding ``bases``."""
    return [(f"record{index}", bases) for index in range(count)]


class TestTheFloor:
    def test_an_input_under_the_floor_falls_back_to_uniform(self) -> None:
        bases = "AAAACCGT" * 100
        assert len(bases) < BACKGROUND_FLOOR
        assert resolve_background(None, [("chrTest", bases)])[0] == UNIFORM_BACKGROUND

    def test_an_input_over_the_floor_is_derived_from(self) -> None:
        found, _ = resolve_background(None, [("chrTest", "AAAACCGT" * 2000)])
        assert found == (0.5, 0.25, 0.125, 0.125)

    def test_the_floor_is_ten_thousand_unambiguous_bases(self) -> None:
        assert BACKGROUND_FLOOR == 10_000
        under, _ = resolve_background(None, [("chrTest", "AC" * (BACKGROUND_FLOOR // 2 - 1))])
        over, _ = resolve_background(None, [("chrTest", "AC" * (BACKGROUND_FLOOR // 2))])
        assert under == UNIFORM_BACKGROUND
        assert over != UNIFORM_BACKGROUND

    def test_ambiguous_bases_do_not_count_towards_the_floor(self) -> None:
        # A run of N is not a composition, and must not be mistaken for enough of one.
        found, _ = resolve_background(None, [("chrTest", "N" * (2 * BACKGROUND_FLOOR))])
        assert found == UNIFORM_BACKGROUND

    def test_soft_masked_bases_count_as_the_bases_they_are(self) -> None:
        # A scan discards Soft-masking, so the composition it derives must not see it.
        masked, _ = resolve_background(None, [("chrTest", "aaaaccgt" * 2000)])
        upper, _ = resolve_background(None, [("chrTest", "AAAACCGT" * 2000)])
        assert masked == upper

    def test_many_small_records_reach_the_floor_together(self) -> None:
        # A peak set is thousands of short sequences; the floor is over the whole input.
        found, _ = resolve_background(None, records("AAAACCGT" * 10, 200))
        assert found == (0.5, 0.25, 0.125, 0.125)


class TestTheModesAndTheExplicitForm:
    def test_the_three_modes_are_auto_uniform_and_derive(self) -> None:
        assert BACKGROUND_MODES == ("auto", "uniform", "derive")

    def test_uniform_is_uniform_however_large_the_input(self) -> None:
        found, _ = resolve_background("uniform", [("chrTest", "AAAACCGT" * 2000)])
        assert found == UNIFORM_BACKGROUND

    def test_derive_derives_below_the_floor_too(self) -> None:
        found, _ = resolve_background("derive", [("chrTest", "AACC")])
        assert found == (0.375, 0.375, 0.125, 0.125)

    def test_auto_is_what_omitting_it_means(self) -> None:
        sequences = [("chrTest", "AAAACCGT" * 2000)]
        assert resolve_background(None, sequences)[0] == resolve_background("auto", sequences)[0]

    def test_an_explicit_four_tuple_is_used_as_given(self) -> None:
        found, _ = resolve_background([0.3, 0.2, 0.2, 0.3], [("chrTest", "AAAACCGT" * 2000)])
        assert found == (0.3, 0.2, 0.2, 0.3)

    def test_a_mode_that_does_not_exist_names_the_ones_that_do(self) -> None:
        with pytest.raises(ValueError, match="auto, uniform, derive"):
            resolve_background("gc", [("chrTest", "ACGT")])  # type: ignore[arg-type]

    @pytest.mark.parametrize(
        ("background", "message"),
        [
            ([0.25, 0.25, 0.25], "4 frequencies"),
            ([0.0, 0.25, 0.25, 0.5], "must all be > 0"),
            ([0.3, 0.3, 0.3, 0.3], "must sum to 1"),
        ],
    )
    def test_an_explicit_background_is_still_validated(
        self, background: list[float], message: str
    ) -> None:
        with pytest.raises(ValueError, match=message):
            resolve_background(background, [("chrTest", "ACGT")])


class TestTheSourceSurvives:
    def test_the_records_read_while_deciding_come_back_in_order(self) -> None:
        given_records = records("ACGTACGT" * 200, 20)
        _, rest = resolve_background(None, iter(given_records))
        assert list(rest) == given_records

    def test_a_generator_source_is_drained_exactly_once(self) -> None:
        given_records = records("ACGTACGT" * 400, 10)
        pulled: list[str] = []

        def source() -> Iterator[tuple[str, str]]:
            for name, bases in given_records:
                pulled.append(name)
                yield name, bases

        _, rest = resolve_background(None, source())
        assert len(pulled) < len(given_records)  # deciding stopped as soon as it had enough
        assert list(rest) == given_records
        assert pulled == [name for name, _ in given_records]  # each record yielded once

    @pytest.mark.parametrize(
        ("background", "expected"),
        [("uniform", UNIFORM_BACKGROUND), ([0.3, 0.2, 0.2, 0.3], (0.3, 0.2, 0.2, 0.3))],
    )
    def test_nothing_is_read_when_the_background_needs_no_input(
        self, background: object, expected: tuple[float, ...]
    ) -> None:
        def source() -> Iterator[tuple[str, str]]:
            raise AssertionError("the source was read")
            yield  # pragma: no cover - unreachable; it is here to make this a generator

        found, rest = resolve_background(background, source())  # type: ignore[arg-type]
        assert found == expected
        with pytest.raises(AssertionError, match="the source was read"):
            list(rest)

    def test_a_file_of_nothing_but_n_is_not_buffered_whole(self) -> None:
        # The floor counts unambiguous bases, so without a cap on how much is examined an
        # all-N file would be pulled into memory entire — which no scan may do.
        record_length = 10_000
        pulled: list[str] = []

        def source() -> Iterator[tuple[str, str]]:
            for index in range(1000):
                pulled.append(f"record{index}")
                yield f"record{index}", "N" * record_length

        found, rest = resolve_background(None, source())
        assert found == UNIFORM_BACKGROUND
        assert len(pulled) <= _SAMPLE_CAP // record_length + 1
        assert len(list(rest)) == 1000


class TestQuantisation:
    def test_a_background_already_on_the_grid_is_untouched(self) -> None:
        assert quantise_background([0.3, 0.2, 0.2, 0.3]) == (0.3, 0.2, 0.2, 0.3)
        assert quantise_background([0.25, 0.25, 0.25, 0.25]) == UNIFORM_BACKGROUND

    def test_a_background_off_the_grid_is_rounded_onto_it(self) -> None:
        assert quantise_background([1 / 3, 1 / 3, 1 / 6, 1 / 6]) == (0.333, 0.333, 0.167, 0.167)

    def test_two_backgrounds_within_the_grid_become_one(self) -> None:
        # What lets two peak sets from one genome share a threshold cache entry.
        assert quantise_background([0.3001, 0.1999, 0.2, 0.3]) == quantise_background(
            [0.2999, 0.2001, 0.2, 0.3]
        )

    @given(background=frequencies)
    def test_the_result_is_always_four_positive_frequencies_summing_to_one(
        self, background: list[float]
    ) -> None:
        found = quantise_background(background)
        assert len(found) == 4
        assert all(value > 0 for value in found)
        assert sum(found) == pytest.approx(1.0)

    @given(background=frequencies)
    def test_every_frequency_lands_on_a_thousandth(self, background: list[float]) -> None:
        for value in quantise_background(background):
            assert round(value * 1000) == pytest.approx(value * 1000)

    @given(background=frequencies)
    def test_quantising_is_idempotent(self, background: list[float]) -> None:
        once = quantise_background(background)
        assert quantise_background(once) == once

    @given(background=frequencies)
    def test_nothing_moves_further_than_the_noise_the_floor_allows(
        self, background: list[float]
    ) -> None:
        # The floor is set where the standard error on a frequency is about 0.004; rounding
        # onto the grid has to stay inside that, or it would be the larger perturbation.
        for wanted, found in zip(background, quantise_background(background), strict=True):
            assert abs(found - wanted) < 0.004


class TestDerivation:
    def test_counts_become_the_frequencies_they_are(self) -> None:
        assert derive_background([250, 250, 250, 250]) == UNIFORM_BACKGROUND

    def test_a_base_never_seen_still_gets_a_positive_frequency(self) -> None:
        # A frequency of zero makes every log-odds against that base infinite.
        found = derive_background([1000, 0, 0, 1000])
        assert all(value > 0 for value in found)
        assert sum(found) == pytest.approx(1.0)

    def test_the_pseudocount_is_invisible_at_the_floor(self) -> None:
        # One extra observation per base moves a frequency by under a tenth of the grid.
        assert derive_background([BACKGROUND_FLOOR // 4] * 4) == UNIFORM_BACKGROUND

    @given(
        counts=st.lists(st.integers(min_value=0, max_value=10**7), min_size=4, max_size=4).filter(
            lambda values: sum(values) > 0
        )
    )
    def test_a_derived_background_is_always_a_background(self, counts: list[int]) -> None:
        found = derive_background(counts)
        assert len(found) == 4
        assert all(value > 0 for value in found)
        assert sum(found) == pytest.approx(1.0)
