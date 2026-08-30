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
    def test_the_floor_decides_uniform_or_derived_and_what_counts_towards_it(self) -> None:
        bases = "AAAACCGT" * 100
        assert len(bases) < BACKGROUND_FLOOR
        assert resolve_background(None, [("chrTest", bases)])[0] == UNIFORM_BACKGROUND
        over, _ = resolve_background(None, [("chrTest", "AAAACCGT" * 2000)])
        assert over == (0.5, 0.25, 0.125, 0.125)

        assert BACKGROUND_FLOOR == 10_000
        under, _ = resolve_background(None, [("chrTest", "AC" * (BACKGROUND_FLOOR // 2 - 1))])
        exactly, _ = resolve_background(None, [("chrTest", "AC" * (BACKGROUND_FLOOR // 2))])
        assert under == UNIFORM_BACKGROUND
        assert exactly != UNIFORM_BACKGROUND

        # A run of N is not a composition, and must not count towards the floor.
        found, _ = resolve_background(None, [("chrTest", "N" * (2 * BACKGROUND_FLOOR))])
        assert found == UNIFORM_BACKGROUND

        # A scan discards soft-masking, so the composition it derives must not see it.
        masked, _ = resolve_background(None, [("chrTest", "aaaaccgt" * 2000)])
        assert masked == over

        # A peak set is thousands of short sequences; the floor is over the whole input.
        many, _ = resolve_background(None, records("AAAACCGT" * 10, 200))
        assert many == over


class TestTheModesAndTheExplicitForm:
    def test_the_three_modes_behave_as_named(self) -> None:
        assert BACKGROUND_MODES == ("auto", "uniform", "derive")
        sequences = [("chrTest", "AAAACCGT" * 2000)]
        assert resolve_background("uniform", sequences)[0] == UNIFORM_BACKGROUND
        derived = resolve_background("derive", [("chrTest", "AACC")])[0]
        assert derived == (0.375, 0.375, 0.125, 0.125)
        assert resolve_background(None, sequences)[0] == resolve_background("auto", sequences)[0]
        explicit = resolve_background([0.3, 0.2, 0.2, 0.3], sequences)[0]
        assert explicit == (0.3, 0.2, 0.2, 0.3)

    def test_an_invalid_background_names_the_problem(self) -> None:
        cases: list[tuple[object, str]] = [
            ("gc", "auto, uniform, derive"),
            ([0.25, 0.25, 0.25], "4 frequencies"),
            ([0.0, 0.25, 0.25, 0.5], "must all be > 0"),
            ([0.3, 0.3, 0.3, 0.3], "must sum to 1"),
        ]
        for background, message in cases:
            with pytest.raises(ValueError, match=message):
                resolve_background(background, [("chrTest", "ACGT")])  # type: ignore[arg-type]


class TestTheSourceSurvives:
    def test_the_source_is_preserved_and_drained_lazily(self) -> None:
        given_records = records("ACGTACGT" * 200, 20)
        _, rest = resolve_background(None, iter(given_records))
        assert list(rest) == given_records

        pulled: list[str] = []

        def source() -> Iterator[tuple[str, str]]:
            for name, bases in given_records:
                pulled.append(name)
                yield name, bases

        _, rest2 = resolve_background(None, source())
        assert len(pulled) < len(given_records)  # deciding stopped as soon as it had enough
        assert list(rest2) == given_records
        assert pulled == [name for name, _ in given_records]  # each record yielded once

        # A file of nothing but N is not buffered whole: the floor counts unambiguous
        # bases, so without a cap on how much is examined an all-N file would be pulled
        # into memory entire — which no scan may do.
        record_length = 10_000
        pulled_n: list[str] = []

        def n_source() -> Iterator[tuple[str, str]]:
            for index in range(1000):
                pulled_n.append(f"record{index}")
                yield f"record{index}", "N" * record_length

        found, rest3 = resolve_background(None, n_source())
        assert found == UNIFORM_BACKGROUND
        assert len(pulled_n) <= _SAMPLE_CAP // record_length + 1
        assert len(list(rest3)) == 1000

    def test_nothing_is_read_when_the_background_needs_no_input(self) -> None:
        cases: list[tuple[object, tuple[float, ...]]] = [
            ("uniform", UNIFORM_BACKGROUND),
            ([0.3, 0.2, 0.2, 0.3], (0.3, 0.2, 0.2, 0.3)),
        ]
        for background, expected in cases:

            def source() -> Iterator[tuple[str, str]]:
                raise AssertionError("the source was read")
                yield  # pragma: no cover - unreachable; it is here to make this a generator

            found, rest = resolve_background(background, source())  # type: ignore[arg-type]
            assert found == expected
            with pytest.raises(AssertionError, match="the source was read"):
                list(rest)


class TestQuantisation:
    def test_quantisation_rounds_onto_the_thousandth_grid(self) -> None:
        assert quantise_background([0.3, 0.2, 0.2, 0.3]) == (0.3, 0.2, 0.2, 0.3)
        assert quantise_background([0.25, 0.25, 0.25, 0.25]) == UNIFORM_BACKGROUND
        assert quantise_background([1 / 3, 1 / 3, 1 / 6, 1 / 6]) == (0.333, 0.333, 0.167, 0.167)
        # What lets two peak sets from one genome share a threshold cache entry.
        assert quantise_background([0.3001, 0.1999, 0.2, 0.3]) == quantise_background(
            [0.2999, 0.2001, 0.2, 0.3]
        )

    @given(background=frequencies)
    def test_quantisation_is_a_stable_idempotent_background(self, background: list[float]) -> None:
        found = quantise_background(background)
        assert len(found) == 4
        assert all(value > 0 for value in found)
        assert sum(found) == pytest.approx(1.0)
        for value in found:
            assert round(value * 1000) == pytest.approx(value * 1000)
        assert quantise_background(found) == found  # idempotent
        # The floor is set where the standard error on a frequency is about 0.004; rounding
        # onto the grid has to stay inside that, or it would be the larger perturbation.
        for wanted, got in zip(background, found, strict=True):
            assert abs(got - wanted) < 0.004


class TestDerivation:
    def test_counts_become_frequencies_with_a_pseudocount_for_safety(self) -> None:
        assert derive_background([250, 250, 250, 250]) == UNIFORM_BACKGROUND
        # A base never seen still gets a positive frequency: a frequency of zero would make
        # every log-odds against that base infinite.
        found = derive_background([1000, 0, 0, 1000])
        assert all(value > 0 for value in found)
        assert sum(found) == pytest.approx(1.0)
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
