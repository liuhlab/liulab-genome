"""Tests for genome.region — the 0-based, half-open coordinate primitives.

No native binaries are needed: this is pure coordinate logic. Hypothesis pins
the invariants (parse round-trips, ``len`` == span, the 0→1-based boundary).
"""

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.region import Region, parse_region

# Chromosome-name strategy: anything the parser treats as a name, minus the ``:``
# that separates name from coordinates.
_chrom = st.text(
    alphabet="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-",
    min_size=1,
    max_size=12,
)
_coord = st.integers(min_value=0, max_value=10**9)


# --- parse_region ---


def test_parse_region_basic() -> None:
    assert parse_region("chr1:0-10") == ("chr1", 0, 10)


def test_parse_region_tolerates_thousands_separators() -> None:
    assert parse_region("chr1:1,000-2,000") == ("chr1", 1000, 2000)


def test_parse_region_strips_surrounding_whitespace() -> None:
    assert parse_region("  chr1:0-10\n") == ("chr1", 0, 10)


@pytest.mark.parametrize("text", ["chrM", "GL000009.2", "2", "scaffold_17"])
def test_parse_region_bare_chromosome_has_no_coords(text: str) -> None:
    assert parse_region(text) == (text, None, None)


@pytest.mark.parametrize("bad", ["chr1:abc", "chr1:0-", "chr1:-5", "chr1:0:10", "chr1:1.5-2"])
def test_parse_region_malformed_locus_raises(bad: str) -> None:
    with pytest.raises(ValueError, match="malformed region"):
        parse_region(bad)


@given(chrom=_chrom, start=_coord, span=st.integers(min_value=0, max_value=10**6))
def test_parse_region_roundtrips(chrom: str, start: int, span: int) -> None:
    end = start + span
    assert parse_region(f"{chrom}:{start}-{end}") == (chrom, start, end)


# --- Region ---


def test_region_defaults_to_unknown_strand() -> None:
    r = Region("chr1", 0, 10)
    assert (r.chrom, r.start, r.end, r.strand) == ("chr1", 0, 10, ".")


def test_region_len_and_length_agree() -> None:
    r = Region("chr1", 10, 25)
    assert len(r) == r.length == 15


def test_region_empty_interval_is_valid() -> None:
    assert len(Region("chr1", 5, 5)) == 0


def test_region_str_is_zero_based() -> None:
    assert str(Region("chr2", 100, 200)) == "chr2:100-200"


def test_region_from_string_attaches_strand() -> None:
    assert Region.from_string("chr2:100-200", strand="-") == Region("chr2", 100, 200, "-")


def test_region_from_string_requires_coordinates() -> None:
    with pytest.raises(ValueError, match="no coordinates"):
        Region.from_string("chrM")


@pytest.mark.parametrize(("start", "end"), [(-1, 5), (5, 3)])
def test_region_invalid_coordinates_raise(start: int, end: int) -> None:
    with pytest.raises(ValueError, match="must be >="):
        Region("chr1", start, end)


def test_region_invalid_strand_raises() -> None:
    with pytest.raises(ValueError, match="strand"):
        Region("chr1", 0, 10, strand="x")


@given(chrom=_chrom, start=_coord, span=st.integers(min_value=0, max_value=10**6))
def test_region_len_equals_span(chrom: str, start: int, span: int) -> None:
    assert len(Region(chrom, start, start + span)) == span


@given(
    chrom=_chrom,
    start=_coord,
    span=st.integers(min_value=0, max_value=10**6),
    strand=st.sampled_from(["+", "-", "."]),
)
def test_region_string_roundtrips(chrom: str, start: int, span: int, strand: str) -> None:
    r = Region(chrom, start, start + span, strand)
    # str(r) drops the strand; from_string re-attaches the same one.
    assert Region.from_string(str(r), strand=strand) == r
