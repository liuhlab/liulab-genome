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
    # Thousands separators and surrounding whitespace are tolerated too.
    assert parse_region("chr1:1,000-2,000") == ("chr1", 1000, 2000)
    assert parse_region("  chr1:0-10\n") == ("chr1", 0, 10)


def test_parse_region_bare_chromosome_has_no_coords() -> None:
    # One representative each: a plain name, an accession with a dot, a bare digit
    # that could be mistaken for a coordinate, and an underscore-separated scaffold.
    for text in ("chrM", "GL000009.2", "2", "scaffold_17"):
        assert parse_region(text) == (text, None, None)


def test_parse_region_malformed_locus_raises() -> None:
    # One representative each: no separator, missing end, missing start, wrong
    # separator, and a non-integer coordinate.
    for bad in ("chr1:abc", "chr1:0-", "chr1:-5", "chr1:0:10", "chr1:1.5-2"):
        with pytest.raises(ValueError, match="malformed region"):
            parse_region(bad)


@given(chrom=_chrom, start=_coord, span=st.integers(min_value=0, max_value=10**6))
def test_parse_region_roundtrips(chrom: str, start: int, span: int) -> None:
    end = start + span
    assert parse_region(f"{chrom}:{start}-{end}") == (chrom, start, end)


# --- Region ---


def test_region_basic_properties() -> None:
    unstranded = Region("chr1", 0, 10)
    assert (unstranded.chrom, unstranded.start, unstranded.end, unstranded.strand) == (
        "chr1",
        0,
        10,
        ".",
    )
    r = Region("chr1", 10, 25)
    assert len(r) == r.length == 15
    assert len(Region("chr1", 5, 5)) == 0  # an empty half-open interval is valid
    assert str(Region("chr2", 100, 200)) == "chr2:100-200"  # 0-based in, 0-based out


def test_region_from_string() -> None:
    assert Region.from_string("chr2:100-200", strand="-") == Region("chr2", 100, 200, "-")
    with pytest.raises(ValueError, match="no coordinates"):
        Region.from_string("chrM")


def test_region_construction_validates_coordinates_and_strand() -> None:
    for start, end in ((-1, 5), (5, 3)):
        with pytest.raises(ValueError, match="must be >="):
            Region("chr1", start, end)
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
