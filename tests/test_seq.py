"""Tests for genome.seq — concrete examples plus hypothesis property tests."""

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.seq import DNA, RNA, _Seq

# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

dna_text = st.text(alphabet="ACGTacgt", min_size=0, max_size=64)
rna_text = st.text(alphabet="ACGUacgu", min_size=0, max_size=64)


# ---------------------------------------------------------------------------
# Construction (alphabet is NOT enforced — see seq.py)
# ---------------------------------------------------------------------------


class TestConstruction:
    def test_construction_preserves_case_including_non_alphabet_chars_and_empty(self) -> None:
        assert DNA("") == ""
        assert RNA("") == ""
        assert str(DNA("aTcG")) == "aTcG"
        assert str(RNA("aUcG")) == "aUcG"
        # Alphabet checking is intentionally skipped for performance; any
        # string is accepted and stored as-is.
        assert str(DNA("ATCX")) == "ATCX"
        assert str(DNA("AUCG")) == "AUCG"
        assert str(RNA("ATCG")) == "ATCG"
        with pytest.raises(TypeError, match="abstract"):
            _Seq("ATCG")


# ---------------------------------------------------------------------------
# The alphabet, and the check a caller at the I/O boundary asks for
# ---------------------------------------------------------------------------


class TestOutsideAlphabet:
    def test_outside_alphabet_reports_offenders_case_insensitively_sorted_and_distinct(
        self,
    ) -> None:
        assert DNA.outside_alphabet("ATCX") == ["X"]
        assert DNA.outside_alphabet("ATCG") == []
        assert DNA.outside_alphabet("") == []
        assert DNA.outside_alphabet("aTcG") == []
        # Offenders come back in their own case, sorted and de-duplicated.
        assert DNA.outside_alphabet("aTcx") == ["x"]
        assert DNA.outside_alphabet("XBXB") == ["B", "X"]

    def test_each_class_asks_its_own_alphabet(self) -> None:
        # The classmethod lives on the shared base; the alphabet it compares against is
        # the calling class's, so U offends DNA and T offends RNA.
        assert DNA.outside_alphabet("AUCG") == ["U"]
        assert RNA.outside_alphabet("ATCG") == ["T"]
        assert RNA.outside_alphabet("AUCG") == []
        # Joined sorted, which is how a caller naming the alphabet in a message renders it.
        assert "".join(sorted(DNA.ALPHABET)) == "ACGT"
        assert "".join(sorted(RNA.ALPHABET)) == "ACGU"


class TestPropertiesOutsideAlphabet:
    @given(st.one_of(dna_text, st.text(max_size=32)))
    def test_offenders_are_a_sorted_distinct_subset_and_empty_iff_all_in_alphabet(
        self, s: str
    ) -> None:
        offenders = DNA.outside_alphabet(s)
        assert offenders == sorted(set(offenders))
        assert set(offenders) <= set(s)
        assert (offenders == []) == all(c.upper() in DNA.ALPHABET for c in s)

    @given(st.sampled_from([DNA, RNA]), st.text(max_size=32))
    def test_every_subclass_answers_against_its_own_alphabet(self, cls: type[_Seq], s: str) -> None:
        assert cls.outside_alphabet(s) == sorted({c for c in s if c.upper() not in cls.ALPHABET})

    @given(st.text(max_size=32))
    def test_reporting_an_offender_never_blocks_construction(self, s: str) -> None:
        # The check is offered, never imposed: reporting an offender does not stop the
        # same string from constructing (see docs/adr/0005).
        DNA.outside_alphabet(s)
        assert str(DNA(s)) == s


# ---------------------------------------------------------------------------
# Typed slicing, repr, and interplay with str
# ---------------------------------------------------------------------------


class TestSlicingAndRepr:
    def test_typed_slicing_indexing_repr_and_inherited_str_methods(self) -> None:
        sliced = DNA("ATCG")[1:3]
        assert isinstance(sliced, DNA)
        assert sliced == "TC"
        indexed = DNA("ATCG")[0]
        assert isinstance(indexed, DNA)
        assert indexed == "A"
        assert repr(DNA("ATCG")) == "DNA('ATCG')"
        assert repr(RNA("AUCG")) == "RNA('AUCG')"
        # Documented contract: only __getitem__ and biological methods stay typed.
        upper_result = DNA("aTcG").upper()
        assert type(upper_result) is str
        assert upper_result == "ATCG"


# ---------------------------------------------------------------------------
# DNA/RNA biological transforms
# ---------------------------------------------------------------------------


class TestDNATransforms:
    def test_complement_reverse_complement_and_transcribe(self) -> None:
        assert DNA("ATCG").complement() == DNA("TAGC")
        assert DNA("aTcG").complement() == DNA("tAgC")
        assert DNA("ATCG").reverse_complement() == DNA("CGAT")
        assert DNA("aTcG").reverse_complement() == DNA("CgAt")
        assert DNA("").reverse_complement() == DNA("")
        assert DNA("ATCG").transcribe() == RNA("AUCG")
        assert isinstance(DNA("ATCG").transcribe(), RNA)
        assert DNA("aTcG").transcribe() == RNA("aUcG")


class TestRNATransforms:
    def test_complement_reverse_complement_and_back_transcribe(self) -> None:
        assert RNA("AUCG").complement() == RNA("UAGC")
        assert RNA("AUCG").reverse_complement() == RNA("CGAU")
        assert RNA("AUCG").back_transcribe() == DNA("ATCG")
        assert isinstance(RNA("AUCG").back_transcribe(), DNA)
        assert RNA("aUcG").back_transcribe() == DNA("aTcG")


# ---------------------------------------------------------------------------
# GC content
# ---------------------------------------------------------------------------


class TestGCContent:
    def test_gc_content_pure_mixed_and_empty_for_dna_and_rna(self) -> None:
        assert DNA("GGCC").gc_content == 1.0
        assert DNA("ATAT").gc_content == 0.0
        assert DNA("aTcG").gc_content == 0.5
        assert DNA("").gc_content == 0.0
        assert RNA("").gc_content == 0.0
        assert RNA("GGCC").gc_content == 1.0
        assert RNA("AUAU").gc_content == 0.0


# ---------------------------------------------------------------------------
# Hypothesis property tests
# ---------------------------------------------------------------------------


class TestPropertiesDNA:
    @given(st.text(min_size=0, max_size=32))
    def test_constructor_preserves_case_and_value_for_any_text(self, s: str) -> None:
        # Alphabet is not enforced: any string constructs and is stored as-is.
        assert str(DNA(s)) == s

    @given(dna_text)
    def test_complement_transcription_slicing_and_gc_content_invariants(self, s: str) -> None:
        d = DNA(s)
        assert len(d.complement()) == len(s)
        assert d.complement().complement() == d
        assert d.reverse_complement().reverse_complement() == d
        assert d.transcribe().back_transcribe() == d
        assert 0.0 <= d.gc_content <= 1.0
        for i in range(len(d) + 1):
            for j in range(i, len(d) + 1):
                piece = d[i:j]
                assert isinstance(piece, DNA)
                assert str(piece) == s[i:j]


class TestPropertiesRNA:
    @given(rna_text)
    def test_reverse_complement_is_involution_and_back_transcribe_roundtrips(self, s: str) -> None:
        r = RNA(s)
        assert r.reverse_complement().reverse_complement() == r
        assert r.back_transcribe().transcribe() == r
