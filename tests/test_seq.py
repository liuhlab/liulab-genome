"""Tests for genome.seq — concrete examples plus hypothesis property tests."""

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.seq import DNA, RNA, Protein, _Seq

# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

dna_text = st.text(alphabet="ACGTacgt", min_size=0, max_size=64)
rna_text = st.text(alphabet="ACGUacgu", min_size=0, max_size=64)
protein_text = st.text(
    alphabet="ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy",
    min_size=0,
    max_size=64,
)


# ---------------------------------------------------------------------------
# Construction & validation
# ---------------------------------------------------------------------------


class TestConstruction:
    def test_empty_is_valid(self) -> None:
        assert DNA("") == ""
        assert RNA("") == ""
        assert Protein("") == ""

    def test_basic_construction_preserves_case(self) -> None:
        assert str(DNA("aTcG")) == "aTcG"
        assert str(RNA("aUcG")) == "aUcG"
        assert str(Protein("mKtAy")) == "mKtAy"

    def test_dna_rejects_non_alphabet(self) -> None:
        with pytest.raises(ValueError, match="X"):
            DNA("ATCX")

    def test_dna_rejects_rna_letter_u(self) -> None:
        with pytest.raises(ValueError, match="U"):
            DNA("AUCG")

    def test_rna_rejects_dna_letter_t(self) -> None:
        with pytest.raises(ValueError, match="T"):
            RNA("ATCG")

    def test_protein_rejects_b_o_u_z(self) -> None:
        # B, O, U, Z are not in the 20-standard alphabet
        for ch in "BOUZ":
            with pytest.raises(ValueError, match=ch):
                Protein(f"M{ch}T")

    def test_validation_lists_all_offending_chars_sorted(self) -> None:
        with pytest.raises(ValueError, match="outside alphabet") as exc:
            DNA("ATCXZ_")
        msg = str(exc.value)
        for ch in ("X", "Z", "_"):
            assert ch in msg

    def test_seq_base_is_abstract(self) -> None:
        with pytest.raises(TypeError, match="abstract"):
            _Seq("ATCG")


# ---------------------------------------------------------------------------
# Typed slicing, repr, and interplay with str
# ---------------------------------------------------------------------------


class TestSlicingAndRepr:
    def test_slice_returns_same_subclass(self) -> None:
        s = DNA("ATCG")[1:3]
        assert isinstance(s, DNA)
        assert s == "TC"

    def test_index_returns_same_subclass(self) -> None:
        s = DNA("ATCG")[0]
        assert isinstance(s, DNA)
        assert s == "A"

    def test_repr_shape(self) -> None:
        assert repr(DNA("ATCG")) == "DNA('ATCG')"
        assert repr(RNA("AUCG")) == "RNA('AUCG')"
        assert repr(Protein("MKT")) == "Protein('MKT')"

    def test_inherited_str_methods_return_plain_str(self) -> None:
        # Documented contract: only __getitem__ and biological methods stay typed.
        upper_result = DNA("aTcG").upper()
        assert type(upper_result) is str
        assert upper_result == "ATCG"


# ---------------------------------------------------------------------------
# DNA biological transforms
# ---------------------------------------------------------------------------


class TestDNATransforms:
    def test_complement_uppercase(self) -> None:
        assert DNA("ATCG").complement() == DNA("TAGC")

    def test_complement_preserves_case(self) -> None:
        assert DNA("aTcG").complement() == DNA("tAgC")

    def test_reverse_complement_uppercase(self) -> None:
        assert DNA("ATCG").reverse_complement() == DNA("CGAT")

    def test_reverse_complement_preserves_case(self) -> None:
        assert DNA("aTcG").reverse_complement() == DNA("CgAt")

    def test_reverse_complement_empty(self) -> None:
        assert DNA("").reverse_complement() == DNA("")

    def test_transcribe(self) -> None:
        assert DNA("ATCG").transcribe() == RNA("AUCG")
        assert isinstance(DNA("ATCG").transcribe(), RNA)

    def test_transcribe_preserves_case(self) -> None:
        assert DNA("aTcG").transcribe() == RNA("aUcG")


class TestRNATransforms:
    def test_complement(self) -> None:
        assert RNA("AUCG").complement() == RNA("UAGC")

    def test_reverse_complement(self) -> None:
        assert RNA("AUCG").reverse_complement() == RNA("CGAU")

    def test_back_transcribe(self) -> None:
        assert RNA("AUCG").back_transcribe() == DNA("ATCG")
        assert isinstance(RNA("AUCG").back_transcribe(), DNA)

    def test_back_transcribe_preserves_case(self) -> None:
        assert RNA("aUcG").back_transcribe() == DNA("aTcG")


# ---------------------------------------------------------------------------
# GC content
# ---------------------------------------------------------------------------


class TestGCContent:
    def test_pure_gc(self) -> None:
        assert DNA("GGCC").gc_content == 1.0

    def test_pure_at(self) -> None:
        assert DNA("ATAT").gc_content == 0.0

    def test_mixed(self) -> None:
        assert DNA("aTcG").gc_content == 0.5

    def test_empty(self) -> None:
        assert DNA("").gc_content == 0.0
        assert RNA("").gc_content == 0.0

    def test_rna(self) -> None:
        assert RNA("GGCC").gc_content == 1.0
        assert RNA("AUAU").gc_content == 0.0


# ---------------------------------------------------------------------------
# Hypothesis property tests
# ---------------------------------------------------------------------------


class TestPropertiesDNA:
    @given(dna_text)
    def test_constructor_preserves_case_and_value(self, s: str) -> None:
        assert str(DNA(s)) == s

    @given(dna_text)
    def test_length_preserved_by_complement(self, s: str) -> None:
        assert len(DNA(s).complement()) == len(s)

    @given(dna_text)
    def test_complement_is_involution(self, s: str) -> None:
        d = DNA(s)
        assert d.complement().complement() == d

    @given(dna_text)
    def test_reverse_complement_is_involution(self, s: str) -> None:
        d = DNA(s)
        assert d.reverse_complement().reverse_complement() == d

    @given(dna_text)
    def test_transcribe_back_transcribe_roundtrip(self, s: str) -> None:
        d = DNA(s)
        assert d.transcribe().back_transcribe() == d

    @given(dna_text)
    def test_slice_returns_dna_and_stays_valid(self, s: str) -> None:
        d = DNA(s)
        for i in range(len(d) + 1):
            for j in range(i, len(d) + 1):
                piece = d[i:j]
                assert isinstance(piece, DNA)
                assert str(piece) == s[i:j]

    @given(dna_text)
    def test_gc_content_in_unit_interval(self, s: str) -> None:
        assert 0.0 <= DNA(s).gc_content <= 1.0

    @given(st.text(min_size=1, max_size=32))
    def test_rejects_non_alphabet(self, s: str) -> None:
        # If every char is a valid DNA letter, construction must succeed; else fail.
        all_valid = all(c.upper() in "ACGT" for c in s)
        if all_valid:
            DNA(s)
        else:
            with pytest.raises(ValueError, match="outside alphabet"):
                DNA(s)


class TestPropertiesRNA:
    @given(rna_text)
    def test_reverse_complement_is_involution(self, s: str) -> None:
        r = RNA(s)
        assert r.reverse_complement().reverse_complement() == r

    @given(rna_text)
    def test_back_transcribe_transcribe_roundtrip(self, s: str) -> None:
        r = RNA(s)
        assert r.back_transcribe().transcribe() == r


class TestPropertiesProtein:
    @given(protein_text)
    def test_constructor_preserves_case_and_value(self, s: str) -> None:
        assert str(Protein(s)) == s

    @given(protein_text)
    def test_slice_returns_protein(self, s: str) -> None:
        p = Protein(s)
        if len(p) >= 2:
            assert isinstance(p[1:], Protein)
