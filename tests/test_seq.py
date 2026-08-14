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
# Construction (alphabet is NOT enforced — see seq.py)
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

    def test_non_alphabet_chars_are_accepted_verbatim(self) -> None:
        # Alphabet checking is intentionally skipped for performance; any
        # string is accepted and stored as-is.
        assert str(DNA("ATCX")) == "ATCX"
        assert str(DNA("AUCG")) == "AUCG"
        assert str(RNA("ATCG")) == "ATCG"
        assert str(Protein("MBTOUZ")) == "MBTOUZ"

    def test_seq_base_is_abstract(self) -> None:
        with pytest.raises(TypeError, match="abstract"):
            _Seq("ATCG")


# ---------------------------------------------------------------------------
# The alphabet, and the check a caller at the I/O boundary asks for
# ---------------------------------------------------------------------------


class TestOutsideAlphabet:
    def test_reports_the_offending_character(self) -> None:
        assert DNA.outside_alphabet("ATCX") == ["X"]

    def test_clean_sequence_reports_nothing(self) -> None:
        assert DNA.outside_alphabet("ATCG") == []
        assert DNA.outside_alphabet("") == []

    def test_case_insensitive(self) -> None:
        assert DNA.outside_alphabet("aTcG") == []

    def test_offenders_come_back_in_their_own_case(self) -> None:
        assert DNA.outside_alphabet("aTcx") == ["x"]

    def test_distinct_and_sorted(self) -> None:
        assert DNA.outside_alphabet("XBXB") == ["B", "X"]

    def test_each_class_asks_its_own_alphabet(self) -> None:
        # The classmethod lives on the shared base; the alphabet it compares against is
        # the calling class's, so U offends DNA and T offends RNA.
        assert DNA.outside_alphabet("AUCG") == ["U"]
        assert RNA.outside_alphabet("ATCG") == ["T"]
        assert RNA.outside_alphabet("AUCG") == []
        assert Protein.outside_alphabet("MKTAY") == []
        assert Protein.outside_alphabet("MBTOUZ") == ["B", "O", "U", "Z"]

    def test_alphabet_is_the_class_s_own(self) -> None:
        # Joined sorted, which is how a caller naming the alphabet in a message renders it.
        assert "".join(sorted(DNA.ALPHABET)) == "ACGT"
        assert "".join(sorted(RNA.ALPHABET)) == "ACGU"
        assert "".join(sorted(Protein.ALPHABET)) == "ACDEFGHIKLMNPQRSTVWY"

    def test_construction_still_validates_nothing(self) -> None:
        # The check is offered, never imposed: reporting an offender does not stop the
        # same string from constructing (see docs/adr/0005).
        assert DNA.outside_alphabet("ATCX") == ["X"]
        assert str(DNA("ATCX")) == "ATCX"


class TestPropertiesOutsideAlphabet:
    @given(dna_text)
    def test_in_alphabet_text_reports_nothing(self, s: str) -> None:
        assert DNA.outside_alphabet(s) == []

    @given(st.text(max_size=32))
    def test_offenders_are_a_sorted_distinct_subset_of_the_input(self, s: str) -> None:
        offenders = DNA.outside_alphabet(s)
        assert offenders == sorted(set(offenders))
        assert set(offenders) <= set(s)

    @given(st.text(max_size=32))
    def test_empty_exactly_when_every_character_is_in_the_alphabet(self, s: str) -> None:
        assert (DNA.outside_alphabet(s) == []) == all(c.upper() in DNA.ALPHABET for c in s)

    @given(st.sampled_from([DNA, RNA, Protein]), st.text(max_size=32))
    def test_every_subclass_answers_against_its_own_alphabet(self, cls: type[_Seq], s: str) -> None:
        assert cls.outside_alphabet(s) == sorted({c for c in s if c.upper() not in cls.ALPHABET})

    @given(st.text(max_size=32))
    def test_reporting_an_offender_never_blocks_construction(self, s: str) -> None:
        DNA.outside_alphabet(s)
        assert str(DNA(s)) == s


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
    def test_construction_accepts_any_text_verbatim(self, s: str) -> None:
        # Alphabet is not enforced: any string constructs and is stored as-is.
        assert str(DNA(s)) == s


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
