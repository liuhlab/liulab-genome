"""Typed biological sequence classes.

This module defines :class:`DNA`, :class:`RNA`, and :class:`Protein` — thin,
typed subclasses of :class:`str` that carry biological transforms (complement,
reverse complement, transcription) without losing their type or case.

Notes
-----
- Construction does **not** validate the alphabet: scanning every character is
  prohibitively expensive on large sequences (whole chromosomes), so the
  subclass alphabet is documentation, not a runtime check. Validate at the I/O
  boundary if you need to reject non-alphabet input.
- The stored value is **preserved verbatim**, including case (lowercase carries
  meaning — soft-masking).
- Slicing and indexing return the **same subclass**, not a plain :class:`str`.
- Biological transforms return the correct typed class and preserve case.
- Other inherited :class:`str` methods (``upper``, ``lower``, ``replace``, …)
  return plain :class:`str` by design — if you need a typed variant, wrap
  the result explicitly: ``DNA(my_dna.upper())``.
- IUPAC ambiguity codes (``N`` and friends) are intentionally **out of scope**.

Examples
--------
>>> from genome import DNA, RNA
>>> DNA("ATCG").reverse_complement()
DNA('CGAT')
>>> DNA("aTcG").reverse_complement()        # case preserved
DNA('CgAt')
>>> DNA("ATCG").transcribe()
RNA('AUCG')
>>> DNA("GGCC").gc_content
1.0
"""

from __future__ import annotations

from typing import ClassVar, Self

# Module-level translation tables (built once; cheap reuse).
_DNA_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
_RNA_COMPLEMENT = str.maketrans("ACGUacgu", "UGCAugca")
_T_TO_U = str.maketrans("Tt", "Uu")
_U_TO_T = str.maketrans("Uu", "Tt")


class _Seq(str):
    """Private base for typed sequence strings.

    Subclasses declare their expected alphabet via the class variable
    :attr:`_ALPHABET` (uppercase characters only) for documentation; it is
    **not enforced** at construction (see :meth:`__new__`). The original value
    and case are preserved verbatim.

    This class is abstract — instantiating it directly raises :class:`TypeError`.
    """

    _ALPHABET: ClassVar[frozenset[str]] = frozenset()

    def __new__(cls, value: str) -> Self:
        """Return a typed instance wrapping ``value`` (case preserved).

        The subclass alphabet (:attr:`_ALPHABET`) documents the *expected*
        characters but is **not enforced**: alphabet checking is intentionally
        skipped because scanning every character is prohibitively expensive on
        large sequences (whole chromosomes). Callers that need to reject
        non-alphabet input must validate at the I/O boundary.

        Parameters
        ----------
        value
            The sequence string. May be empty. Stored verbatim.

        Returns
        -------
        Self
            An instance of the calling subclass, storing ``value`` verbatim
            (case preserved).

        Raises
        ------
        TypeError
            If called on :class:`_Seq` directly rather than a concrete subclass.
        """
        if cls is _Seq:
            raise TypeError("_Seq is abstract; instantiate DNA, RNA, or Protein instead.")
        return str.__new__(cls, value)

    @classmethod
    def _unchecked(cls, value: str) -> Self:
        """Build an instance without re-validating.

        Reserved for internal transforms whose output alphabet is guaranteed
        valid by construction (e.g. a complement of a valid DNA is valid DNA).
        Do not call from outside this module.
        """
        return str.__new__(cls, value)

    def __getitem__(self, key: int | slice) -> Self:  # type: ignore[override]
        """Return a slice/index as the same subclass, not plain ``str``."""
        return type(self)._unchecked(str.__getitem__(self, key))

    def __repr__(self) -> str:
        """Return e.g. ``DNA('ATCG')``."""
        return f"{type(self).__name__}({str.__repr__(self)})"


class DNA(_Seq):
    """A DNA sequence over the four canonical bases.

    Alphabet: ``A``, ``C``, ``G``, ``T`` (case-insensitive at construction;
    case preserved in the stored value).

    Examples
    --------
    >>> DNA("ATCG")
    DNA('ATCG')
    >>> DNA("aTcG")[1:3]                      # slicing stays typed
    DNA('Tc')
    >>> DNA("ATCG").reverse_complement()
    DNA('CGAT')
    >>> DNA("ATCG").transcribe()
    RNA('AUCG')
    >>> DNA("GGCC").gc_content
    1.0
    """

    _ALPHABET: ClassVar[frozenset[str]] = frozenset("ACGT")

    def complement(self) -> DNA:
        """Return the Watson-Crick complement, case preserved.

        Returns
        -------
        DNA
            A new :class:`DNA` of the same length with ``A↔T`` and ``C↔G``
            swapped. Lowercase letters complement to lowercase.

        Examples
        --------
        >>> DNA("ATCG").complement()
        DNA('TAGC')
        >>> DNA("aTcG").complement()
        DNA('tAgC')
        """
        return DNA._unchecked(self.translate(_DNA_COMPLEMENT))

    def reverse_complement(self) -> DNA:
        """Return the reverse complement, case preserved.

        Returns
        -------
        DNA
            A new :class:`DNA` of the same length: complement then reverse.

        Examples
        --------
        >>> DNA("ATCG").reverse_complement()
        DNA('CGAT')
        >>> DNA("aTcG").reverse_complement()
        DNA('CgAt')
        >>> DNA("").reverse_complement()
        DNA('')
        """
        return DNA._unchecked(self.translate(_DNA_COMPLEMENT)[::-1])

    def transcribe(self) -> RNA:
        """Transcribe DNA to RNA by replacing ``T``/``t`` with ``U``/``u``.

        Returns
        -------
        RNA
            A new :class:`RNA` of equal length; other bases (and case) unchanged.

        Examples
        --------
        >>> DNA("ATCG").transcribe()
        RNA('AUCG')
        >>> DNA("aTcG").transcribe()
        RNA('aUcG')
        """
        return RNA._unchecked(self.translate(_T_TO_U))

    @property
    def gc_content(self) -> float:
        """Fraction of bases that are ``G`` or ``C`` (case-insensitive).

        Returns
        -------
        float
            Value in ``[0.0, 1.0]``. Defined as ``0.0`` for the empty sequence.

        Examples
        --------
        >>> DNA("GGCC").gc_content
        1.0
        >>> DNA("ATAT").gc_content
        0.0
        >>> DNA("aTcG").gc_content
        0.5
        >>> DNA("").gc_content
        0.0
        """
        return _gc_fraction(self)


class RNA(_Seq):
    """An RNA sequence over the four canonical bases.

    Alphabet: ``A``, ``C``, ``G``, ``U`` (case-insensitive at construction;
    case preserved in the stored value).

    Examples
    --------
    >>> RNA("AUCG").reverse_complement()
    RNA('CGAU')
    >>> RNA("AUCG").back_transcribe()
    DNA('ATCG')
    """

    _ALPHABET: ClassVar[frozenset[str]] = frozenset("ACGU")

    def complement(self) -> RNA:
        """Return the complement, case preserved (``A↔U``, ``C↔G``).

        Examples
        --------
        >>> RNA("AUCG").complement()
        RNA('UAGC')
        """
        return RNA._unchecked(self.translate(_RNA_COMPLEMENT))

    def reverse_complement(self) -> RNA:
        """Return the reverse complement, case preserved.

        Examples
        --------
        >>> RNA("AUCG").reverse_complement()
        RNA('CGAU')
        """
        return RNA._unchecked(self.translate(_RNA_COMPLEMENT)[::-1])

    def back_transcribe(self) -> DNA:
        """Reverse-transcribe to DNA by replacing ``U``/``u`` with ``T``/``t``.

        Examples
        --------
        >>> RNA("AUCG").back_transcribe()
        DNA('ATCG')
        >>> RNA("aUcG").back_transcribe()
        DNA('aTcG')
        """
        return DNA._unchecked(self.translate(_U_TO_T))

    @property
    def gc_content(self) -> float:
        """Fraction of bases that are ``G`` or ``C`` (case-insensitive).

        ``0.0`` for the empty sequence.

        Examples
        --------
        >>> RNA("GGCC").gc_content
        1.0
        >>> RNA("").gc_content
        0.0
        """
        return _gc_fraction(self)


class Protein(_Seq):
    """A protein sequence over the 20 standard amino acids.

    Alphabet: ``ACDEFGHIKLMNPQRSTVWY`` (case-insensitive at construction;
    case preserved in the stored value). No biological transforms are defined.

    Examples
    --------
    >>> Protein("MKTAY")
    Protein('MKTAY')
    >>> Protein("MKTAY")[1:3]
    Protein('KT')
    """

    _ALPHABET: ClassVar[frozenset[str]] = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _gc_fraction(seq: str) -> float:
    """Return the GC fraction of ``seq`` (case-insensitive); ``0.0`` if empty."""
    if not seq:
        return 0.0
    upper = seq.upper()
    return (upper.count("G") + upper.count("C")) / len(seq)
