"""What a gene symbol is here, and why the two directions are not mirror images.

A pure module — it opens nothing, downloads nothing and knows no **Xref source** — and a
peer of :mod:`genome.xref.ids` and :mod:`genome.xref.evidence`: that one says how an
identifier is spelled, the second what may be asked about the evidence behind a pair, and
this one what a **Symbol match** is.

**A symbol is not an id, and it is not answered like one.** Away from the hub a **Gene id
stem** yields the authority's single current approved symbol, which is labelling a plot and
is one-to-one by the authority's own construction. Toward the hub a symbol matches
approved, previous **and** alias spellings, answers with every stem any of them names, and
records which kind of match each was — so ambiguity is the return type rather than an edge
case. The two verbs are therefore :meth:`~genome.xref.xref.XrefSet.from_stems` for the
first and :meth:`~genome.xref.xref.XrefSet.match_symbols` for the second, and
:meth:`~genome.xref.xref.XrefSet.to_stems` refuses the symbol **Namespace** rather than
answering it half-way: it would match approved spellings only, which is the failure this
whole thing exists to prevent — 31 of 801 EpiFactors rows spell their gene the way HGNC
spelled it years ago, and an approved-only match drops exactly those.

**Matching is exact by default.** The species is fixed by the set, so ``Brca1`` asked of a
human set is a mouse spelling asked of the wrong authority and says so by matching nothing.
Case-insensitive matching is opt-in, and the insensitive path still answers with **every**
gene matched rather than picking one.

**The kind rides on the match.** A source that publishes typed previous and alias spellings
carries all three kinds; one that publishes an undifferentiated synonym list carries
:data:`APPROVED` alone, because labelling an untyped spelling ``previous`` or ``alias``
would put a kind on a claim no publisher made. Which kinds a set carries, and why the
others are missing, ride back on every answer.

Examples
--------
>>> from genome.xref.symbols import SYMBOL_KINDS, fold_symbol, normalise_symbol
>>> SYMBOL_KINDS
('approved', 'previous', 'alias')
>>> normalise_symbol("  ARNTL  ")
'ARNTL'
>>> fold_symbol("Brca1") == fold_symbol("BRCA1")
True
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType

#: The authority's current spelling — the one it would print on a figure axis today, and
#: the only one every **Xref source** that carries symbols at all is able to publish.
APPROVED = "approved"

#: A spelling the authority itself used and has since retired. The kind that matters: a
#: table five years old spells 31 of EpiFactors' 801 rows this way and nothing else finds
#: them.
PREVIOUS = "previous"

#: A spelling the authority records as a name the gene also goes by, and never approved.
ALIAS = "alias"

#: The three kinds, most authoritative first, which is the order a match list is sorted in
#: so that an approved hit reads before a previous one for the same symbol.
SYMBOL_KINDS: tuple[str, ...] = (APPROVED, PREVIOUS, ALIAS)

#: Which **Namespace** each kind is stored under in the set's plain gzipped TSV. The kind
#: is carried by the namespace column rather than by a fourth column of its own: a column
#: that says *what sort of row this is* is the level discriminator this design refuses
#: everywhere else, and these three spellings read for themselves in a shell.
#:
#: Only :data:`APPROVED`'s spelling is a **Namespace** a caller may name — a previous and an
#: alias spelling are the same identifier system as an approved one with a different
#: standing, not systems of their own — so the other two are stored names and nothing more.
KIND_NAMESPACES: Mapping[str, str] = MappingProxyType(
    {APPROVED: "symbol", PREVIOUS: "previous_symbol", ALIAS: "alias_symbol"}
)

#: The three stored spellings, in :data:`SYMBOL_KINDS` order.
SYMBOL_NAMESPACES: tuple[str, ...] = tuple(KIND_NAMESPACES[kind] for kind in SYMBOL_KINDS)

#: The kind each stored **Namespace** carries — :data:`KIND_NAMESPACES` read backwards.
NAMESPACE_KINDS: Mapping[str, str] = MappingProxyType(
    {namespace: kind for kind, namespace in KIND_NAMESPACES.items()}
)


class SymbolDirectionError(ValueError):
    """A symbol was asked for the way an id is, and the two directions are not the same.

    Raised by :meth:`~genome.xref.xref.XrefSet.to_stems` when the **Namespace** named is
    the symbol one. A :class:`ValueError` rather than a
    :class:`~genome.xref.xref.NamespaceNotCarriedError`: the set does carry symbols, and
    what is wrong is the verb. Answering it anyway would match approved spellings only and
    silently drop every row that spells its gene the way the authority used to, which is
    the measured failure — so the message names
    :meth:`~genome.xref.xref.XrefSet.match_symbols` instead.

    Examples
    --------
    >>> human = XrefSet("Homo sapiens", "hgnc")                     # doctest: +SKIP
    >>> human.to_stems(["ARNTL"], "symbol")                         # doctest: +SKIP
    Traceback (most recent call last):
    SymbolDirectionError: ...
    """


def normalise_symbol(symbol: str) -> str:
    r"""Return ``symbol`` in the one spelling a match is looked up under.

    Surrounding whitespace and nothing else. A symbol is **not** put through
    :func:`~genome.xref.ids.normalise_id`: that drops everything after the first ``.``,
    which is right for a versioned gene id and wrong for a symbol — WormBase names
    thousands of genes by their sequence name, and ``Y110A7A.10`` stemmed to ``Y110A7A``
    is a different gene's spelling or none at all.

    Idempotent, and total: an empty or all-whitespace symbol comes back empty and matches
    nothing.

    Parameters
    ----------
    symbol : str
        The symbol, as a caller typed it or a publisher wrote it.

    Returns
    -------
    str
        The symbol with surrounding whitespace removed. Case is untouched, because case is
        significant: exact matching is the default and a mouse-cased spelling asked of a
        human set is the wrong authority's, not a typo to absorb.

    Examples
    --------
    >>> normalise_symbol("  TP53\n")
    'TP53'
    >>> normalise_symbol("Y110A7A.10")
    'Y110A7A.10'
    >>> normalise_symbol(normalise_symbol(" daf-16 "))
    'daf-16'
    """
    return symbol.strip()


def fold_symbol(symbol: str) -> str:
    """Return the key ``symbol`` is matched under when case is not to be significant.

    :func:`normalise_symbol` and then :meth:`str.casefold`, which is the full-Unicode fold
    rather than :meth:`str.lower`. Used only on the opt-in path, and used on **both** sides
    of it — the caller's spelling and the authority's — so that folding is a property of
    the lookup rather than of what was stored.

    Parameters
    ----------
    symbol : str
        The symbol, in any spelling.

    Returns
    -------
    str
        Its case-folded key. Empty for an empty or all-whitespace symbol.

    Examples
    --------
    >>> fold_symbol("Brca1")
    'brca1'
    >>> fold_symbol(" p53 ") == fold_symbol("P53")
    True
    """
    return normalise_symbol(symbol).casefold()
