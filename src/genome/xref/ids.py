"""What an identifier is here: which **Namespace** it belongs to, and how it is spelled.

A pure module — it opens nothing, downloads nothing and knows no **Xref source**. Every
identifier that reaches an **Xref set**, from a publisher's file or from a caller's list,
comes through :func:`normalise_id` first, and every hub id through
:func:`gene_id_stem`, so the two sides of every join are spelled the same way.

**This is the most error-prone detail in the landscape.** Publishers spell versions
inconsistently and even inconsistently within one row — NCBI's ``gene2ensembl`` writes the
gene id bare and the transcript id versioned — and joining a versioned id to a bare one
returns zero matches with no error at all. So a version suffix is dropped on ingest, on
both sides, always: ``ENSG00000141510.18`` and ``ENSG00000141510`` are one identifier
here. Worm ids are unversioned as published, so a **Gene id stem** and a WormBase gene id
are the same string and the reduction is the identity.

The second inconsistency is the CURIE prefix. One **Namespace**'s ids are conventionally
written with theirs and another's without — nobody writes an HGNC id as ``1100`` and
nobody writes an Ensembl gene as ``ENSEMBL:ENSG00000141510`` — so each namespace declares
the prefix its canonical spelling carries (:data:`NAMESPACE_PREFIX`) and the prefixes it
is *found* under in published files (:data:`_NAMESPACE_ALIASES`). Both are accepted on the
way in and one comes back out, which is what lets two **Xref source**s be compared at all:
an answer from Alliance and an answer from Ensembl are two answers (ADR-0017), and they
would not even be legible side by side if each spelled the same id its publisher's way.

Examples
--------
>>> from genome.xref.ids import gene_id_stem, normalise_id
>>> gene_id_stem("ENSG00000141510.18")
'ENSG00000141510'
>>> normalise_id("1100", "hgnc")
'HGNC:1100'
>>> normalise_id("NCBI_Gene:672", "entrez")
'672'
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from types import MappingProxyType

#: The hub every other **Namespace** is a spoke off (ADR-0017). Its ids reduce to **Gene
#: id stem**s, and its string shape is per-species: ``ENSG…``, ``ENSMUSG…``, and for worm
#: the WormBase gene id, which *is* the Ensembl stable gene id.
ENSEMBL = "ensembl"

#: NCBI's Entrez GeneID — a bare integer, which is how a GEO series spells a gene.
ENTREZ = "entrez"

#: A UniProtKB accession, which is what a mass-spec run hands back.
UNIPROT = "uniprot"

#: The three species authorities. One per species and never interchangeable: an **Xref
#: set** carries the one its own species has and raises for the other two.
HGNC, MGI, WORMBASE = "hgnc", "mgi", "wormbase"

#: The gene symbol — what a paper's supplementary table actually holds. A **Namespace**
#: like the rest and answered unlike the rest: away from the hub it is the authority's one
#: current approved spelling, and toward the hub it is
#: :meth:`~genome.xref.xref.XrefSet.match_symbols` rather than the ordinary verb, because a
#: symbol matches previous and alias spellings too and each match carries which kind it
#: was. See :mod:`genome.xref.symbols`.
SYMBOL = "symbol"

#: Every **Namespace** this package knows how to spell, in the order an answer lists them.
#: A **Namespace** an **Xref source** does not carry raises naming the ones it does, and
#: what a set carries is read off the set rather than from here — this is the vocabulary,
#: not a promise about any one source.
NAMESPACES: tuple[str, ...] = (ENSEMBL, ENTREZ, UNIPROT, HGNC, MGI, WORMBASE, SYMBOL)

#: The prefix each **Namespace**'s canonical spelling carries — the way a practitioner
#: writes the id, not the way a CURIE would. HGNC's and MGI's ids are written with theirs
#: and nobody writes them without; the other five are written bare and nobody writes them
#: with. There is no rule behind this, only usage, which is why it is a table.
#:
#: :data:`SYMBOL` is here for completeness and is never reached: a symbol goes through
#: :func:`~genome.xref.symbols.normalise_symbol` instead, since :func:`gene_id_stem` would
#: cut ``Y110A7A.10`` down to a spelling WormBase gives no gene.
NAMESPACE_PREFIX: Mapping[str, str] = MappingProxyType(
    {
        ENSEMBL: "",
        ENTREZ: "",
        UNIPROT: "",
        HGNC: "HGNC:",
        MGI: "MGI:",
        WORMBASE: "",
        SYMBOL: "",
    }
)

#: What separates an id from its version, everywhere in this package. The annotation half
#: splits on the same character, which is what makes a stem produced here joinable to the
#: gene ids one **Annotation** actually spells.
_VERSION_SEPARATOR = "."

#: Every CURIE prefix each **Namespace** is *found* under in a published file, lower-cased
#: — a superset of :data:`NAMESPACE_PREFIX`, since a spelling nobody writes by hand is
#: still one a file uses. ``MGI:MGI:87904`` is real and is why stripping repeats.
_NAMESPACE_ALIASES: Mapping[str, frozenset[str]] = MappingProxyType(
    {
        ENSEMBL: frozenset({"ensembl"}),
        ENTREZ: frozenset({"ncbi_gene", "ncbigene", "entrezgene", "entrez", "geneid"}),
        UNIPROT: frozenset({"uniprotkb", "uniprot"}),
        HGNC: frozenset({"hgnc"}),
        MGI: frozenset({"mgi"}),
        WORMBASE: frozenset({"wb", "wormbase"}),
    }
)

#: What a CURIE prefix looks like, so that stripping one never eats an id that merely
#: holds a colon. Anchored, and the remainder must be non-empty.
_CURIE = re.compile(r"^([A-Za-z][A-Za-z0-9_.-]*):(?=.)")


def gene_id_stem(gene_id: str) -> str:
    """Return ``gene_id`` with its version dropped — its **Gene id stem**.

    Everything before the first ``.``, and the whole id when it carries none, which is what
    makes an unversioned id its own stem. The **Annotation** half of this package reduces a
    GTF's gene ids by the same rule, which is what lets a stem answered here be handed
    straight to :meth:`~genome.annotation.registry.AnnotationRegistry.resolve_gene_ids`.

    Idempotent: the stem of a stem is that stem.

    Parameters
    ----------
    gene_id : str
        A gene id, versioned or not.

    Returns
    -------
    str
        Its stem.

    Examples
    --------
    >>> gene_id_stem("ENSG00000141510.18")
    'ENSG00000141510'
    >>> gene_id_stem("WBGene00000001")
    'WBGene00000001'
    >>> gene_id_stem(gene_id_stem("ENSMUSG00000059552.16"))
    'ENSMUSG00000059552'
    """
    stem, separator, _version = gene_id.partition(_VERSION_SEPARATOR)
    return stem if separator else gene_id


def normalise_id(xref_id: str, namespace: str) -> str:
    r"""Return ``xref_id`` in ``namespace``'s canonical spelling, version dropped.

    The one door every identifier comes through, a publisher's and a caller's alike. Four
    things happen, in this order: surrounding whitespace goes; every CURIE prefix the
    namespace is published under is stripped, repeatedly, case-insensitively; the version
    suffix goes (:func:`gene_id_stem`) **and whitespace goes again**; and the namespace's
    own canonical prefix goes back on. So ``HGNC:1100``, ``hgnc:1100`` and ``1100`` are one
    identifier, and ``ENSEMBL:ENSG00000141510.18`` and ``ENSG00000141510`` are another.

    Idempotent, which is the property that matters: the same id read from a file and typed
    by a caller must land on the same string, or the join returns nothing and says nothing.
    The second strip is what makes that true rather than nearly true — a version separator
    hides trailing whitespace behind it, so ``"7157\r."`` stems to ``"7157\r"`` on the
    first pass and only reaches ``"7157"`` on the second, and two spellings of one id that
    settle on different strings after a different number of passes join to nothing.

    Parameters
    ----------
    xref_id : str
        The identifier, in whatever spelling it arrived in.
    namespace : str
        The **Namespace** it belongs to, one of :data:`NAMESPACES`. An identifier does not
        say which it is, so it is named rather than sniffed. A namespace this module does
        not know is not an error here: the id is stripped of nothing and its version is
        still dropped, because the caller of this function has already checked which
        namespaces its set carries and this is not the place to check it twice.

    Returns
    -------
    str
        The canonical spelling. Empty for an id that is empty or whitespace.

    Examples
    --------
    >>> normalise_id("HGNC:1100", "hgnc")
    'HGNC:1100'
    >>> normalise_id("1100", "hgnc")
    'HGNC:1100'
    >>> normalise_id("MGI:MGI:88276", "mgi")
    'MGI:88276'
    >>> normalise_id("UniProtKB:P38398", "uniprot")
    'P38398'
    >>> normalise_id("  ENSEMBL:ENSG00000141510.18  ", "ensembl")
    'ENSG00000141510'
    >>> normalise_id("7157\r.", "entrez")
    '7157'
    """
    body = xref_id.strip()
    aliases = _NAMESPACE_ALIASES.get(namespace, frozenset())
    while True:
        found = _CURIE.match(body)
        if found is None or found.group(1).lower() not in aliases:
            break
        body = body[found.end() :]
    # Stripped again after the version goes: the separator hides trailing whitespace behind
    # it, and an id that only settles on its second pass is not one spelling but two.
    body = gene_id_stem(body).strip()
    return f"{NAMESPACE_PREFIX.get(namespace, '')}{body}" if body else ""
