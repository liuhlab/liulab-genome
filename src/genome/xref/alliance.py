r"""The Alliance of Genome Resources as an **Xref source** — one file, three species.

The first source, and the **Default xref source** for all three of the lab's species.
``GENECROSSREFERENCE_COMBINED`` is one gene-level file covering ten organisms across every
**Namespace** this package carries — Ensembl, Entrez, UniProt and each species'
authority — pinned by a machine-readable version endpoint and a published md5.

**That md5 is of the unpacked TSV and not of the ``.tsv.gz`` that is served**, which is
the convention this package already holds itself to (ADR-0006) and is also a trap: hashing
the bytes as they arrive mismatches every time, and the served file's own md5 is the S3
ETag instead. :func:`genome.xref.xref._unpacked_lines` decompresses first and hashes what
comes out. Do not "fix" that back.

A pure reader: it is handed lines and hands back triples, opens nothing and downloads
nothing. :mod:`genome.xref.xref` owns fetching, unpacking, hashing and writing the slice;
this module owns *what the lines mean*.

**The hub is derived, not published.** Alliance keys every row by its own gene curie —
``HGNC:1100``, ``MGI:98919``, ``WB:WBGene00000001`` — and the **Gene id stem** is the
``ENSEMBL:`` cross-reference on that gene. So a gene the Alliance lists with no Ensembl
cross-reference has no hub and contributes nothing: 3,904 of 44,569 human genes and 10,668
of 88,150 mouse genes in release 9.0.0, and **none at all** of the 46,926 worm genes.

**The three species' hops have three different shapes**, which is the argument for this
being an object a caller opens rather than a step performed invisibly. For worm the hop is
the identity — all 46,926 ``WB:WBGene…`` genes carry ``ENSEMBL:WBGene…``, the same string,
so a **Gene id stem** and a WormBase gene id are one thing. For mouse it is a real join
onto ``ENSMUSG…``. For human it reaches Ensembl only through HGNC, and 2,535 human genes
carry more than one ``ENSEMBL:`` cross-reference, so 6.2% of HGNC ids name two stems or
more and nothing here picks one.

Three traps this reader exists to encode:

- **The authority id is the** ``GeneID`` **column and never the authority-prefixed
  cross-reference.** Worm's ``WB:``-prefixed cross-reference rows carry *symbols* —
  ``WB:WBGene00000001 → WB:aap-1`` — so reading the authority off the cross-reference
  column would key 47,156 worm rows by gene symbol while calling them WormBase gene ids.
- **The file ships duplicate rows, and a whole-row ``uniq`` removes none of them.** The
  same pair recurs once per web page the Alliance links it from, and the two columns that
  differ are the ones this reader drops: ``HGNC:1100 → NCBI_Gene:672`` appears under
  ``generic_cross_reference`` and again under ``gene/other_expression``, so the rows are
  distinct and the *pairs* are not. Counted on release 9.0.0: **2,659,704 rows reduce to
  1,811,267 distinct** ``(GeneID, GlobalCrossReferenceID, TaxonID)``, 31.9% redundant.
  Deduplication is therefore on the key and never on the row, and it happens here rather
  than in a caller, so the answer counts genes rather than links.
- **A cross-reference's prefix says nothing about the species.** ``RGD:`` is the single
  most frequent cross-reference prefix on human rows. Species is column five and only
  column five.

Examples
--------
>>> from genome.xref.alliance import ALLIANCE, read_alliance
>>> lines = [
...     "# a comment the file opens with",
...     "GeneID\tGlobalCrossReferenceID\tCrossReferenceCompleteURL"
...     "\tResourceDescriptorPage\tTaxonID",
...     "WB:WBGene00000001\tENSEMBL:WBGene00000001\thttps://e\tgeneric\tNCBITaxon:6239",
...     "WB:WBGene00000001\tUniProtKB:G5EDP9\thttps://u\tgeneric\tNCBITaxon:6239",
...     "WB:WBGene00000001\tWB:aap-1\thttps://s\tgene/spell\tNCBITaxon:6239",
... ]
>>> read_alliance(lines, ncbi_taxid=6239, origin="example")
(('ensembl', 'WBGene00000001', 'WBGene00000001'), ('uniprot', 'G5EDP9', 'WBGene00000001'), ('wormbase', 'WBGene00000001', 'WBGene00000001'))
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from types import MappingProxyType

from genome.xref.ids import ENSEMBL, ENTREZ, HGNC, MGI, UNIPROT, WORMBASE, normalise_id

#: What this **Xref source** is called, in the curated table, on every answer and in the
#: directory its sets are filed under.
ALLIANCE = "alliance"

#: The publisher's own column names, in file order. Read by name is not enough here — the
#: file carries no header until fourteen comment lines in — so the header line is found by
#: its first field and then held to this exactly, and a re-spelled column raises rather
#: than shifting every value one place left.
ALLIANCE_COLUMNS: tuple[str, ...] = (
    "GeneID",
    "GlobalCrossReferenceID",
    "CrossReferenceCompleteURL",
    "ResourceDescriptorPage",
    "TaxonID",
)

#: How the file spells a taxon in its fifth column, e.g. ``NCBITaxon:9606``.
TAXON_PREFIX = "NCBITaxon:"

#: The CURIE prefix of the ``GeneID`` column, per species, and the **Namespace** it is the
#: species authority for. This is what makes the authority a first-class namespace: the
#: Alliance's own key *is* the authority's id, so it needs no cross-reference row to be
#: reachable — which is as well, since the rows carrying that prefix are symbols for worm.
AUTHORITY_NAMESPACES: Mapping[str, str] = MappingProxyType(
    {"HGNC": HGNC, "MGI": MGI, "WB": WORMBASE}
)

#: The cross-reference prefixes this reader carries, and the **Namespace** each is. Every
#: other prefix in the file is skipped: PANTHER, RGD, ExpressionAtlas_gene, RNAcentral and
#: RefSeq are not gene-level identifier systems this package answers in, and a level
#: discriminator in the data would make the illegal state representable.
XREF_NAMESPACES: Mapping[str, str] = MappingProxyType(
    {"ENSEMBL": ENSEMBL, "NCBI_Gene": ENTREZ, "UniProtKB": UNIPROT}
)

#: Which of the two comment characters opens a header line the file expects to be skipped.
_COMMENT = "#"


class AllianceFileError(ValueError):
    """The Alliance file is not the file this reader reads.

    A bad *file*, not a bad call: a missing or re-spelled header, a row with the wrong
    number of fields, or a ``GeneID`` under a prefix no species authority claims. Each
    means the publisher changed the file's shape, and the message names the file and what
    was wrong with it so the fix is a reader change rather than a guess.

    Examples
    --------
    >>> try:
    ...     read_alliance(["nothing like a header"], ncbi_taxid=9606, origin="x")
    ... except AllianceFileError as error:
    ...     print("GeneID" in str(error))
    True
    """


def read_alliance(
    lines: Iterable[str], *, ncbi_taxid: int, origin: str
) -> tuple[tuple[str, str, str], ...]:
    r"""Read the Alliance cross-reference file into one species' ``(namespace, id, stem)``.

    One streaming pass over the whole file, keeping only the rows whose ``TaxonID`` is the
    species asked for. Every identifier goes through
    :func:`~genome.xref.ids.normalise_id` on the way in, so what comes out is spelled the
    one way this package spells it, versionless and with each **Namespace**'s own
    conventional prefix.

    The result is **deduplicated and sorted**, so the file's repeated rows collapse and two
    machines reading one release write byte-identical bytes.

    Parameters
    ----------
    lines : iterable of str
        The file's lines, comments and header included, newline or not.
    ncbi_taxid : int
        The species to keep, matched against ``TaxonID`` as ``NCBITaxon:<taxid>``. Read
        from the curated metadata row and never inferred from a cross-reference's prefix,
        which says nothing about species.
    origin : str
        Where the lines came from; named in every message, since a shape change here is
        fixed by editing this reader.

    Returns
    -------
    tuple of tuple of str
        ``(namespace, xref_id, gene_id_stem)`` triples, sorted and unique. Empty for a
        taxon the file does not carry, which the caller reads as *this is the wrong file*.

    Raises
    ------
    AllianceFileError
        If the header is missing or re-spelled, a data row has the wrong number of fields,
        or a ``GeneID`` carries a prefix no species authority claims.

    Examples
    --------
    >>> header = "\t".join(ALLIANCE_COLUMNS)
    >>> rows = [
    ...     "# comment",
    ...     header,
    ...     "HGNC:11998\tENSEMBL:ENSG00000141510\tu\tgeneric\tNCBITaxon:9606",
    ...     "HGNC:11998\tNCBI_Gene:7157\tu\tgeneric\tNCBITaxon:9606",
    ...     "HGNC:11998\tNCBI_Gene:7157\tu\tgene/other_expression\tNCBITaxon:9606",
    ...     "MGI:98834\tENSEMBL:ENSMUSG00000059552\tu\tgeneric\tNCBITaxon:10090",
    ... ]
    >>> read_alliance(rows, ncbi_taxid=9606, origin="example")
    (('ensembl', 'ENSG00000141510', 'ENSG00000141510'), ('entrez', '7157', 'ENSG00000141510'), ('hgnc', 'HGNC:11998', 'ENSG00000141510'))
    """
    taxon = f"{TAXON_PREFIX}{ncbi_taxid}"
    seen_header = False
    # One entry per gene of this species: the stems its ENSEMBL rows name, and every other
    # carried cross-reference on it. Held rather than streamed straight out because a
    # gene's hub may be written after the cross-reference that hangs off it.
    stems: dict[str, set[str]] = {}
    spokes: dict[str, set[tuple[str, str]]] = {}
    for line in lines:
        row = line.rstrip("\n")
        if not row or row.startswith(_COMMENT):
            continue
        fields = row.split("\t")
        if not seen_header:
            _check_header(fields, origin=origin)
            seen_header = True
            continue
        if len(fields) != len(ALLIANCE_COLUMNS):
            raise AllianceFileError(
                f"{origin} has a row of {len(fields)} fields where every row of this file "
                f"has {len(ALLIANCE_COLUMNS)}: {row[:120]!r}. The publisher changed the "
                f"file's shape, so this reader must change with it."
            )
        gene, xref, _url, _page, row_taxon = fields
        if row_taxon != taxon:
            continue
        namespace = _authority_of(gene, origin=origin)
        stems.setdefault(gene, set())
        spoke = spokes.setdefault(gene, {(namespace, normalise_id(gene, namespace))})
        prefix, _, _body = xref.partition(":")
        carried = XREF_NAMESPACES.get(prefix)
        if carried is None:
            # Every prefix this package does not answer in, and — for worm — the
            # authority-prefixed rows, which are symbols and not gene ids.
            continue
        if carried == ENSEMBL:
            stems[gene].add(normalise_id(xref, ENSEMBL))
        else:
            spoke.add((carried, normalise_id(xref, carried)))
    if not seen_header:
        raise AllianceFileError(
            f"{origin} carries no header line: this file opens with comment lines and then "
            f"one naming {', '.join(ALLIANCE_COLUMNS)}. Nothing was read from it."
        )
    triples: set[tuple[str, str, str]] = set()
    for gene, hubs in stems.items():
        for stem in hubs:
            # The hub is its own spoke: asking which stem an Ensembl id names is answered
            # by that id, never by the other Ensembl ids the same gene happens to carry —
            # that would be a foreign-to-foreign hop composed through the hub (ADR-0017).
            triples.add((ENSEMBL, stem, stem))
            for namespace, identifier in spokes[gene]:
                triples.add((namespace, identifier, stem))
    return tuple(sorted(triples))


def _check_header(fields: list[str], *, origin: str) -> None:
    """Hold the first non-comment line to the publisher's own column names."""
    if tuple(fields) != ALLIANCE_COLUMNS:
        raise AllianceFileError(
            f"{origin} names its columns {fields} where this reader reads "
            f"{list(ALLIANCE_COLUMNS)}. Columns are read by position once the header has "
            f"proved itself, so a re-spelled one would shift every value along silently."
        )


def _authority_of(gene: str, *, origin: str) -> str:
    """Return the **Namespace** the Alliance's own gene curie belongs to."""
    prefix, _, _body = gene.partition(":")
    namespace = AUTHORITY_NAMESPACES.get(prefix)
    if namespace is None:
        raise AllianceFileError(
            f"{origin} keys a gene by {gene!r}, and the species authorities this reader "
            f"knows are {', '.join(sorted(AUTHORITY_NAMESPACES))}. Every gene of a species "
            f"this package prepares is keyed by one of those; add the new prefix here "
            f"rather than letting an unnamed authority through unlabelled."
        )
    return namespace
