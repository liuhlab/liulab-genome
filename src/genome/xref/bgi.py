r"""The Alliance's per-species gene submission as an **Xref source** — mouse and worm.

The fourth source, and the one that gives mouse and worm a current approved symbol at all.
``BGI`` — Basic Gene Information — is each model-organism database's own gene submission,
served by the Alliance at a versioned path with a published md5, one file per contributing
database. Mouse's is MGI's; worm's is WormBase's, and its own header states
``"release": "WS298"`` — the final WormBase release, which is the one the lab has
registered.

**It is here because the two authorities themselves cannot be pinned or cannot be
fetched.** MGI keeps no dated archive and its files carry no date stamp at all, so it is
not an eligible **Xref source** (ADR-0018). ``downloads.wormbase.org`` answers **403** to
an automated client — the release directories, not only the blog — measured from two
networks. The Alliance holds a dated, immutable, checksummed copy of both submissions and
is in any case the ongoing publisher for *C. elegans* now that WS298 is final, so what
this source carries is those authorities' own assertions, retrievable a year from now.

**Approved spellings only, and that is a decision rather than a shortfall.** Each record
carries one ``symbol`` and one ``synonyms`` list, and the list is **undifferentiated**:
worm's ``daf-16`` record puts the sequence names ``R13H8.1`` and ``CELE_R13H8.1`` in the
same list as ``daf-17``, a name the gene genuinely went by, and nothing in the file says
which is which. Reading them would mean labelling each one ``previous`` or ``alias`` on
this package's own authority, which is precisely the claim it never makes (ADR-0017). So
:data:`BGI_SYMBOL_LIMIT` rides back on every answer this source produces, saying so.

**The file is JSON, and it is read a record at a time.** One object with a ``metaData``
header and a ``data`` array; the mouse file unpacks to 72 MB and the worm file to 74 MB,
so neither is decoded whole. :func:`read_bgi` finds the array, then peels one gene object
off the front of a small rolling buffer with :meth:`json.JSONDecoder.raw_decode` — which
also means the publishers' two different pretty-printers, two-space for MGI and
three-space for WormBase, are both simply JSON and neither is parsed by eye.

**The published md5 is of the unpacked JSON**, the same convention the Alliance's
cross-reference file uses and the opposite of Ensembl's — verified by decompressing both
files and hashing: ``a834098c9505ec7fb4a0151480a90734`` for the mouse submission and
``4a45ce6beb26dd0dc8c053e5b2e1a835`` for the worm one, each matching the file-management
API's ``md5Sum`` only after ``gzip -d``.

Examples
--------
>>> from genome.xref.bgi import ALLIANCE_BGI, read_bgi
>>> import json
>>> record = {
...     "basicGeneticEntity": {
...         "primaryId": "WB:WBGene00000912",
...         "taxonId": "NCBITaxon:6239",
...         "synonyms": ["daf-17", "R13H8.1"],
...         "crossReferences": [{"id": "ENSEMBL:WBGene00000912"}, {"id": "NCBI_Gene:172981"}],
...     },
...     "symbol": "daf-16",
... }
>>> lines = json.dumps({"metaData": {}, "data": [record]}, indent=1).splitlines(keepends=True)
>>> read_bgi(lines, ncbi_taxid=6239, origin="example")
(('ensembl', 'WBGene00000912', 'WBGene00000912'), ('entrez', '172981', 'WBGene00000912'), ('symbol', 'daf-16', 'WBGene00000912'), ('wormbase', 'WBGene00000912', 'WBGene00000912'))
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Iterator
from typing import Any

from genome.xref.alliance import AUTHORITY_NAMESPACES, TAXON_PREFIX, XREF_NAMESPACES
from genome.xref.evidence import EvidenceNotRecordedError
from genome.xref.ids import ENSEMBL, normalise_id
from genome.xref.symbols import APPROVED, KIND_NAMESPACES, normalise_symbol

#: What this **Xref source** is called, in the curated table, on every answer and in the
#: directory its sets are filed under. Named for the publisher that serves it and the
#: product it is, because ``alliance`` already names that publisher's *other* file and the
#: two assert different things about different genes.
ALLIANCE_BGI = "alliance_bgi"

#: The key the array of genes sits under, and the only part of the file that is read.
DATA_KEY = "data"

#: Where a record's identifiers live, inside each gene object.
ENTITY_KEY = "basicGeneticEntity"

#: The record's own key — the species authority's id, as a curie: ``MGI:98919``,
#: ``WB:WBGene00000001``. Read from here and never from a cross-reference row, for the
#: reason the Alliance reader beside this one records: worm's authority-prefixed
#: cross-references carry *symbols*.
PRIMARY_ID_KEY = "primaryId"

#: The record's current approved symbol, and the only spelling this source publishes with a
#: kind attached to it.
SYMBOL_KEY = "symbol"

#: What rides back on every answer this source produces, saying which kinds of **Symbol
#: match** it cannot make and why. The explanation is part of the behaviour: a caller who
#: asks a mouse set about a spelling MGI retired gets nothing back, and *why* nothing came
#: back is the difference between a gene that is absent and a gene that is spelled another
#: way.
BGI_SYMBOL_LIMIT = (
    "this source publishes one current approved symbol per gene beside an undifferentiated "
    "synonyms list, and a spelling in that list does not say whether the authority retired "
    "it or merely records it — WormBase files the sequence name R13H8.1 and the genuine "
    "former name daf-17 in the same list — so matching on it would put a kind on a claim no "
    "publisher made. The typed previous and alias spellings belong to the species "
    "authorities themselves: MGI for mouse, which keeps no dated archive and is therefore "
    "not an eligible xref source (ADR-0018), and WormBase for worm, whose download host "
    "answers 403 to an automated client. Human has all three kinds, from the hgnc source."
)


class BgiFileError(ValueError):
    """The Alliance gene submission is not the file this reader reads.

    A bad *file*, not a bad call: no ``data`` array, JSON that does not parse, or a gene
    keyed by a prefix no species authority claims. The message names the file and what was
    wrong with it, so the fix is a reader change rather than a guess.

    Examples
    --------
    >>> try:
    ...     read_bgi(["{}"], ncbi_taxid=10090, origin="x")
    ... except BgiFileError as error:
    ...     print("data" in str(error))
    True
    """


def read_bgi(
    lines: Iterable[str],
    *,
    ncbi_taxid: int,
    origin: str,
    evidence: tuple[str, ...] = (),
) -> tuple[tuple[str, str, str], ...]:
    r"""Read one Alliance gene submission into ``(namespace, id, stem)`` triples.

    One streaming pass over the ``data`` array, keeping the records whose ``taxonId`` is
    the species asked for. Identifiers go through :func:`~genome.xref.ids.normalise_id` and
    the symbol through :func:`~genome.xref.symbols.normalise_symbol`, so what comes out is
    spelled the one way this package spells it. The result is **deduplicated and sorted**,
    so two machines reading one submission write byte-identical bytes.

    Parameters
    ----------
    lines : iterable of str
        The file's lines, in order. Newlines may be kept or stripped; the JSON is
        re-assembled either way.
    ncbi_taxid : int
        The species to keep, matched against each record's ``taxonId`` as
        ``NCBITaxon:<taxid>``. A submission is already one species, so this is a check
        rather than a filter — but it is applied rather than assumed, since which file
        arrived is decided by the URL the curated row pins.
    origin : str
        Where the lines came from; named in every message.
    evidence : tuple of str, optional
        Must be empty. **This file grades nothing**, so a filter cannot be honoured and is
        refused rather than ignored.

    Returns
    -------
    tuple of tuple of str
        ``(namespace, xref_id, gene_id_stem)`` triples, sorted and unique, over the
        ``ensembl``, ``entrez``, ``uniprot`` and species-authority **Namespace**s and the
        approved-symbol one. Empty for a taxon the file does not carry, which the caller
        reads as *this is the wrong file*.

    Raises
    ------
    BgiFileError
        If the file carries no ``data`` array, does not parse as JSON, or keys a gene by a
        prefix no species authority claims.
    genome.xref.evidence.EvidenceNotRecordedError
        If ``evidence`` names anything at all.

    Examples
    --------
    >>> import json
    >>> mouse = {
    ...     "basicGeneticEntity": {
    ...         "primaryId": "MGI:98834",
    ...         "taxonId": "NCBITaxon:10090",
    ...         "crossReferences": [{"id": "ENSEMBL:ENSMUSG00000059552"}],
    ...     },
    ...     "symbol": "Trp53",
    ... }
    >>> text = json.dumps({"data": [mouse]}, indent=2)
    >>> read_bgi(text.splitlines(), ncbi_taxid=10090, origin="example")
    (('ensembl', 'ENSMUSG00000059552', 'ENSMUSG00000059552'), ('mgi', 'MGI:98834', 'ENSMUSG00000059552'), ('symbol', 'Trp53', 'ENSMUSG00000059552'))
    >>> read_bgi(text.splitlines(), ncbi_taxid=6239, origin="example")
    ()
    """
    if evidence:
        raise EvidenceNotRecordedError(
            f"the {ALLIANCE_BGI} gene submission records no evidence type, so the filter "
            f"{'/'.join(evidence)} cannot be applied to it: a record carries a gene's ids, "
            f"its symbol and its synonyms, and nothing that grades how any of them was "
            f"arrived at. Construct the set without an evidence filter, or name a source "
            f"whose file records one — the ensembl per-species TSV grades every row."
        )
    taxon = f"{TAXON_PREFIX}{ncbi_taxid}"
    triples: set[tuple[str, str, str]] = set()
    for record in _genes(lines, origin=origin):
        entity = record.get(ENTITY_KEY)
        if not isinstance(entity, dict) or entity.get("taxonId") != taxon:
            continue
        gene = str(entity.get(PRIMARY_ID_KEY, ""))
        namespace = _authority_of(gene, origin=origin)
        stems = {
            normalise_id(str(reference.get("id", "")), ENSEMBL)
            for reference in _references(entity)
            if str(reference.get("id", "")).partition(":")[0] == "ENSEMBL"
        }
        if not stems:
            # No hub, so nothing to hang a namespace off — the same silence the Alliance
            # cross-reference reader keeps for a gene it lists with no ENSEMBL row.
            continue
        spokes = {(namespace, normalise_id(gene, namespace))}
        for reference in _references(entity):
            identifier = str(reference.get("id", ""))
            carried = XREF_NAMESPACES.get(identifier.partition(":")[0])
            if carried is not None and carried != ENSEMBL:
                spokes.add((carried, normalise_id(identifier, carried)))
        symbol = normalise_symbol(str(record.get(SYMBOL_KEY, "")))
        if symbol:
            spokes.add((KIND_NAMESPACES[APPROVED], symbol))
        for stem in stems:
            # The hub is its own spoke, and never the other stems the same gene carries —
            # that would be a foreign-to-foreign hop composed through it (ADR-0017).
            triples.add((ENSEMBL, stem, stem))
            for carried, identifier in spokes:
                triples.add((carried, identifier, stem))
    return tuple(sorted(triples))


def _references(entity: dict[str, Any]) -> list[dict[str, Any]]:
    """Return the record's cross-reference objects, or none when it carries the key empty."""
    references = entity.get("crossReferences")
    return [item for item in references if isinstance(item, dict)] if references else []


def _authority_of(gene: str, *, origin: str) -> str:
    """Return the **Namespace** the submission's own gene curie belongs to."""
    namespace = AUTHORITY_NAMESPACES.get(gene.partition(":")[0])
    if namespace is None:
        raise BgiFileError(
            f"{origin} keys a gene by {gene!r}, and the species authorities this reader "
            f"knows are {', '.join(sorted(AUTHORITY_NAMESPACES))}. Every gene of a species "
            f"this package prepares is keyed by one of those; add the new prefix rather "
            f"than letting an unnamed authority through unlabelled."
        )
    return namespace


def _genes(lines: Iterable[str], *, origin: str) -> Iterator[dict[str, Any]]:
    """Yield one gene object at a time off the front of a rolling buffer.

    The array is found by its key and then peeled one object at a time, so the 72 MB the
    mouse submission unpacks to is never held: the buffer holds one pretty-printed record
    and whatever of the next has arrived.
    """
    decoder = json.JSONDecoder()
    buffer = ""
    opened = False
    for line in lines:
        buffer += line
        if not opened:
            key = buffer.find(f'"{DATA_KEY}"')
            bracket = buffer.find("[", key) if key >= 0 else -1
            if bracket < 0:
                continue
            buffer, opened = buffer[bracket + 1 :], True
        while True:
            buffer = buffer.lstrip(" \t\r\n,")
            if not buffer.startswith("{"):
                break
            try:
                record, end = decoder.raw_decode(buffer)
            except ValueError:
                # The record has not arrived whole yet; the next line finishes it.
                break
            buffer = buffer[end:]
            yield record
    if not opened:
        raise BgiFileError(
            f"{origin} carries no {DATA_KEY!r} array: an Alliance gene submission is one "
            f"JSON object with a metaData header and a {DATA_KEY!r} array of genes beside "
            f"it. Nothing was read from it."
        )
