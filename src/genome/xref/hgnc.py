r"""HGNC's quarterly archive as an **Xref source** — human, and the only typed symbols.

The third source, and the one that answers the question the rest cannot: **which gene did
this table mean when it spelled the symbol that way?** HGNC publishes, beside every
approved symbol, the spellings it has retired (``prev_symbol``) and the spellings the gene
also goes by (``alias_symbol``), each typed — so a match can say which kind it was rather
than the caller guessing. No other publisher surveyed for this package types them
(ADR-0018), which is why human has all three kinds of **Symbol match** and mouse and worm
have one.

**The failure this exists to prevent is measured.** Of EpiFactors v2.0's 801 human rows,
all 801 carry an HGNC id and all 801 map — and **31 still spell the gene by a symbol HGNC
has since retired**: ``ARNTL`` for ``BMAL1``, ``C11orf30`` for ``EMSY``, ``ACINUS``
for ``ACIN1``. A symbol join that knows only approved spellings mis-keys or drops exactly
those 31, silently, which is what everyone in the lab has been doing by hand.

**The reader parses by header name and never by position**, because the schema has drifted:
the header is **52** columns wide in ``hgnc_complete_set_2020-07-01.txt`` and **54** in
``…_2026-07-07.txt``. A reader indexing by position reads the wrong column on the older
snapshots, so this one reads the header, finds :data:`HGNC_COLUMNS` in it by name, and
raises if one is missing. Reordering a column therefore cannot change the answer.

**The pin names a file from the archive listing and never one built from a date.** The
archive's dates are irregular — ``2024-07-02``, ``2025-01-06``, ``2025-10-07``, and both
``2026-07-03`` and ``2026-07-07`` — so a URL assembled from *the first of the quarter*
404s about half the time. List the bucket and copy a name out of it:
``https://storage.googleapis.com/storage/v1/b/public-download-files/o?prefix=hgnc/archive/archive/quarterly/tsv``.

**The path has a doubled segment and the remembered one is dead.**
``ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt`` is a live 404;
HGNC serves from a Google Cloud Storage bucket whose archive path repeats itself —
``…/hgnc/archive/archive/quarterly/tsv/…``, both segments, not a typo here.

**A multi-valued cell is quoted and pipe-separated, and a single value is bare.**
``alias_symbol`` reads ``"MOP3|JAP3|PASD3|bHLHe5|ARNTL1"`` with the quotes in the file, and
``AC3`` without them. Stripping the quotes before splitting is not tidying: leaving them
on would key the namespace by ``"MOP3`` and ``ARNTL1"``, neither of which anybody types.

**The published checksum covers the served bytes, which here are also the unpacked ones.**
The file is plain text, not gzipped, so the bucket's own md5 — from the same listing the
pin is read out of — is the digest of exactly what a curated row must pin, and the two
conventions this directory already holds side by side (the Alliance's over the unpacked
bytes, Ensembl's over the served ones) coincide for once. Do not carry that coincidence to
a source whose file is compressed.

Examples
--------
>>> from genome.xref.hgnc import HGNC_ARCHIVE, read_hgnc
>>> lines = [
...     "\t".join(("hgnc_id", "symbol", "alias_symbol", "prev_symbol", "ensembl_gene_id",
...                "entrez_id", "uniprot_ids")),
...     "\t".join(("HGNC:701", "BMAL1", '"MOP3|ARNTL1"', "ARNTL", "ENSG00000133794",
...                "406", "O00327")),
... ]
>>> for triple in read_hgnc(lines, ncbi_taxid=9606, origin="example"):
...     print(triple)
('alias_symbol', 'ARNTL1', 'ENSG00000133794')
('alias_symbol', 'MOP3', 'ENSG00000133794')
('ensembl', 'ENSG00000133794', 'ENSG00000133794')
('entrez', '406', 'ENSG00000133794')
('hgnc', 'HGNC:701', 'ENSG00000133794')
('previous_symbol', 'ARNTL', 'ENSG00000133794')
('symbol', 'BMAL1', 'ENSG00000133794')
('uniprot', 'O00327', 'ENSG00000133794')
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from types import MappingProxyType

from genome.xref.evidence import EvidenceNotRecordedError
from genome.xref.ids import ENSEMBL, ENTREZ, HGNC, UNIPROT, normalise_id
from genome.xref.symbols import ALIAS, APPROVED, KIND_NAMESPACES, PREVIOUS, normalise_symbol

#: What this **Xref source** is called, in the curated table, on every answer and in the
#: directory its sets are filed under. Spelled like the HGNC **Namespace**
#: (:data:`~genome.xref.ids.HGNC`) and meaning something else — this is the publisher, that
#: is the identifier system — exactly as ``ensembl`` is both a source and the hub namespace.
HGNC_ARCHIVE = "hgnc"

#: The column holding the gene's own HGNC id, which is this file's key.
HGNC_ID_COLUMN = "hgnc_id"

#: The column holding the **Gene id stem**. A row with this blank has no hub and is
#: skipped: 2,682 of the 45,019 rows of the ``2026-07-07`` archive, mostly pseudogenes and
#: non-coding RNAs Ensembl does not carry.
ENSEMBL_COLUMN = "ensembl_gene_id"

#: Every column this reader reads, by the name HGNC gives it. A file missing one of these
#: raises; a file carrying fifty others is read anyway, which is what lets one reader read
#: both the 52-column and the 54-column schema.
HGNC_COLUMNS: tuple[str, ...] = (
    HGNC_ID_COLUMN,
    "symbol",
    "prev_symbol",
    "alias_symbol",
    ENSEMBL_COLUMN,
    "entrez_id",
    "uniprot_ids",
)

#: Which **Namespace** each id-bearing column is read into.
ID_COLUMNS: Mapping[str, str] = MappingProxyType(
    {HGNC_ID_COLUMN: HGNC, ENSEMBL_COLUMN: ENSEMBL, "entrez_id": ENTREZ, "uniprot_ids": UNIPROT}
)

#: Which kind of **Symbol match** each symbol-bearing column publishes. This is the whole
#: reason HGNC is here: the kinds arrive typed and none has to be inferred.
SYMBOL_COLUMNS: Mapping[str, str] = MappingProxyType(
    {"symbol": APPROVED, "prev_symbol": PREVIOUS, "alias_symbol": ALIAS}
)

#: Every column whose cell may hold several values. Quoted and pipe-separated when it does,
#: bare when it holds one.
MULTI_VALUED: frozenset[str] = frozenset({"prev_symbol", "alias_symbol", "uniprot_ids"})

#: What HGNC wraps a multi-valued cell in.
_QUOTE = '"'

#: What separates the values inside one.
_SEPARATOR = "|"

#: This source carries every kind of **Symbol match** there is, so nothing rides back on an
#: answer explaining what is missing.
HGNC_SYMBOL_LIMIT: str | None = None


class HgncFileError(ValueError):
    r"""The HGNC archive file is not the file this reader reads.

    A bad *file*, not a bad call: no header line, or a header that names none of the
    columns this reader needs. Not raised for a column this reader ignores, and not raised
    for a re-ordered one — reading by name is what makes both harmless, and a schema that
    has already gone from 52 columns to 54 will move again.

    Examples
    --------
    >>> try:
    ...     read_hgnc(["nothing\tlike\ta\theader"], ncbi_taxid=9606, origin="x")
    ... except HgncFileError as error:
    ...     print("prev_symbol" in str(error))
    True
    """


def read_hgnc(
    lines: Iterable[str],
    *,
    ncbi_taxid: int,
    origin: str,
    evidence: tuple[str, ...] = (),
) -> tuple[tuple[str, str, str], ...]:
    r"""Read one HGNC quarterly archive file into ``(namespace, id, stem)`` triples.

    One streaming pass. The header is read by name into a column index, every later row is
    read through that index, and both sides of each pair go through
    :func:`~genome.xref.ids.normalise_id` — or, for a symbol,
    :func:`~genome.xref.symbols.normalise_symbol`, which leaves the case alone and does not
    cut the spelling at a dot. The result is **deduplicated and sorted**, so two machines
    reading one archive file write byte-identical bytes.

    Symbols are stored under three **Namespace**s rather than one, so that the kind of
    spelling survives into the set on disk and out of
    :meth:`~genome.xref.xref.XrefSet.match_symbols` — see
    :data:`~genome.xref.symbols.KIND_NAMESPACES`.

    Parameters
    ----------
    lines : iterable of str
        The file's lines, header included, newline or not.
    ncbi_taxid : int
        Accepted for the reader protocol and **not consulted**: HGNC names human genes and
        nothing else, and which file arrived is decided by the URL the curated row pins.
    origin : str
        Where the lines came from; named in every message.
    evidence : tuple of str, optional
        Must be empty. **This file grades nothing** — a row is a gene, its spellings and
        its cross-references, with no column saying how any of them was arrived at — so a
        filter cannot be honoured and is refused rather than ignored.

    Returns
    -------
    tuple of tuple of str
        ``(namespace, xref_id, gene_id_stem)`` triples, sorted and unique, over the
        ``ensembl``, ``entrez``, ``uniprot`` and ``hgnc`` **Namespace**s and the three
        symbol ones. Empty for a file whose every row lacks a **Gene id stem**, which the
        caller reads as *this is the wrong file*.

    Raises
    ------
    HgncFileError
        If the file has no header, or its header names none of :data:`HGNC_COLUMNS`.
    genome.xref.evidence.EvidenceNotRecordedError
        If ``evidence`` names anything at all.

    Examples
    --------
    >>> header = "\t".join(HGNC_COLUMNS)
    >>> rows = [
    ...     header,
    ...     "HGNC:11998\tTP53\t\t\"p53|LFS1\"\tENSG00000141510\t7157\tP04637",
    ...     "HGNC:22\tAAVS1\t\t\t\t\t",
    ... ]
    >>> read_hgnc(rows, ncbi_taxid=9606, origin="example")[:3]
    (('alias_symbol', 'LFS1', 'ENSG00000141510'), ('alias_symbol', 'p53', 'ENSG00000141510'), ('ensembl', 'ENSG00000141510', 'ENSG00000141510'))
    >>> [triple for triple in read_hgnc(rows, ncbi_taxid=9606, origin="e") if "AAVS1" in triple]
    []
    """
    if evidence:
        raise EvidenceNotRecordedError(
            f"the {HGNC_ARCHIVE} archive file records no evidence type, so the filter "
            f"{'/'.join(evidence)} cannot be applied to it: a row carries a gene's "
            f"spellings and its cross-references, and nothing that grades how any of them "
            f"was arrived at. Construct the set without an evidence filter, or name a "
            f"source whose file records one — the ensembl per-species TSV grades every row."
        )
    columns: dict[str, int] | None = None
    triples: set[tuple[str, str, str]] = set()
    for line in lines:
        row = line.rstrip("\n")
        if not row:
            continue
        fields = row.split("\t")
        if columns is None:
            columns = _header_index(fields, origin=origin)
            continue
        stem = normalise_id(_cell(fields, columns, ENSEMBL_COLUMN), ENSEMBL)
        if not stem:
            # No hub, so nothing to hang a namespace off — 2,682 rows of the pinned
            # archive, and the same silence the Alliance reader keeps for a gene it lists
            # with no ENSEMBL cross-reference.
            continue
        triples.add((ENSEMBL, stem, stem))
        for column, namespace in ID_COLUMNS.items():
            for value in _values(_cell(fields, columns, column), column):
                identifier = normalise_id(value, namespace)
                if identifier:
                    triples.add((namespace, identifier, stem))
        for column, kind in SYMBOL_COLUMNS.items():
            for value in _values(_cell(fields, columns, column), column):
                spelling = normalise_symbol(value)
                if spelling:
                    triples.add((KIND_NAMESPACES[kind], spelling, stem))
    if columns is None:
        raise HgncFileError(
            f"{origin} carries no header line: this file opens with one naming "
            f"{', '.join(HGNC_COLUMNS)} among its columns. Nothing was read from it."
        )
    return tuple(sorted(triples))


def _header_index(fields: list[str], *, origin: str) -> dict[str, int]:
    """Return where each column this reader needs sits, or say which one is missing.

    By name, never by position: HGNC's header has already gone from 52 columns to 54, so a
    reader that counted would read the wrong column on the older snapshots and say nothing
    about it.
    """
    found = {name: index for index, name in enumerate(fields)}
    missing = [name for name in HGNC_COLUMNS if name not in found]
    if missing:
        raise HgncFileError(
            f"{origin} names its columns {fields[:8]}… and this reader needs "
            f"{', '.join(missing)}, which it does not name. Columns are read by name here "
            f"rather than by position, so a re-ordered or added column is harmless and a "
            f"renamed one is not: either the file is not an HGNC archive file, or the "
            f"publisher renamed a column and this reader must change with it."
        )
    return {name: found[name] for name in HGNC_COLUMNS}


def _cell(fields: list[str], columns: dict[str, int], name: str) -> str:
    """Return one column's cell, empty for a row that stops before it."""
    index = columns[name]
    return fields[index] if index < len(fields) else ""


def _values(cell: str, column: str) -> tuple[str, ...]:
    """Return the values in one cell — several when the column allows it, else the one.

    A multi-valued cell arrives quoted and pipe-separated (``"MOP3|JAP3"``) when it holds
    more than one and bare (``AC3``) when it holds one, so the quotes come off before the
    split rather than riding into an id nobody would type.
    """
    text = cell.strip()
    if not text:
        return ()
    if column not in MULTI_VALUED:
        return (text,)
    if len(text) > 1 and text.startswith(_QUOTE) and text.endswith(_QUOTE):
        text = text[1:-1]
    return tuple(value for value in (part.strip() for part in text.split(_SEPARATOR)) if value)
