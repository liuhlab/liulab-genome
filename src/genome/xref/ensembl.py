r"""Ensembl's per-species TSV dumps as an **Xref source** — one file per species.

The second source, and **not the equal of the first**. It is here because it pins cleanly:
``Homo_sapiens.GRCh38.<release>.entrez.tsv.gz`` answers 200 for every release probed from
**88** through **116**, and each directory ships a ``CHECKSUMS`` file beside it (ADR-0018).
It is *not* a **Default xref source** for anything, because of what follows.

**The fan-out, which is the reason to choose this source deliberately.** Measured on human
release 116 against NCBI's own ``gene2ensembl``: the two agree on only **57.6%** of the
gene-level (GeneID, ENSG) pairs they assert between them — Ensembl asserts 36,824, NCBI
38,577, and 27,554 are common. The cause is method rather than release skew. NCBI's mapping
is a sequence match at a published overlap threshold and comes out all but one-to-one, no
GeneID naming more than two stems. Ensembl's fans out by two orders of magnitude at the
tail: **72 stems for one GeneID** (``79166``) and **208 GeneIDs for one stem**
(``ENSG00000278233``), both exact. So an answer from here is wider than an answer from the
Alliance for the same id, and the width is the publisher's assertion rather than an error —
which is precisely why the two are two answers and never one merged one (ADR-0017).

**The empty-filter trap.** Ensembl grades each cross-reference in an ``info_type`` column,
so the obvious quality filter is ``DIRECT``. Every human ``EntrezGene`` row in release 116
carries :data:`DEPENDENT` and **not one** carries :data:`DIRECT` — 552,633 rows, zero
direct — and mouse is the same at 358,853. Filtering to ``DIRECT`` therefore yields an
empty set rather than a smaller one, so it raises
:class:`~genome.xref.evidence.EmptyEvidenceFilterError` naming what the release does carry.
See :mod:`genome.xref.evidence`.

**The published checksum covers the served bytes, and the Alliance's covers the unpacked
ones.** Do not carry either convention across. Ensembl's ``CHECKSUMS`` is BSD ``sum``
format over the ``.tsv.gz`` exactly as downloaded — ``28782 5952`` for release 116's human
dump, verified against the file — while the Alliance publishes an md5 of the TSV *inside*
its gzip. A BSD ``sum`` is a 16-bit checksum and no integrity check at all for a 6 MB file,
so the curated row pins an md5 of the unpacked bytes computed here on the pinned release
instead, and the attribution beside the table records Ensembl's own value and what it
covers.

Three properties of the file this reader exists to encode:

- **Only the** ``EntrezGene`` **rows are cross-references.** The dump also carries
  :data:`TRANSCRIPT_NAME_DB_NAME` rows — 644,668 of them for human — whose ``xref`` column
  holds a *transcript name* such as ``KU-MEL-3-201``, not a GeneID. The file's own
  ``README_entrez.tsv`` warns that it "contains all Ensembl external database names which
  started with entrez so duplication of hits is possible". Reading them as Entrez ids would
  key the namespace by transcript labels; they are a different assertion and are dropped.
- **The file is transcript- and protein-grained.** 552,633 human rows collapse to 36,824
  gene-level pairs. Deduplication is on the pair and happens here, so an answer counts
  genes rather than transcripts.
- **The file is already one species**, so ``ncbi_taxid`` is accepted for the reader
  protocol's sake and never consulted: which species arrived is decided by the URL the
  curated row pins and held to it by that row's checksum, not sniffed out of an id prefix.

One name to keep straight: this source is spelled ``ensembl`` and so is the hub
**Namespace** (:data:`~genome.xref.ids.ENSEMBL`). They are the same publisher named twice —
once as *who asserted this* and once as *which identifier system* — and they are not
interchangeable. :data:`ENSEMBL_TSV` is the source.

Examples
--------
>>> from genome.xref.ensembl import ENSEMBL_TSV, read_ensembl
>>> lines = [
...     "\t".join(
...         (
...             "gene_stable_id", "transcript_stable_id", "protein_stable_id", "xref",
...             "db_name", "info_type", "source_identity", "xref_identity", "linkage_type",
...         )
...     ),
...     "ENSG00000141510\tENST00000269305\tENSP00000269305\t7157"
...     "\tEntrezGene\tDEPENDENT\t-\t-\t-",
...     "ENSG00000141510\tENST00000445888\tENSP00000391127\t7157"
...     "\tEntrezGene\tDEPENDENT\t-\t-\t-",
...     "ENSG00000141510\tENST00000269305\t-\tTP53-201"
...     "\tEntrezGene_trans_name\tMISC\t-\t-\t-",
... ]
>>> read_ensembl(lines, ncbi_taxid=9606, origin="example")
(('ensembl', 'ENSG00000141510', 'ENSG00000141510'), ('entrez', '7157', 'ENSG00000141510'))
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable

from genome.xref.evidence import check_evidence_filter
from genome.xref.ids import ENSEMBL, ENTREZ, normalise_id

#: What this **Xref source** is called, in the curated table, on every answer and in the
#: directory its sets are filed under. Spelled like the hub **Namespace** and meaning
#: something else: this is the publisher, that is the identifier system.
ENSEMBL_TSV = "ensembl"

#: The publisher's own column names, in file order. The first line of the file is the
#: header — there are no comment lines — and it is held to this exactly, so a re-spelled
#: column raises rather than shifting every value one place along.
ENSEMBL_COLUMNS: tuple[str, ...] = (
    "gene_stable_id",
    "transcript_stable_id",
    "protein_stable_id",
    "xref",
    "db_name",
    "info_type",
    "source_identity",
    "xref_identity",
    "linkage_type",
)

#: The ``db_name`` whose rows are Entrez cross-references, and the only rows read.
ENTREZ_DB_NAME = "EntrezGene"

#: The ``db_name`` whose rows carry a *transcript name* rather than a GeneID, and are
#: dropped. A different assertion, and the trap this reader exists to sidestep.
TRANSCRIPT_NAME_DB_NAME = "EntrezGene_trans_name"

#: Ensembl's two ``info_type`` gradings for a cross-reference. ``DIRECT`` means the
#: publisher matched the entities themselves; ``DEPENDENT`` means the link was inherited
#: through another one. **No human or mouse** ``EntrezGene`` **row of release 116 is
#: ``DIRECT``**, which is why these are named here rather than left as strings in a filter.
DIRECT, DEPENDENT = "DIRECT", "DEPENDENT"

#: How the file spells "no value" in a column, e.g. an unnamed protein.
_ABSENT = "-"


class EnsemblTsvFileError(ValueError):
    """The Ensembl TSV is not the file this reader reads.

    A bad *file*, not a bad call: a missing or re-spelled header, or a row with the wrong
    number of fields. Each means the publisher changed the file's shape, so the message
    names the file and what was wrong with it and the fix is a reader change rather than a
    guess.

    Examples
    --------
    >>> try:
    ...     read_ensembl(["nothing like a header"], ncbi_taxid=9606, origin="x")
    ... except EnsemblTsvFileError as error:
    ...     print("gene_stable_id" in str(error))
    True
    """


def read_ensembl(
    lines: Iterable[str],
    *,
    ncbi_taxid: int,
    origin: str,
    evidence: tuple[str, ...] = (),
) -> tuple[tuple[str, str, str], ...]:
    r"""Read one Ensembl per-species TSV into ``(namespace, id, stem)`` triples.

    One streaming pass, keeping the :data:`ENTREZ_DB_NAME` rows and dropping every other
    ``db_name``. Both sides of each pair go through
    :func:`~genome.xref.ids.normalise_id` on the way in, so a versioned gene id and a bare
    one are one identifier. The result is **deduplicated and sorted**, which is what
    collapses the file's transcript- and protein-grained rows to gene-level pairs and what
    makes two machines reading one release write byte-identical bytes.

    Parameters
    ----------
    lines : iterable of str
        The file's lines, header included, newline or not.
    ncbi_taxid : int
        Accepted for the reader protocol and **not consulted**: this file is already one
        species, and which one is decided by the URL the curated row pins.
    origin : str
        Where the lines came from; named in every message.
    evidence : tuple of str, optional
        The ``info_type`` gradings to keep, as
        :func:`~genome.xref.evidence.normalise_evidence` spells them. Empty keeps every
        row, which is the only filter that answers anything for the releases shipped here.

    Returns
    -------
    tuple of tuple of str
        ``(namespace, xref_id, gene_id_stem)`` triples, sorted and unique, over the
        ``ensembl`` and ``entrez`` **Namespace**s. Empty for a file with no such rows,
        which the caller reads as *this is the wrong file*.

    Raises
    ------
    EnsemblTsvFileError
        If the header is missing or re-spelled, or a row has the wrong number of fields.
    genome.xref.evidence.EmptyEvidenceFilterError
        If ``evidence`` keeps none of the rows — which is what ``DIRECT`` does to every
        release shipped here. The message names what the release actually carries.

    Examples
    --------
    >>> header = "\t".join(ENSEMBL_COLUMNS)
    >>> rows = [
    ...     header,
    ...     "ENSG00000170858.7\tENST1\t-\t79166\tEntrezGene\tDEPENDENT\t-\t-\t-",
    ...     "ENSG00000293273\tENST2\t-\t79166\tEntrezGene\tDEPENDENT\t-\t-\t-",
    ... ]
    >>> read_ensembl(rows, ncbi_taxid=9606, origin="example")
    (('ensembl', 'ENSG00000170858', 'ENSG00000170858'), ('ensembl', 'ENSG00000293273', 'ENSG00000293273'), ('entrez', '79166', 'ENSG00000170858'), ('entrez', '79166', 'ENSG00000293273'))
    """
    seen_header = False
    graded: Counter[str] = Counter()
    triples: set[tuple[str, str, str]] = set()
    for line in lines:
        row = line.rstrip("\n")
        if not row:
            continue
        fields = row.split("\t")
        if not seen_header:
            _check_header(fields, origin=origin)
            seen_header = True
            continue
        if len(fields) != len(ENSEMBL_COLUMNS):
            raise EnsemblTsvFileError(
                f"{origin} has a row of {len(fields)} fields where every row of this file "
                f"has {len(ENSEMBL_COLUMNS)}: {row[:120]!r}. The publisher changed the "
                f"file's shape, so this reader must change with it."
            )
        gene, _transcript, _protein, xref, db_name, info_type = fields[:6]
        if db_name != ENTREZ_DB_NAME:
            # Every other db_name, and in particular the transcript-name rows, whose xref
            # column holds a transcript label rather than a GeneID.
            continue
        graded[info_type] += 1
        if evidence and info_type not in evidence:
            continue
        stem = normalise_id(gene, ENSEMBL)
        identifier = normalise_id(xref, ENTREZ)
        if not stem or not identifier or stem == _ABSENT or identifier == _ABSENT:
            continue
        # The hub is its own spoke: asking which stem an Ensembl id names is answered by
        # that id and never by the other stems the same GeneID happens to reach, which
        # would be a foreign-to-foreign hop composed through the hub (ADR-0017).
        triples.add((ENSEMBL, stem, stem))
        triples.add((ENTREZ, identifier, stem))
    if not seen_header:
        raise EnsemblTsvFileError(
            f"{origin} carries no header line: this file opens with one naming "
            f"{', '.join(ENSEMBL_COLUMNS)}. Nothing was read from it."
        )
    check_evidence_filter(wanted=evidence, seen=graded, origin=origin)
    return tuple(sorted(triples))


def _check_header(fields: list[str], *, origin: str) -> None:
    """Hold the first line to the publisher's own column names."""
    if tuple(fields) != ENSEMBL_COLUMNS:
        raise EnsemblTsvFileError(
            f"{origin} names its columns {fields} where this reader reads "
            f"{list(ENSEMBL_COLUMNS)}. Columns are read by position once the header has "
            f"proved itself, so a re-spelled one would shift every value along silently."
        )
