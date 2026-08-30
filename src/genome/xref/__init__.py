"""Which foreign identifiers name a gene, and which genes a foreign identifier names.

The Xref context, and a peer of :mod:`genome.tf` rather than a part of it. An **Xref set**
is one species, one **Xref source** and one pinned **Release**: constructing it fetches the
publisher's file once into the **Data dir**, slices it to that species as a plain gzipped
TSV, and re-reads it on every construction after. It holds no coordinates, opens no
**Genome** and belongs to no **Assembly** — an identifier is a name and not a place.

Two directions and only two (ADR-0017): :meth:`~genome.xref.xref.XrefSet.to_stems` toward
the hub and :meth:`~genome.xref.xref.XrefSet.from_stems` away from it — with
:meth:`~genome.xref.xref.XrefSet.match_symbols` the first of them asked about a gene
symbol, which is answered unlike an id and so has a verb of its own. Nothing here converts
one foreign **Namespace** directly into another, and nothing merges two publishers — a
query reads exactly one set, so two sources that disagree are two answers rather than a
contradiction to resolve.

Four sources ship, and they are **not equals**. :data:`~genome.xref.alliance.ALLIANCE` is
the **Default xref source** for all three species when the question is identifiers;
:data:`~genome.xref.ensembl.ENSEMBL_TSV` is selectable and pins its own numbered
**Release**, and its mapping fans out to 72 stems for one GeneID where the other's is
near-one-to-one. Which one answers is a scientific choice — read
:class:`~genome.xref.xref.XrefSet`'s ``source`` parameter before making it.

The other two carry **symbols**, which the first two do not.
:data:`~genome.xref.hgnc.HGNC_ARCHIVE` is human's, from a pinned quarterly archive file,
and is the only source that publishes previous and alias spellings *typed*;
:data:`~genome.xref.bgi.ALLIANCE_BGI` is mouse's and worm's, and carries the current
approved symbol alone. So a symbol is matched with
:meth:`~genome.xref.xref.XrefSet.match_symbols`, whose answer says on every hit which kind
of spelling matched and, on the answer as a whole, which kinds that source could not match
and why.

**A default is therefore per species and per question** (ADR-0021), and the question is
named where the source is filled in and nowhere else:
:meth:`~genome.xref.xref.XrefSet.for_symbols` is the constructor that reaches one of those
two, :meth:`~genome.xref.xref.XrefSet.for_namespace` the same fill-in for a caller holding
a **Namespace** rather than a verb, and the plain constructor reaches the identifier
default. They sit side by side on purpose — a set built for one publisher is never answered
out of another's bytes, so ``XrefSet("Homo sapiens")`` matches no symbol at all and raises
naming the source that does.

Putting an answer into your own annotation's gene ids is the call that already existed,
:meth:`~genome.io.gtf.AnnotationRegistry.resolve_gene_ids`, which is why the hub is the
**Gene id stem**: what comes out of here goes straight in there.

Examples
--------
>>> from genome.xref import XrefSet, xref_table
>>> sorted({record.species for record in xref_table()})
['Caenorhabditis elegans', 'Homo sapiens', 'Mus musculus']
>>> human = XrefSet("Homo sapiens")                                # doctest: +SKIP
>>> human.to_stems(["HGNC:11998"], "hgnc").gene_id_stems           # doctest: +SKIP
['ENSG00000141510']
"""

from genome.io.prepared import xref_data_dir
from genome.io.results import ResolvedStems, ResolvedSymbols, ResolvedXrefIds, SymbolMatch
from genome.xref.alliance import ALLIANCE, AllianceFileError
from genome.xref.bgi import ALLIANCE_BGI, BgiFileError
from genome.xref.ensembl import ENSEMBL_TSV, EnsemblTsvFileError
from genome.xref.evidence import (
    EmptyEvidenceFilterError,
    EvidenceNotRecordedError,
    normalise_evidence,
)
from genome.xref.hgnc import HGNC_ARCHIVE, HgncFileError
from genome.xref.ids import (
    ENSEMBL,
    ENTREZ,
    HGNC,
    MGI,
    NAMESPACES,
    SYMBOL,
    UNIPROT,
    WORMBASE,
    gene_id_stem,
    normalise_id,
)
from genome.xref.metadata import (
    NoXrefSetError,
    XrefMetadata,
    lookup_xref,
    xref_releases,
    xref_sources,
    xref_species,
    xref_table,
)
from genome.xref.symbols import (
    ALIAS,
    APPROVED,
    PREVIOUS,
    SYMBOL_KINDS,
    SymbolDirectionError,
    fold_symbol,
    normalise_symbol,
)
from genome.xref.xref import (
    NamespaceNotCarriedError,
    XrefSet,
    XrefSetNotDownloadedError,
    XrefTableError,
    xref_prepare_command,
    xref_set_dir,
    xref_slice_name,
)

__all__ = [
    "ALIAS",
    "ALLIANCE",
    "ALLIANCE_BGI",
    "APPROVED",
    "ENSEMBL",
    "ENSEMBL_TSV",
    "ENTREZ",
    "HGNC",
    "HGNC_ARCHIVE",
    "MGI",
    "NAMESPACES",
    "PREVIOUS",
    "SYMBOL",
    "SYMBOL_KINDS",
    "UNIPROT",
    "WORMBASE",
    "AllianceFileError",
    "BgiFileError",
    "EmptyEvidenceFilterError",
    "EnsemblTsvFileError",
    "EvidenceNotRecordedError",
    "HgncFileError",
    "NamespaceNotCarriedError",
    "NoXrefSetError",
    "ResolvedStems",
    "ResolvedSymbols",
    "ResolvedXrefIds",
    "SymbolDirectionError",
    "SymbolMatch",
    "XrefMetadata",
    "XrefSet",
    "XrefSetNotDownloadedError",
    "XrefTableError",
    "fold_symbol",
    "gene_id_stem",
    "lookup_xref",
    "normalise_evidence",
    "normalise_id",
    "normalise_symbol",
    "xref_data_dir",
    "xref_prepare_command",
    "xref_releases",
    "xref_set_dir",
    "xref_slice_name",
    "xref_sources",
    "xref_species",
    "xref_table",
]
