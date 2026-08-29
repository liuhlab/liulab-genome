"""Which foreign identifiers name a gene, and which genes a foreign identifier names.

The Xref context, and a peer of :mod:`genome.tf` rather than a part of it. An **Xref set**
is one species, one **Xref source** and one pinned **Release**: constructing it fetches the
publisher's file once into the **Data dir**, slices it to that species as a plain gzipped
TSV, and re-reads it on every construction after. It holds no coordinates, opens no
**Genome** and belongs to no **Assembly** — an identifier is a name and not a place.

Two verbs and only two (ADR-0017): :meth:`~genome.xref.xref.XrefSet.to_stems` toward the
hub and :meth:`~genome.xref.xref.XrefSet.from_stems` away from it. Nothing here converts
one foreign **Namespace** directly into another, and nothing merges two publishers — a
query reads exactly one set, so two sources that disagree are two answers rather than a
contradiction to resolve.

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

from genome.io.registration import xref_data_dir
from genome.io.results import ResolvedStems, ResolvedXrefIds
from genome.xref.alliance import ALLIANCE, AllianceFileError
from genome.xref.ids import (
    ENSEMBL,
    ENTREZ,
    HGNC,
    MGI,
    NAMESPACES,
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
    "ALLIANCE",
    "ENSEMBL",
    "ENTREZ",
    "HGNC",
    "MGI",
    "NAMESPACES",
    "UNIPROT",
    "WORMBASE",
    "AllianceFileError",
    "NamespaceNotCarriedError",
    "NoXrefSetError",
    "ResolvedStems",
    "ResolvedXrefIds",
    "XrefMetadata",
    "XrefSet",
    "XrefSetNotDownloadedError",
    "XrefTableError",
    "gene_id_stem",
    "lookup_xref",
    "normalise_id",
    "xref_data_dir",
    "xref_prepare_command",
    "xref_releases",
    "xref_set_dir",
    "xref_slice_name",
    "xref_sources",
    "xref_species",
    "xref_table",
]
