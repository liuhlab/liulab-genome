"""Which genes a published census judges transcription factors, and how it classifies them.

The gene half of the TF context, keyed by gene where :mod:`genome.tf.motif` is keyed
by motif. One **TF gene table** per species ships inside the package —
:mod:`genome.tf.gene.census` reads them — and every verdict in one travels with the
census that reached it: this package decides nothing about what a transcription
factor is.

**A census is keyed by Gene id stem and an annotation is not**, so
:mod:`genome.tf.gene.annotation` is what puts one into a registered **Annotation**'s own
gene ids. It crosses :meth:`~genome.io.gtf.AnnotationRegistry.resolve_gene_ids` and adds
nothing to the registry to do it: which species selects the census, and what a row of it
says, are this half's business and not the Annotation context's.

**Attribution.** Human is Lambert *et al.* 2018 (PMID 29425488),
https://humantfs.ccbr.utoronto.ca/; mouse is AnimalTFDB 4.0 (PMID 36268869),
https://guolab.wchscu.cn/. Cite the publisher when you use a census;
:meth:`~genome.tf.gene.census.CensusProvenance.attribution` renders the line.

Examples
--------
>>> from genome.tf.gene import tf_gene_table
>>> tf_gene_table("Homo sapiens").provenance.publisher
'Lambert et al. 2018'
>>> print(tf_gene_table("Mus musculus").provenance.attribution())
AnimalTFDB 4.0 (PMID 36268869) — https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/Mus_musculus_TF
"""

from genome.shipped import species_slug
from genome.tf.gene.annotation import (
    TFGene,
    TFGeneList,
    resolve_tf_genes,
    tf_gene_list,
)
from genome.tf.gene.census import (
    CENSUS_METADATA_RESOURCE,
    CENSUS_SUBDIR,
    CENSUS_SUFFIX,
    FALSE_CELL,
    TRUE_CELL,
    UNIFORM_COLUMNS,
    CensusProvenance,
    TFGeneTable,
    TFGeneTableError,
    census_metadata,
    census_species,
    tf_gene_table,
)
from genome.tf.species import NoTFCensusError

__all__ = [
    "CENSUS_METADATA_RESOURCE",
    "CENSUS_SUBDIR",
    "CENSUS_SUFFIX",
    "FALSE_CELL",
    "TRUE_CELL",
    "UNIFORM_COLUMNS",
    "CensusProvenance",
    "NoTFCensusError",
    "TFGene",
    "TFGeneList",
    "TFGeneTable",
    "TFGeneTableError",
    "census_metadata",
    "census_species",
    "resolve_tf_genes",
    "species_slug",
    "tf_gene_list",
    "tf_gene_table",
]
