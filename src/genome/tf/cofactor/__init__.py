"""Which genes a publisher lists as transcription cofactors, and how it classifies them.

The third part of the TF context and a peer of :mod:`genome.tf.gene` rather than a part
of it: the gene half answers whether a gene is a **TF gene** and of what **DBD family**,
this one answers whether it is a **Transcription cofactor** and of what class. It is
keyed the same way, by **Gene id stem**. One **Cofactor table** per species ships inside
the package — :mod:`genome.tf.cofactor.table` reads them — and every classification
travels with the publisher that reached it.
:mod:`genome.tf.cofactor.annotation` puts one into a registered **Annotation**'s own gene
ids, crossing :meth:`~genome.io.annotation.registry.AnnotationRegistry.resolve_gene_ids` exactly as the
census half does and adding nothing to the registry to do it.

**Human membership is this package's own** and is the one thing here that is: the human
table is the union of two publishers' lists, 1,466 genes that neither of them publishes
(ADR-0016). Mouse and worm relay one publisher unchanged.

**Attribution.** AnimalTFDB 4.0 (PMID 36268869), https://guolab.wchscu.cn/, lists
cofactors for every species that ships; EpiFactors v2.0 (PMID 36350659),
https://epifactors.autosome.org/, lists human's beside it; and a pinned dated HGNC
monthly archive (PMID 41287213), https://www.genenames.org/, supplies the **Gene id
stem** of every gene EpiFactors names. Cite the publishers whose table you use;
:meth:`~genome.tf.cofactor.table.CofactorProvenance.attribution` renders the line.

**Worm ships although no publisher has released a worm TF census**, so a worm assembly
answers here while the TF gene half raises for it. That is the publishers' shape and not
a defect; ``src/genome/data/tf_cofactor/ATTRIBUTION.md`` says so beside the files.

Examples
--------
>>> from genome.tf.cofactor import cofactor_table
>>> len(cofactor_table("Mus musculus"))
970
>>> human = cofactor_table("Homo sapiens")
>>> len(human), sorted({row[3] for row in human.rows})
(1466, ['animaltfdb', 'both', 'epifactors'])
>>> print(cofactor_table("Caenorhabditis elegans").provenance.attribution())
AnimalTFDB 4.0 (PMID 36268869) — https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/Caenorhabditis_elegans_Cof
"""

from genome.shipped import species_slug
from genome.tf.cofactor.annotation import (
    TFCofactor,
    TFCofactorList,
    resolve_tf_cofactors,
    tf_cofactor_list,
)
from genome.tf.cofactor.table import (
    ANIMALTFDB,
    BOTH,
    CITED_SOURCES,
    COFACTOR_FORMAT,
    COFACTOR_METADATA_FORMAT,
    COFACTOR_METADATA_RESOURCE,
    COFACTOR_SOURCE_METADATA_FORMAT,
    COFACTOR_SOURCE_METADATA_RESOURCE,
    COFACTOR_SUBDIR,
    COFACTOR_SUFFIX,
    EPIFACTORS,
    FALSE_CELL,
    HGNC,
    SOURCES,
    TRUE_CELL,
    UNIFORM_COLUMNS,
    CofactorProvenance,
    CofactorSource,
    CofactorTable,
    CofactorTableError,
    cofactor_metadata,
    cofactor_species,
    cofactor_table,
    parse_cofactor_table,
)
from genome.tf.species import NoCofactorTableError

__all__ = [
    "ANIMALTFDB",
    "BOTH",
    "CITED_SOURCES",
    "COFACTOR_FORMAT",
    "COFACTOR_METADATA_FORMAT",
    "COFACTOR_METADATA_RESOURCE",
    "COFACTOR_SOURCE_METADATA_FORMAT",
    "COFACTOR_SOURCE_METADATA_RESOURCE",
    "COFACTOR_SUBDIR",
    "COFACTOR_SUFFIX",
    "EPIFACTORS",
    "FALSE_CELL",
    "HGNC",
    "SOURCES",
    "TRUE_CELL",
    "UNIFORM_COLUMNS",
    "CofactorProvenance",
    "CofactorSource",
    "CofactorTable",
    "CofactorTableError",
    "NoCofactorTableError",
    "TFCofactor",
    "TFCofactorList",
    "cofactor_metadata",
    "cofactor_species",
    "cofactor_table",
    "parse_cofactor_table",
    "resolve_tf_cofactors",
    "species_slug",
    "tf_cofactor_list",
]
