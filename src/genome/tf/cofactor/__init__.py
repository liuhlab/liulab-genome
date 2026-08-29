"""Which genes a publisher lists as transcription cofactors, and how it classifies them.

The third part of the TF context and a peer of :mod:`genome.tf.gene` rather than a part
of it: the gene half answers whether a gene is a **TF gene** and of what **DBD family**,
this one answers whether it is a **Transcription cofactor** and of what class. It is
keyed the same way, by **Gene id stem**. One **Cofactor table** per species ships inside
the package — :mod:`genome.tf.cofactor.table` reads them — and membership and
classification both travel with the publisher that reached them.

**Attribution.** Every table shipping today is AnimalTFDB 4.0 (PMID 36268869),
https://guolab.wchscu.cn/. Cite the publisher when you use one;
:meth:`~genome.tf.cofactor.table.CofactorProvenance.attribution` renders the line.

**Worm ships although no publisher has released a worm TF census**, so a worm assembly
answers here while the TF gene half raises for it. That is the publishers' shape and not
a defect; ``src/genome/data/tf_cofactor/ATTRIBUTION.md`` says so beside the files.

Examples
--------
>>> from genome.tf.cofactor import cofactor_table
>>> len(cofactor_table("Mus musculus"))
970
>>> print(cofactor_table("Caenorhabditis elegans").provenance.attribution())
AnimalTFDB 4.0 (PMID 36268869) — https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/Caenorhabditis_elegans_Cof
"""

from genome.tf.cofactor.table import (
    ANIMALTFDB,
    BOTH,
    COFACTOR_METADATA_RESOURCE,
    COFACTOR_SOURCE_METADATA_RESOURCE,
    COFACTOR_SUBDIR,
    COFACTOR_SUFFIX,
    EPIFACTORS,
    FALSE_CELL,
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
    species_slug,
)

__all__ = [
    "ANIMALTFDB",
    "BOTH",
    "COFACTOR_METADATA_RESOURCE",
    "COFACTOR_SOURCE_METADATA_RESOURCE",
    "COFACTOR_SUBDIR",
    "COFACTOR_SUFFIX",
    "EPIFACTORS",
    "FALSE_CELL",
    "SOURCES",
    "TRUE_CELL",
    "UNIFORM_COLUMNS",
    "CofactorProvenance",
    "CofactorSource",
    "CofactorTable",
    "CofactorTableError",
    "cofactor_metadata",
    "cofactor_species",
    "cofactor_table",
    "parse_cofactor_table",
    "species_slug",
]
