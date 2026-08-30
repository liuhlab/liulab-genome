"""Which genes in another species a gene is homologous to, on Ensembl Compara's trees.

The Orthology context, and a peer of :mod:`genome.tf` rather than a part of it. A
**Homology set** is anchored to a species pair and a pinned **Release** — never to an
**Assembly** or an **Annotation** — is downloaded once into the **Data dir** as a sibling
of the assembly tree, and answers with no **Genome** open, the shape ``motif/`` already
established. It reads bulk files fetched once and never a remote API, so a pipeline built
on it does not fail intermittently the way one built on BioMart does.

**Orthology is served and never consumed** (ADR-0019). This answers a user's cross-species
question and nothing else: no **TF gene table**, no **Cofactor table**, no list this
package publishes is derived through homology, and no answer is ever silently
species-mapped. A species with no census still has none — what changed is that a user can
now cross the line themselves, deliberately, with the publisher's **Homology type** in
hand.

**Attribution.** Ensembl Compara, Herrero J *et al.*, *Database (Oxford)* 2016:bav096
(PMID 26896847), from <https://ftp.ensembl.org/pub/>. Every link, every **Homology type**
and every confidence field is the publisher's; this package computes none of them and
ranks nothing. :meth:`~genome.homology.metadata.HomologyMetadata.attribution` renders the line to
print, and ``src/genome/data/homology/ATTRIBUTION.md`` carries the same facts beside the
provenance table.

Examples
--------
>>> from genome.homology import homology_table, homology_data_dir
>>> {row.publisher for row in homology_table()}
{'Ensembl Compara'}
>>> import os
>>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
>>> homology_data_dir()
PosixPath('/scratch/liulab/homology')
>>> del os.environ["LIULAB_DATA"]
"""

from genome.homology.annotation import resolve_homologs
from genome.homology.compara import (
    COMPARA_COLUMNS,
    COMPARA_DUMP,
    COMPARA_SUBDIR,
    DEFAULT_RELEASE,
    ENSEMBL_BASE_URL,
    NULL_CELL,
    QUALITY_SCORE_COLUMNS,
    RECORD_KIND,
    ComparaFileError,
    ComparaPartitionError,
    HomologySet,
    HomologySetNotDownloadedError,
    NoHomologyPairError,
    UnknownHomologySpeciesError,
    VersionedGeneIdError,
    check_pair,
    check_release,
    check_species,
    check_stem,
    compara_url,
    homology_prepare_command,
    pair_name,
    set_dir,
    slice_filename,
    source_filename,
)
from genome.homology.metadata import (
    HOMOLOGY_METADATA_RESOURCE,
    METADATA_COLUMNS,
    HomologyMetadata,
    HomologyMetadataError,
    homology_metadata,
    homology_releases,
    homology_species,
    homology_table,
    read_metadata,
)
from genome.io.prepared import HOMOLOGY_SUBDIR, homology_data_dir

__all__ = [
    "COMPARA_COLUMNS",
    "COMPARA_DUMP",
    "COMPARA_SUBDIR",
    "DEFAULT_RELEASE",
    "ENSEMBL_BASE_URL",
    "HOMOLOGY_METADATA_RESOURCE",
    "HOMOLOGY_SUBDIR",
    "METADATA_COLUMNS",
    "NULL_CELL",
    "QUALITY_SCORE_COLUMNS",
    "RECORD_KIND",
    "ComparaFileError",
    "ComparaPartitionError",
    "HomologyMetadata",
    "HomologyMetadataError",
    "HomologySet",
    "HomologySetNotDownloadedError",
    "NoHomologyPairError",
    "UnknownHomologySpeciesError",
    "VersionedGeneIdError",
    "check_pair",
    "check_release",
    "check_species",
    "check_stem",
    "compara_url",
    "homology_data_dir",
    "homology_metadata",
    "homology_prepare_command",
    "homology_releases",
    "homology_species",
    "homology_table",
    "pair_name",
    "read_metadata",
    "resolve_homologs",
    "set_dir",
    "slice_filename",
    "source_filename",
]
