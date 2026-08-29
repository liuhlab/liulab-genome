"""Transcription factors and the motifs they bind.

Two halves, kept apart because they are keyed differently. :mod:`genome.tf.gene`
is about the *gene* — which genes a published census judges transcription
factors, and the DNA-binding-domain family it classifies each one under — and so
it is keyed by gene. :mod:`genome.tf.motif` is about the *sequence a factor
recognises* — the count matrix and where it occurs in an assembly — and so it is
keyed by motif. The mapping between the two is many-to-many, which is why neither
half owns it. Cofactors are out of scope and nothing here answers for them.

The join itself is :mod:`genome.tf.link`, which lives here rather than in either
half because it imports both and neither imports it. It is re-exported from this
package for that reason: it belongs to the pair, so it is addressed at the level
the pair is.

Examples
--------
>>> from genome.tf import motif_links
>>> motif_links("TP53", "Homo sapiens").motif_ids
('MA0106.3',)
"""

from genome.tf.link import (
    COMPLEX,
    LINK_COLUMNS,
    LINK_SUBDIR,
    LINK_SUFFIX,
    LINK_TAX_GROUP,
    MONOMER,
    RELEASE_PREFIX,
    VALUE_SEPARATOR,
    GeneNotAssessedError,
    MotifLink,
    MotifLinks,
    MotifLinkTable,
    MotifLinkTableError,
    NoMotifLinkTableError,
    VersionedGeneIdError,
    motif_link_table,
    motif_links,
    parse_motif_link_table,
    shipped_link_tables,
)

__all__ = [
    "COMPLEX",
    "LINK_COLUMNS",
    "LINK_SUBDIR",
    "LINK_SUFFIX",
    "LINK_TAX_GROUP",
    "MONOMER",
    "RELEASE_PREFIX",
    "VALUE_SEPARATOR",
    "GeneNotAssessedError",
    "MotifLink",
    "MotifLinkTable",
    "MotifLinkTableError",
    "MotifLinks",
    "NoMotifLinkTableError",
    "VersionedGeneIdError",
    "motif_link_table",
    "motif_links",
    "parse_motif_link_table",
    "shipped_link_tables",
]
