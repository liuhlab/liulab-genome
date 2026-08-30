"""Annotations — a GTF filed under one **Assembly**, and everything asked of it after.

I/O boundary package. A reference assembly may carry several gene annotations (GENCODE,
RefSeq, WormBase, …). Each is registered under a **Registered name** and lives in its own
directory beside the assembly's sequence files::

    <LIULAB_DATA>/genome/<assembly>/gtf/<name>/
        <name>.gtf          # the annotation, kept decompressed
        <name>.db           # the gffutils SQLite database built from it
        .completion.json    # the record saying all of that finished
        .work/              # the disposable working area a fetch downloads into

:class:`~genome.annotation.registry.AnnotationRegistry` is the way in, and it is **one
class with one interface**. Bound once to one assembly — its name, its **Assembly dir** and
its ``chrom.sizes`` — it holds every annotation that assembly has and answers everything
about them. Its implementation is spread over four modules and it calls across them:

- :mod:`~genome.annotation.registration` puts an annotation on disk: the fetch, the
  placement, the **Chromosome** check, the repair-command strings, the **Completion
  marker**, the **Merged annotation** a **Chimera** build derives, and the two registrars
  addressed by assembly name.
- :mod:`~genome.annotation.registry` holds the class itself, the three scans and the
  **Default annotation** rule it settles at construction, and the by-assembly-name
  questions ``genome annotations``, ``genome gene-list`` and ``genome gene-categories``
  ask.
- :mod:`~genome.annotation.stems` resolves a **Gene id stem** against an annotation's
  own gene ids — the seam the Xref, Orthology and TF contexts all cross, and the only one.
- :mod:`~genome.annotation.database` is the `gffutils` adapter: the build, and the read
  that yields gene ids a row at a time. Nothing else in the package imports the library.

Two more modules sit beside those four because they are about the same thing and nothing
else: :mod:`~genome.annotation.metadata` is the curated table of what the lab supports for
each assembly, keyed by assembly plus **Registered name**, and
:mod:`~genome.annotation.curated` reads the shipped **Curated gene list** — which of an
annotation's genes are in a hand-curated category (ADR-0011).

**Nothing here knows what a caller's gene ids are a list of.** Which species selects a
shipped file, and what a row of one says, are facts this package has no stake in, so it
holds no import of any of them and gains no method when a fourth topic arrives.

What a registration answers with is re-exported here beside what produces it, so a caller
holding one imports it from the package rather than from the module that happens to define
it (ADR-0022).

Examples
--------
>>> from pathlib import Path
>>> from genome.annotation import annotation_dir
>>> annotation_dir(Path("/data/genome/sacCer3"), "ensgene_v101").name
'ensgene_v101'
"""

from genome.annotation.curated import (
    GeneCategoryNotDeclaredError,
    NoGeneCategoriesError,
)
from genome.annotation.metadata import AnnotationMetadata
from genome.annotation.registration import (
    UNCHECKED_CALLER_OVERRIDE,
    UNCHECKED_NO_CHROM_SIZES,
    ChromosomeMismatchError,
    GtfAnnotation,
    MergeSource,
    RegisteredAnnotation,
    annotation_dir,
    annotation_register_command,
    chromosome_check_summary,
    discard_merged_annotation,
    register_annotation,
    register_gtf,
    register_merged_gtf,
)

# The three names the registry absorbed — `list_annotations`, `list_broken_annotations`
# and `default_annotation` — are deliberately absent. Each is reached through
# `AnnotationRegistry`, nothing in the package calls them by name, and what a caller wants
# from an assembly's annotations is the registry rather than the three steps it takes
# internally. They stay importable from `genome.annotation.registry` for the tests that
# hold each rule on its own.
from genome.annotation.registry import (
    AnnotationNotRegisteredError,
    AnnotationRegistry,
    AnnotationStatus,
    AnnotationStatusRow,
    BrokenAnnotation,
    GeneList,
    GeneListSource,
    annotation_status,
    gene_list,
    gene_lists,
)
from genome.annotation.stems import NoGeneFeaturesError, ResolvedGeneIds

__all__ = [
    "UNCHECKED_CALLER_OVERRIDE",
    "UNCHECKED_NO_CHROM_SIZES",
    "AnnotationMetadata",
    "AnnotationNotRegisteredError",
    "AnnotationRegistry",
    "AnnotationStatus",
    "AnnotationStatusRow",
    "BrokenAnnotation",
    "ChromosomeMismatchError",
    "GeneCategoryNotDeclaredError",
    "GeneList",
    "GeneListSource",
    "GtfAnnotation",
    "MergeSource",
    "NoGeneCategoriesError",
    "NoGeneFeaturesError",
    "RegisteredAnnotation",
    "ResolvedGeneIds",
    "annotation_dir",
    "annotation_register_command",
    "annotation_status",
    "chromosome_check_summary",
    "discard_merged_annotation",
    "gene_list",
    "gene_lists",
    "register_annotation",
    "register_gtf",
    "register_merged_gtf",
]
