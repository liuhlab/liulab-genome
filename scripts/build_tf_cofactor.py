#!/usr/bin/env python
r"""Rebuild one shipped **Cofactor table** from its publishers' own files.

Run by hand, never at install time and never in CI — the publishers' files are
downloads CI cannot make, which is why the shipped tables are guarded by
``tests/test_tf_cofactor.py`` instead. This script lives outside the package for
the same reason: the wheel carries the table and not the tooling that made it.
**Every curation rule below lives here and none of it lives in the wheel.**

Attribution
-----------
``animaltfdb4_mouse``, ``animaltfdb4_worm`` and the AnimalTFDB half of
``union_human`` — Shen *et al.*, "AnimalTFDB 4.0", *Nucleic Acids Research*
51(D1):D39-D45, 2023 (PMID 36268869). Two files per species, both from
https://guolab.wchscu.cn/: the per-species cofactor list
(``.../Cof_list_final/<Genus_species>_Cof``), which is membership and family, and
the family summary (``.../cof_info_summary.tsv``), which is the category each
family sits in. The all-species bulk file is deliberately not read: the per-species
download is the one whose row count the summary reconciles against.

``union_human`` reads two more. EpiFactors v2.0 — Marakulina *et al.*,
"EpiFactors 2022", *Nucleic Acids Research* 51(D1):D564-D570, 2023 (PMID 36350659),
one CSV from https://epifactors.autosome.org/ — contributes membership and its own
function, target, modification and complex vocabularies. A **dated** HGNC monthly
archive — Seal *et al.*, "Genenames.org: the HGNC and PGNC resources in 2026",
*Nucleic Acids Research* 54(D1):D1098-D1107, 2026 (PMID 41287213) — contributes
identifiers and no membership: it is what makes an EpiFactors row readable as a
**Gene id stem** at all, and it fixes the symbol.

Cite the publishers whose table you use.

The four curation rules
-----------------------
*Membership is unioned; classification is not.* A gene either publisher lists is in
the human table. Its AnimalTFDB columns are filled only if AnimalTFDB listed it and
its EpiFactors columns only if EpiFactors did, and nothing is inferred across the
two in either direction for any pair of values — a row saying ``both`` asserts
agreement on membership and on nothing else, and ADR-0014 applies here unchanged.
Human membership is therefore this package's own verdict and nobody else's
(ADR-0016); mouse and worm relay one publisher unaltered.

*EpiFactors joins to Ensembl through its HGNC id, never through its symbol.* All 801
of its rows carry one and all 801 reach the archive, where matching on the symbol
would mis-key or drop the 31 rows EpiFactors still spells by a name HGNC has since
retired — ``ACINU`` for ``ACIN1``, ``ARNTL`` for ``BMAL1``, ``C11orf30`` for
``EMSY``. The archive is a **pinned dated file** and never the rolling current one,
so the 442 stems only it can supply are reproducible; its dates are irregular, so
the pin names a file read from the archive listing rather than one built from a date.

*A gene EpiFactors gives two rows collapses into one, its cells unioned.* Five do —
``ALKBH1``, ``HSPA1A``, ``HSPA1B``, ``NAT10`` and ``PTBP1`` — and they are not
duplicate rows: one carries a histone-modification annotation and the other an
RNA-modification one. Values are deduplicated within a cell as well as across the
two rows, since a table is one row per stem. The cost is real and is stated in
``ATTRIBUTION.md``: for those five the pairing between a function and its own
modification is lost.

*The symbol is HGNC's current approved spelling* for every human row, reached
through the id and never through a name — the EpiFactors rows through their HGNC id,
the AnimalTFDB rows through the stem that is already the table's key. Mouse and worm
keep their publisher's own spelling, where HGNC does not apply.

What it does
------------
Reads its publishers' files, holds their headers to the columns each release
published, joins each gene's category onto its family, unions the two human lists,
and writes a gzipped TSV plus one row in the species provenance table and one per
source in the ragged one. Four properties are the point:

*It fails loudly on an unrecognised header.* A publisher who re-spells ``Family``
would otherwise have that column silently dropped, and nothing downstream would
say so. Every published name is listed here, so a re-spelling is one missing name
and one unexpected name, and both are named in the error. HGNC's archive is held to
:func:`read_named_columns` instead — it publishes fifty-odd columns and gains more
between releases, so a new one there is the publisher growing and only a *missing*
one is the re-spelling worth breaking on.

*It fails loudly on a family with no category.* The publisher's own two files
spell some families differently, so the join goes through the hand-written
:data:`FAMILY_SPELLINGS` below. A family that survives that map and is still not in
the summary breaks the build rather than blanking a column, which is what a future
release renaming a family must do.

*It fails loudly on an identifier it cannot resolve.* An EpiFactors row whose HGNC
id the archive does not name, one the archive gives no Ensembl id, a stem HGNC names
under two symbols, and a published value that already contains
:data:`VALUE_SEPARATOR` are each a refusal naming the rows, never a blanked cell.

*It writes byte-stable output.* Rows keep the publishers' own order — AnimalTFDB's
list, then the genes only EpiFactors names, in its order — cells are written with no
quoting, the line terminator is ``\n``, and gzip is given ``mtime=0``, so re-running
on unchanged inputs produces no diff.

Usage
-----
``python scripts/build_tf_cofactor.py <table> <the species cofactor list>
<cof_info_summary.tsv>``, where ``<table>`` names one of :data:`RECIPES`.

``union_human`` needs its other two publishers as well::

    python scripts/build_tf_cofactor.py union_human Homo_sapiens_Cof \
        cof_info_summary.tsv --epifactors EpiGenes_main.csv \
        --hgnc hgnc_complete_set_2026-08-07.txt

Add ``--data-dir`` to write somewhere other than ``src/genome/data/tf_cofactor``.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# The one function this script shares with the package rather than copying. Everything
# else below is the writer's own declaration of what it writes.
from genome.tf.gene.census import species_slug

#: Where the shipped tables live, relative to the repository root.
DATA_SUBDIR = Path("src") / "genome" / "data" / "tf_cofactor"

#: What one table file is called: the species slug, then this. Kept in step with
#: ``genome.tf.cofactor.table.COFACTOR_SUFFIX``, which is what enumerates them.
COFACTOR_SUFFIX = ".cofactor_table.tsv.gz"

#: The provenance table keyed by species, and its columns in table order.
METADATA_FILE = "cofactor_metadata.tsv"
METADATA_COLUMNS = ("species", "ncbi_taxid", "file", "sha256")

#: The provenance table keyed by species *and* source, and its columns in table order.
#: Two tables and not one, because one row cannot describe a file built from three
#: publishers and joining them positionally inside a cell breaks quietly. This one is
#: deliberately ragged: one row per species here, three for a species built from three.
SOURCE_METADATA_FILE = "cofactor_source_metadata.tsv"
SOURCE_METADATA_COLUMNS = ("species", "source", "publisher", "version", "pubmed_id", "source_url")

#: The four columns every shipped table carries under the same name, in the same place.
UNIFORM_COLUMNS = ("gene_id_stem", "symbol", "is_cofactor", "source")

#: AnimalTFDB's own two columns, under a namespaced spelling, after the uniform four.
#: A second publisher's group is more columns here and one more source row — never a
#: change to the format.
ANIMALTFDB_COLUMNS = ("animaltfdb_family", "animaltfdb_category")

#: How a shipped table spells its cofactor flag.
TRUE_CELL, FALSE_CELL = "yes", "no"

#: The values of the uniform ``source`` column, from the closed vocabulary
#: ``genome.tf.cofactor`` validates on read. ``both`` says two publishers listed the
#: gene — agreement on membership only, and never on how either of them classified it.
ANIMALTFDB, EPIFACTORS, BOTH = "animaltfdb", "epifactors", "both"

#: A source that supplies identifiers and lists no gene. It never spells a row's
#: ``source`` and appears only in the provenance table, which is a wider vocabulary for
#: exactly this reason: human's 442 EpiFactors-only stems exist only because HGNC said
#: so, and a source that earns a stem earns a citation.
HGNC = "hgnc"

#: How a cell spells more than one value. The separator every multi-valued cell in this
#: package already uses — ``interpro_ids`` in a **TF gene table**, a **Motif link**'s
#: partners, ``genome.tf.link.VALUE_SEPARATOR`` — so that a caller never has to remember
#: which column uses which. :func:`published_values` refuses a publisher's value that
#: already contains it rather than writing a cell that splits into the wrong values.
VALUE_SEPARATOR = ";"

#: The header the per-species cofactor list publishes, in order.
GENE_LIST_COLUMNS = ("Species", "Symbol", "Ensembl", "Family", "Entrez_ID")

#: The header EpiFactors' main gene table publishes, in order. Held whole, as
#: AnimalTFDB's two files are: it is one curated release rather than a rolling file, so
#: a column that is not on this list is a re-spelling and worth breaking on.
EPIFACTORS_GENE_COLUMNS = (
    "Id",
    "HGNC_symbol",
    "Status",
    "HGNC_ID",
    "HGNC_name",
    "GeneID",
    "UniProt_AC",
    "UniProt_ID",
    "Domain",
    "MGI_symbol",
    "MGI_ID",
    "UniProt_AC_Mm",
    "UniProt_ID_Mm",
    "GeneTag",
    "GeneDesc",
    "Function",
    "Modification",
    "PMID_function",
    "Complex_name",
    "Target",
    "Specific_target",
    "Product",
    "UniProt_ID_target",
    "PMID_target",
    "Comment",
)

#: Which of EpiFactors' published columns ship, and under what namespaced name. The
#: order here is the order they appear in the shipped table, which is this package's
#: and not the publisher's — the publisher files ``Target`` after ``Complex_name``.
EPIFACTORS_FIELDS: Mapping[str, str] = {
    "Function": "epifactors_function",
    "Target": "epifactors_target",
    "Modification": "epifactors_modification",
    "Complex_name": "epifactors_complex_name",
}

#: EpiFactors' own four columns, under a namespaced spelling, after AnimalTFDB's two.
EPIFACTORS_COLUMNS = tuple(EPIFACTORS_FIELDS.values())

#: How EpiFactors spells *nothing recorded*: the literal ``#``, on every one of those
#: four columns and never a blank or an ``NA``. It reaches the shipped table as an empty
#: cell, which is how every table here says that a publisher recorded nothing.
EPIFACTORS_NOTHING = "#"

#: How EpiFactors separates two values inside one cell: a comma **and a space**, never a
#: bare comma. Two of its complex names carry a comma inside themselves —
#: ``COMPASS-like MLL1,2`` and ``COMPASS-like MLL3,4`` — so splitting on the comma alone
#: invents complexes called ``2`` and ``4``. Its other three columns are safe either way
#: and are split by the same rule, because one rule per file is the readable one.
EPIFACTORS_SEPARATOR = ", "

#: How an EpiFactors row spells the HGNC id the archive keys on: the bare number the
#: publisher records, under HGNC's own prefix.
HGNC_ID_PREFIX = "HGNC:"

#: The three columns this script reads out of HGNC's monthly archive: the id EpiFactors
#: joins on, the approved symbol every human row's ``symbol`` is, and the Ensembl gene id
#: that is the **Gene id stem**. The archive publishes fifty-odd more and is not held to
#: a whole header — see :func:`read_named_columns`.
HGNC_COLUMNS = ("hgnc_id", "symbol", "ensembl_gene_id")

#: The header the family summary publishes, in order. ``categorie`` is the publisher's
#: own spelling and is written here as published.
SUMMARY_COLUMNS = ("Scientific_name_paste", "family", "family_count", "categorie")

#: The five spellings that differ between AnimalTFDB's gene list and its own family
#: summary, mapping the gene list's spelling to the summary's. Hand-written, because
#: the publisher provides no crosswalk between its two files.
#:
#: Each rule is proved by the publisher's own arithmetic, measured on the 2026-08-29
#: downloads: filter the summary to the species, map every family the gene list uses
#: through this table, and sum ``family_count`` over the *distinct* summary families
#: that result. That sum equals the number of genes the list carries, exactly:
#:
#: ===================  =====  =====================  ====================  ==========
#: species              genes  families in gene list  distinct summary keys  summed
#: ===================  =====  =====================  ====================  ==========
#: Homo sapiens          1024                     85                     82        1024
#: Mus musculus           970                     84                     81         970
#: C. elegans             317                     57                     56         317
#: ===================  =====  =====================  ====================  ==========
#:
#: Mouse needs all five rules; worm needs only ``Other_Co-activator/repressors`` and
#: ``MYC``. :func:`reconcile` re-runs that arithmetic on every build, so a release that
#: re-spells a family, moves a gene between families or re-counts one is a failed build
#: rather than a quietly wrong category column.
FAMILY_SPELLINGS: Mapping[str, str] = {
    "Lysine methyltransferase": "Lysine methyltransferase family",
    "Histone lysine methyltransferase": "Lysine methyltransferase family",
    "Other_Co-activator/repressors": "Other_Co-activator_repressors",
    "MYB": "Others",
    "MYC": "Others",
}


class CofactorSourceError(ValueError):
    """The cofactor table asked for cannot be built.

    Either the publisher's files are not the files the recipe was written against, or
    the two of them no longer reconcile with each other.
    """


@dataclass(frozen=True)
class Source:
    """One source's row of the ragged provenance table: what to cite for part of a table.

    A source that only supplies identifiers gets one of these too. HGNC lists no gene
    and so never spells a row's ``source``, but the stems of 442 human genes exist
    because it said so, which is a contribution to cite rather than an implementation
    detail to bury.

    Attributes
    ----------
    source : str
        Which source the row is about, spelled as the shipped ``source`` column spells
        one where the two vocabularies meet. Never ``both``, which describes a row of a
        table rather than anybody who published one.
    publisher : str
        Who published it, and who is to be cited for it.
    version : str
        The publisher's own release identifier — a version for a released database,
        the archive's own date for a dated file.
    pubmed_id : int
        The PubMed id of the paper to cite.
    source_url : str
        Where that publisher's own file was downloaded from.
    """

    source: str
    publisher: str
    version: str
    pubmed_id: int
    source_url: str


@dataclass(frozen=True)
class Recipe:
    """How one species' publishers' downloads become a shipped **Cofactor table**.

    Attributes
    ----------
    species : str
        The species as the assembly metadata table spells it, e.g. ``"Mus musculus"``.
        Its slug is the table's file name.
    ncbi_taxid : int
        NCBI taxonomy id for that species.
    animaltfdb_species : str
        The publisher's own spelling of the species, which both of its files key on —
        ``Species`` in the gene list and ``Scientific_name_paste`` in the summary.
    publisher : str
        Who published the AnimalTFDB list every recipe starts from.
    version : str
        That publisher's own release identifier.
    pubmed_id : int
        The PubMed id of the paper to cite for it.
    source_url : str
        Where the per-species cofactor list is downloaded from. The one URL that
        provenance row carries, because it is the one that names this species; the
        summary file is the same release and is named in ``ATTRIBUTION.md``.
    summary_url : str
        Where the family summary is downloaded from.
    unions : tuple of Source
        Every further source this species' table is built from, in the order their
        provenance rows are written. Empty for a species one publisher answers for;
        EpiFactors and HGNC for human, whose membership is a union and therefore this
        package's own (ADR-0016).
    """

    species: str
    ncbi_taxid: int
    animaltfdb_species: str
    publisher: str
    version: str
    pubmed_id: int
    source_url: str
    summary_url: str
    unions: tuple[Source, ...] = ()

    @property
    def slug(self) -> str:
        """Return the species slug this table's file is named by."""
        return species_slug(self.species)

    @property
    def file_name(self) -> str:
        """Return the table's file name."""
        return f"{self.slug}{COFACTOR_SUFFIX}"

    @property
    def sources(self) -> tuple[Source, ...]:
        """Return every source this table is built from, AnimalTFDB's first."""
        animaltfdb = Source(
            source=ANIMALTFDB,
            publisher=self.publisher,
            version=self.version,
            pubmed_id=self.pubmed_id,
            source_url=self.source_url,
        )
        return (animaltfdb, *self.unions)

    @property
    def unions_epifactors(self) -> bool:
        """Return whether this table unions EpiFactors' list in beside AnimalTFDB's."""
        return any(source.source == EPIFACTORS for source in self.unions)


#: The summary file, which every recipe reads and which names all 186 species it covers.
SUMMARY_URL = "https://guolab.wchscu.cn/AnimalTFDB4_static/download/cof_info_summary.tsv"

#: The per-species cofactor list, whose name is the publisher's own species spelling.
_COF_LIST_URL = "https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/{species}_Cof"

ANIMALTFDB4_MOUSE = Recipe(
    species="Mus musculus",
    ncbi_taxid=10090,
    animaltfdb_species="Mus_musculus",
    publisher="AnimalTFDB",
    version="4.0",
    pubmed_id=36268869,
    source_url=_COF_LIST_URL.format(species="Mus_musculus"),
    summary_url=SUMMARY_URL,
)

ANIMALTFDB4_WORM = Recipe(
    species="Caenorhabditis elegans",
    ncbi_taxid=6239,
    animaltfdb_species="Caenorhabditis_elegans",
    publisher="AnimalTFDB",
    version="4.0",
    pubmed_id=36268869,
    source_url=_COF_LIST_URL.format(species="Caenorhabditis_elegans"),
    summary_url=SUMMARY_URL,
)

#: EpiFactors' own release, read from its versioned main gene table rather than from
#: whatever ``latest`` resolves to, so a rebuild names the release it read.
EPIFACTORS_URL = "https://epifactors.autosome.org/public_data/v2.0/EpiGenes_main.csv"

#: The HGNC archive this build is pinned to: a **dated monthly file**, never the rolling
#: current one, because the stems of 442 human genes come from nowhere else and a moving
#: crosswalk makes them unreproducible. The archive's dates are irregular, so the pin is
#: a file name read from the archive listing —
#: ``https://storage.googleapis.com/storage/v1/b/public-download-files/o?prefix=hgnc/archive/archive/monthly/tsv/hgnc_complete_set_``
#: — and never one constructed from a date. Changing it is a re-run and a reviewable diff.
HGNC_ARCHIVE_DATE = "2026-08-07"
HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/"
    f"hgnc_complete_set_{HGNC_ARCHIVE_DATE}.txt"
)

UNION_HUMAN = Recipe(
    species="Homo sapiens",
    ncbi_taxid=9606,
    animaltfdb_species="Homo_sapiens",
    publisher="AnimalTFDB",
    version="4.0",
    pubmed_id=36268869,
    source_url=_COF_LIST_URL.format(species="Homo_sapiens"),
    summary_url=SUMMARY_URL,
    unions=(
        Source(
            source=EPIFACTORS,
            publisher="EpiFactors",
            version="v2.0",
            pubmed_id=36350659,
            source_url=EPIFACTORS_URL,
        ),
        Source(
            source=HGNC,
            publisher="HGNC",
            version=HGNC_ARCHIVE_DATE,
            pubmed_id=41287213,
            source_url=HGNC_URL,
        ),
    ),
)

#: Every table this script knows how to build, by the name given on the command line.
#: Named for what makes the table rather than for the species alone: mouse and worm relay
#: one publisher, and human is a union this package publishes.
RECIPES: dict[str, Recipe] = {
    "animaltfdb4_mouse": ANIMALTFDB4_MOUSE,
    "animaltfdb4_worm": ANIMALTFDB4_WORM,
    "union_human": UNION_HUMAN,
}


def read_source(path: Path, columns: tuple[str, ...], *, separator: str = "\t") -> pd.DataFrame:
    """Read one of the publisher's files, refusing it unless its header is the one expected.

    Every cell is text: the publisher's own missing-value spellings are not read as
    NaN here, since which of them a column uses is a fact about that column worth
    seeing rather than one pandas should guess at. EpiFactors publishes a CSV whose
    cells are quoted, hence ``separator``; the reader is otherwise the same one.
    """
    frame = pd.read_csv(path, sep=separator, dtype=str, keep_default_na=False)
    frame.columns = pd.Index([str(name).strip() for name in frame.columns])
    check_header(tuple(frame.columns), columns, origin=str(path))
    return frame


def read_named_columns(path: Path, columns: tuple[str, ...]) -> pd.DataFrame:
    """Read the named columns of a file that publishes many more, refusing a missing one.

    The counterpart of :func:`read_source` for a file this script reads a handful of
    columns out of. HGNC's monthly archive publishes fifty-odd and gains more between
    releases, so holding it to a whole header would break the build every time the
    publisher added something this script does not read. A *missing* column is the
    re-spelling worth breaking on, and only that is checked.

    Raises
    ------
    CofactorSourceError
        If the file does not carry one of ``columns``.
    """
    header = tuple(str(name).strip() for name in pd.read_csv(path, sep="\t", nrows=0).columns)
    missing = [name for name in columns if name not in header]
    if missing:
        raise CofactorSourceError(
            f"{path} does not carry the columns {missing}, which this script reads out of it. "
            f"It publishes {len(header)} columns and is allowed to publish more than these — a "
            f"missing one is the publisher re-spelling something, so update this script's column "
            f"list to the release you are building from and say in the commit which names moved."
        )
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, usecols=list(columns))
    frame.columns = pd.Index([str(name).strip() for name in frame.columns])
    return frame


def check_header(seen: tuple[str, ...], expected: tuple[str, ...], *, origin: str) -> None:
    """Hold ``seen`` to the columns this script was written against, naming every difference.

    Raises
    ------
    CofactorSourceError
        If a published column is missing or an unexpected one appears. A publisher
        who re-spells a column shows up as both, and both are named — dropping the
        column silently is the outcome this exists to prevent.
    """
    missing = [name for name in expected if name not in seen]
    unexpected = [name for name in seen if name not in expected]
    if not missing and not unexpected:
        return
    parts = []
    if missing:
        parts.append(f"columns this script expects and {origin} does not carry: {missing}")
    if unexpected:
        parts.append(f"columns {origin} carries and this script does not name: {unexpected}")
    raise CofactorSourceError(
        f"{origin} is not the file this script was written against — "
        + "; ".join(parts)
        + ". A publisher who re-spells a column would otherwise have it silently dropped, "
        "so update this script's published column lists to match the release you are "
        "building from, and say in the commit which columns moved."
    )


def species_rows(frame: pd.DataFrame, column: str, recipe: Recipe, *, origin: str) -> pd.DataFrame:
    """Return the rows of ``frame`` the publisher keys to this recipe's species.

    Raises
    ------
    CofactorSourceError
        If the file names that species nowhere, which is what building from the wrong
        species' download, or from a release that re-spelled the species, looks like.
    """
    kept = frame[frame[column].str.strip() == recipe.animaltfdb_species]
    if kept.empty:
        raise CofactorSourceError(
            f"{origin} carries no row whose {column!r} is {recipe.animaltfdb_species!r}, so "
            f"nothing in it describes {recipe.species}. Check that this is the file for that "
            f"species, and that the release still spells the species that way."
        )
    return kept


def category_by_family(summary: pd.DataFrame, recipe: Recipe, *, origin: str) -> dict[str, str]:
    """Return this species' family-to-category map, as the summary spells its families.

    Raises
    ------
    CofactorSourceError
        If the summary names one family twice for this species, which would make the
        category depend on which row won.
    """
    rows = species_rows(summary, "Scientific_name_paste", recipe, origin=origin)
    families = rows["family"].str.strip()
    repeated = sorted(families[families.duplicated()].unique())
    if repeated:
        raise CofactorSourceError(
            f"{origin} names the families {repeated} more than once for {recipe.species}, so "
            f"which category each of them sits in depends on which row is read last. Reconcile "
            f"the publisher's file before building from it."
        )
    return dict(zip(families, rows["categorie"].str.strip(), strict=True))


def counts_by_family(summary: pd.DataFrame, recipe: Recipe, *, origin: str) -> dict[str, int]:
    """Return this species' per-family gene counts, as the summary publishes them."""
    rows = species_rows(summary, "Scientific_name_paste", recipe, origin=origin)
    try:
        counts = [int(count) for count in rows["family_count"]]
    except ValueError as error:
        raise CofactorSourceError(
            f"{origin} holds a 'family_count' that is not a number for {recipe.species}. "
            f"That column is what proves the family spelling map right, so it has to be "
            f"countable — reconcile the publisher's file before building from it."
        ) from error
    return dict(zip(rows["family"].str.strip(), counts, strict=True))


def map_family(family: str) -> str:
    """Return the summary's spelling of one gene-list family, mapped where the two differ."""
    return FAMILY_SPELLINGS.get(family, family)


def check_families(families: Sequence[str], category: Mapping[str, str], *, origin: str) -> None:
    """Refuse to build when a family survives the spelling map with no category.

    Raises
    ------
    CofactorSourceError
        If any family the gene list uses reaches no row of the summary. Blanking the
        category column for those genes is the silent outcome this exists to prevent:
        a release that renames a family must break the build.
    """
    unmapped = sorted({family for family in families if map_family(family) not in category})
    if not unmapped:
        return
    raise CofactorSourceError(
        f"{origin} classifies genes under the families {unmapped}, and the family summary "
        f"gives no category for them even after this script's spelling map. Every family must "
        f"reach a category, or those genes ship with the column blank and nothing says so. Add "
        f"the release's new spelling to FAMILY_SPELLINGS in this script — with the arithmetic "
        f"that proves it, as the five already there carry — or fix the family it should map to."
    )


def reconcile(
    families: Sequence[str], counts: Mapping[str, int], *, genes: int, origin: str
) -> None:
    """Re-run the arithmetic that proves the family spelling map right, on this build's own files.

    Summing the summary's ``family_count`` over the *distinct* summary families the
    gene list maps onto must give exactly the number of genes the list carries. It
    does for all three species AnimalTFDB is read for here, which is what makes five
    hand-written spellings a reading of the publisher's files rather than a guess at
    them.

    Raises
    ------
    CofactorSourceError
        If the two numbers disagree. A wrong spelling shows up here even when every
        family happens to find *some* category, so this is the check that catches a
        family mapped onto the wrong neighbour rather than onto nothing.
    """
    keys = sorted({map_family(family) for family in families})
    total = sum(counts[key] for key in keys)
    if total == genes:
        return
    raise CofactorSourceError(
        f"{origin} lists {genes} genes, and the family summary counts {total} across the "
        f"{len(keys)} families they map onto. Those two reconciling exactly is what proves "
        f"this script's family spelling map is reading the publisher's files rather than "
        f"guessing at them, so a difference means a family is mapped onto the wrong summary "
        f"row, or the release has re-counted one. Re-measure before changing FAMILY_SPELLINGS."
    )


def animaltfdb_rows(
    genes: pd.DataFrame, summary: pd.DataFrame, recipe: Recipe, *, origin: str
) -> pd.DataFrame:
    """Return AnimalTFDB's own contribution — stem, symbol, family, category — in its row order.

    The publisher's order is kept rather than sorted by **Gene id stem**. It is already
    deterministic, which is all byte-stability needs, and the file is ordered by family
    — sorting it scatters exactly the runs that compress.

    The family is the *gene list's* own spelling, never the summary's: the spelling map
    exists to find each gene's category and not to re-spell what the publisher
    classified it under.
    """
    # Positional from here on: the category column below is built as a list, so the
    # kept rows are renumbered rather than carrying the publisher file's own index.
    rows = species_rows(genes, "Species", recipe, origin=origin).reset_index(drop=True)
    family = rows["Family"].str.strip()
    category = category_by_family(summary, recipe, origin=origin)
    check_families(family.tolist(), category, origin=origin)
    reconcile(
        family.tolist(),
        counts_by_family(summary, recipe, origin=origin),
        genes=len(rows),
        origin=origin,
    )
    return pd.DataFrame(
        {
            "gene_id_stem": rows["Ensembl"].str.strip(),
            "symbol": rows["Symbol"].str.strip(),
            "animaltfdb_family": family,
            "animaltfdb_category": [category[map_family(name)] for name in family],
        }
    )


def build_table(
    genes: pd.DataFrame, summary: pd.DataFrame, recipe: Recipe, *, origin: str
) -> pd.DataFrame:
    """Return the shipped table for a species one publisher answers for: the uniform four then AnimalTFDB's two.

    Membership and classification are both AnimalTFDB's here, and the symbol is its own
    spelling: HGNC names human genes and no others, so mouse and worm are relayed
    unaltered rather than corrected against a nomenclature that does not cover them.
    """
    animaltfdb = animaltfdb_rows(genes, summary, recipe, origin=origin)
    table = pd.DataFrame(
        {
            "gene_id_stem": animaltfdb["gene_id_stem"],
            "symbol": animaltfdb["symbol"],
            # Every gene AnimalTFDB lists is one it accepts: it publishes no rejected
            # set, so the flag is supplied and reads `yes` on every row. It is kept
            # anyway — dropping it would make presence in the file the verdict, and a
            # future source could then not record a rejection without a format change.
            "is_cofactor": TRUE_CELL,
            "source": ANIMALTFDB,
            "animaltfdb_family": animaltfdb["animaltfdb_family"],
            "animaltfdb_category": animaltfdb["animaltfdb_category"],
        }
    )
    check_table(table, origin=recipe.file_name)
    return table


def hgnc_stems(hgnc: pd.DataFrame, *, origin: str) -> dict[str, str]:
    """Return the archive's HGNC-id-to-Ensembl-gene-id map, which is what EpiFactors joins through.

    Raises
    ------
    CofactorSourceError
        If the archive names one HGNC id twice, which would make a gene's stem depend
        on which row was read last.
    """
    ids = hgnc["hgnc_id"].str.strip()
    repeated = sorted(ids[ids.duplicated()].unique())
    if repeated:
        raise CofactorSourceError(
            f"{origin} names the HGNC ids {repeated[:10]} more than once, so which Ensembl gene "
            f"id each of them carries depends on which row is read last. An archive keys on that "
            f"id — reconcile the file before building from it."
        )
    return dict(zip(ids, hgnc["ensembl_gene_id"].str.strip(), strict=True))


def hgnc_symbols(hgnc: pd.DataFrame, stems: Sequence[str], *, origin: str) -> dict[str, str]:
    """Return HGNC's current approved symbol for each of ``stems``, refusing anything ambiguous.

    The symbol on every human row, reached through an id and never through a name: an
    EpiFactors gene through its HGNC id, an AnimalTFDB gene through the stem that is
    already the table's key. That is what corrects the 31 rows EpiFactors still spells
    by a retired name.

    Raises
    ------
    CofactorSourceError
        If the archive names no gene at a stem, or names two — either way there is no
        one approved spelling to ship, and picking one would be this package inventing
        a name.
    """
    found: dict[str, set[str]] = {}
    for stem, symbol in zip(
        hgnc["ensembl_gene_id"].str.strip(), hgnc["symbol"].str.strip(), strict=True
    ):
        if stem:
            found.setdefault(stem, set()).add(symbol)
    # Ordered and deduplicated, so both messages below name genes in the table's order.
    wanted = list(dict.fromkeys(stems))
    unknown = [stem for stem in wanted if stem not in found]
    if unknown:
        raise CofactorSourceError(
            f"{origin} names no gene at the Ensembl gene ids {unknown[:10]} "
            f"({len(unknown)} in all), and every human row's symbol is HGNC's approved spelling. "
            f"Either the archive predates those genes or a publisher's stem is stale — pin the "
            f"archive that covers them, or reconcile the publisher's file."
        )
    ambiguous = [stem for stem in wanted if len(found[stem]) > 1]
    if ambiguous:
        raise CofactorSourceError(
            f"{origin} names two genes at each of the Ensembl gene ids {ambiguous[:10]}, so there "
            f"is no one approved symbol to ship for them and choosing one would be this package "
            f"naming a gene. Reconcile the archive before building from it."
        )
    return {stem: next(iter(found[stem])) for stem in wanted}


def published_values(cell: str, *, column: str, origin: str) -> list[str]:
    """Return one EpiFactors cell's values, split on the publisher's own comma-and-space.

    ``#`` is that publisher's *nothing recorded* and returns no values at all.

    Raises
    ------
    CofactorSourceError
        If a value already contains :data:`VALUE_SEPARATOR`, which would make the
        shipped cell split into values the publisher never wrote.
    """
    if cell.strip() == EPIFACTORS_NOTHING:
        return []
    values = [value.strip() for value in cell.split(EPIFACTORS_SEPARATOR)]
    kept = [value for value in values if value]
    carrying = [value for value in kept if VALUE_SEPARATOR in value]
    if carrying:
        raise CofactorSourceError(
            f"{origin} spells the {column!r} values {carrying} with a {VALUE_SEPARATOR!r} inside "
            f"them, and that is what a shipped cell separates two values with. Writing them "
            f"unchanged would make one value read as two. Re-spell them, or change "
            f"VALUE_SEPARATOR here and everywhere else in this package that a multi-valued cell "
            f"is written."
        )
    return kept


def epifactors_stems(
    epifactors: pd.DataFrame, stem_by_id: Mapping[str, str], *, origin: str, hgnc_origin: str
) -> list[str]:
    """Return the **Gene id stem** of every EpiFactors row, in the publisher's own row order.

    The join is on the HGNC id and never on the symbol. That is not a preference: 31 of
    the publisher's rows carry a name HGNC has since retired, so matching on the name
    would key those genes wrongly or drop them, where every row's id resolves.

    Raises
    ------
    CofactorSourceError
        If the archive names no row for an id, or gives one no Ensembl gene id. Either
        is a gene that would silently leave the table.
    """
    keys = [f"{HGNC_ID_PREFIX}{value.strip()}" for value in epifactors["HGNC_ID"]]
    unknown = sorted({key for key in keys if key not in stem_by_id})
    if unknown:
        raise CofactorSourceError(
            f"{hgnc_origin} names no gene under the HGNC ids {unknown[:10]} ({len(unknown)} in "
            f"all), which {origin} lists. Those genes have no Ensembl gene id and would leave the "
            f"table without a word — pin an archive that covers them, or reconcile the "
            f"publisher's file."
        )
    stemless = sorted({key for key in keys if not stem_by_id[key]})
    if stemless:
        raise CofactorSourceError(
            f"{hgnc_origin} gives no Ensembl gene id for the HGNC ids {stemless[:10]} "
            f"({len(stemless)} in all), which {origin} lists. A cofactor table is keyed by "
            f"**Gene id stem**, so a gene without one cannot ship — pin an archive that carries "
            f"them, or reconcile the publisher's file."
        )
    return [stem_by_id[key] for key in keys]


def collapse_epifactors(
    epifactors: pd.DataFrame, stems: Sequence[str], *, origin: str
) -> dict[str, dict[str, list[str]]]:
    """Return EpiFactors' four shipped columns per stem, its double-rowed genes unioned into one.

    Five genes carry two rows each and are not duplicates: one row annotates a histone
    modification and the other an RNA modification. A table is one row per stem, so the
    two are unioned, deduplicated within a cell as well as across the two rows, and
    first-seen order is kept so that two runs write the same bytes. The cost is that for
    those five the pairing between a function and its own modification is gone, which
    ``ATTRIBUTION.md`` says rather than leaving a reader to find out.

    Returns
    -------
    dict of str to dict of str to list of str
        Stem to shipped column name to that gene's values, in first-appearance order.
    """
    cells = {published: list(epifactors[published]) for published in EPIFACTORS_FIELDS}
    collapsed: dict[str, dict[str, list[str]]] = {}
    for index, stem in enumerate(stems):
        gene = collapsed.setdefault(stem, {name: [] for name in EPIFACTORS_COLUMNS})
        for published, name in EPIFACTORS_FIELDS.items():
            for value in published_values(cells[published][index], column=published, origin=origin):
                if value not in gene[name]:
                    gene[name].append(value)
    return collapsed


def row_source(stem: str, *, animaltfdb: Collection[str], epifactors: Collection[str]) -> str:
    """Return which publishers listed one gene, in the uniform ``source`` column's vocabulary."""
    if stem in animaltfdb and stem in epifactors:
        return BOTH
    return ANIMALTFDB if stem in animaltfdb else EPIFACTORS


def build_union_table(
    genes: pd.DataFrame,
    summary: pd.DataFrame,
    epifactors: pd.DataFrame,
    hgnc: pd.DataFrame,
    recipe: Recipe,
    *,
    origin: str,
    epifactors_origin: str,
    hgnc_origin: str,
) -> pd.DataFrame:
    """Return the shipped human table: two publishers' lists unioned, one row per **Gene id stem**.

    **Membership is unioned and classification is not.** A gene either publisher lists
    is a row; its AnimalTFDB columns are filled only if AnimalTFDB listed it and its
    EpiFactors columns only if EpiFactors did, so a blank group says that publisher
    never named the gene. Nothing crosses between the two vocabularies in either
    direction for any pair of values (ADR-0014), and ``source`` reading ``both`` is an
    agreement about membership and about nothing else.

    Row order is AnimalTFDB's own list and then the genes only EpiFactors names, in its
    order — each publisher's order kept where it has one, and deterministic, which is
    all byte-stability needs.
    """
    animaltfdb = animaltfdb_rows(genes, summary, recipe, origin=origin)
    listed = epifactors_stems(
        epifactors,
        hgnc_stems(hgnc, origin=hgnc_origin),
        origin=epifactors_origin,
        hgnc_origin=hgnc_origin,
    )
    epi = collapse_epifactors(epifactors, listed, origin=epifactors_origin)
    family = dict(zip(animaltfdb["gene_id_stem"], animaltfdb["animaltfdb_family"], strict=True))
    category = dict(zip(animaltfdb["gene_id_stem"], animaltfdb["animaltfdb_category"], strict=True))
    stems = [*animaltfdb["gene_id_stem"], *(stem for stem in epi if stem not in family)]
    symbol = hgnc_symbols(hgnc, stems, origin=hgnc_origin)
    table = pd.DataFrame(
        {
            "gene_id_stem": stems,
            # HGNC's approved spelling on every row, which is what corrects the 31 genes
            # EpiFactors still names by a symbol HGNC has retired.
            "symbol": [symbol[stem] for stem in stems],
            "is_cofactor": TRUE_CELL,
            "source": [row_source(stem, animaltfdb=family, epifactors=epi) for stem in stems],
            "animaltfdb_family": [family.get(stem, "") for stem in stems],
            "animaltfdb_category": [category.get(stem, "") for stem in stems],
            **{
                name: [VALUE_SEPARATOR.join(epi.get(stem, {}).get(name, ())) for stem in stems]
                for name in EPIFACTORS_COLUMNS
            },
        }
    )
    check_table(table, origin=recipe.file_name)
    return table


def check_table(table: pd.DataFrame, *, origin: str) -> None:
    """Hold the built table to what a shipped **Cofactor table** promises, or refuse to write it.

    Raises
    ------
    CofactorSourceError
        If the uniform four are not first, a **Gene id stem** is blank or repeated, the
        cofactor flag is spelled a way no table spells one, or a cell carries a tab or a
        newline — the shipped file is read by splitting on those, and quoting it would
        make a plain TSV that is not plainly readable.
    """
    columns = tuple(table.columns)
    if columns[: len(UNIFORM_COLUMNS)] != UNIFORM_COLUMNS:
        raise CofactorSourceError(
            f"{origin} would lead with {columns[: len(UNIFORM_COLUMNS)]} where every cofactor "
            f"table leads with {UNIFORM_COLUMNS}. Fix the column list in this script."
        )
    stems = table["gene_id_stem"]
    if (stems == "").any():
        raise CofactorSourceError(
            f"{origin} would carry a blank gene id stem. A cofactor table is keyed by stem, so "
            f"the publisher's file has to be reconciled first."
        )
    if stems.duplicated().any():
        repeated = sorted(stems[stems.duplicated()].unique())
        raise CofactorSourceError(
            f"{origin} would name these gene id stems more than once: {repeated}. One row per "
            f"gene, so the publisher's file has to be reconciled first."
        )
    flags = set(table["is_cofactor"].unique())
    if not flags <= {TRUE_CELL, FALSE_CELL}:
        raise CofactorSourceError(
            f"{origin} would spell its cofactor flag {sorted(flags)} where a table spells it "
            f"{TRUE_CELL!r} or {FALSE_CELL!r}. Fix the flag in this script."
        )
    for name in columns:
        if table[name].str.contains("[\t\r\n]", regex=True).any():
            raise CofactorSourceError(
                f"{origin} would carry a tab or a newline inside the {name!r} column. A cofactor "
                f"table is read by splitting on those, so such a cell has to be cleaned first."
            )


def render(table: pd.DataFrame) -> bytes:
    r"""Return the shipped table as unpacked TSV bytes — no quoting, ``\n`` throughout."""
    lines = ["\t".join(table.columns)]
    lines.extend("\t".join(row) for row in table.itertuples(index=False, name=None))
    return ("\n".join(lines) + "\n").encode("utf-8")


def write_table(unpacked: bytes, path: Path) -> None:
    """Write ``unpacked`` to ``path`` as gzip with no timestamp, so two runs agree byte for byte."""
    path.write_bytes(gzip.compress(unpacked, compresslevel=9, mtime=0))


def write_metadata(
    path: Path, columns: tuple[str, ...], rows: list[dict[str, str]], *, keys: tuple[str, ...]
) -> None:
    """Rewrite one provenance table with ``rows`` in place of every row they key over.

    ``keys`` is what makes a row the same row — the species for one table, the species
    and the source for the ragged one. The whole file is rewritten sorted by those keys,
    so a rebuild is a diff of the rows that changed and nothing else.
    """
    replaced = {tuple(row[key] for key in keys) for row in rows}
    kept = [
        dict(existing)
        for existing in read_metadata(path, columns)
        if tuple(existing[key] for key in keys) not in replaced
    ]
    kept.extend(rows)
    kept.sort(key=lambda row: tuple(row[key] for key in keys))
    lines = ["\t".join(columns)]
    lines.extend("\t".join(row[name] for name in columns) for row in kept)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_metadata(path: Path, columns: tuple[str, ...]) -> list[dict[str, str]]:
    """Return one provenance table's rows as text, empty when no table is there yet.

    Raises
    ------
    CofactorSourceError
        If the table on disk carries different columns from the ones being written,
        since merging into it would then produce a file neither shape can read.
    """
    if not path.is_file():
        return []
    lines = path.read_text(encoding="utf-8").splitlines()
    header = tuple(lines[0].split("\t"))
    check_header(header, columns, origin=str(path))
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:] if line]


def assemble(
    recipe: Recipe,
    genes_path: Path,
    summary_path: Path,
    *,
    epifactors_path: Path | None,
    hgnc_path: Path | None,
) -> pd.DataFrame:
    """Read this recipe's publishers' files and return the table they make.

    Raises
    ------
    CofactorSourceError
        If a recipe that unions EpiFactors was given neither its file nor the pinned
        HGNC archive, or a recipe that does not union it was given one anyway. Both are
        the wrong recipe name for the files at hand, and building the other table
        silently is what naming them here prevents.
    """
    genes = read_source(genes_path, GENE_LIST_COLUMNS)
    summary = read_source(summary_path, SUMMARY_COLUMNS)
    supplied = epifactors_path is not None and hgnc_path is not None
    if recipe.unions_epifactors != supplied:
        wanted = "needs --epifactors and --hgnc" if recipe.unions_epifactors else "takes neither"
        raise CofactorSourceError(
            f"the recipe for {recipe.species} {wanted}, and the command line said otherwise. "
            f"Only a table this package unions reads a second publisher's list, and only human "
            f"is one — check the recipe name against the files you have."
        )
    if epifactors_path is None or hgnc_path is None:
        return build_table(genes, summary, recipe, origin=str(genes_path))
    return build_union_table(
        genes,
        summary,
        read_source(epifactors_path, EPIFACTORS_GENE_COLUMNS, separator=","),
        read_named_columns(hgnc_path, HGNC_COLUMNS),
        recipe,
        origin=str(genes_path),
        epifactors_origin=str(epifactors_path),
        hgnc_origin=str(hgnc_path),
    )


def report(
    table: pd.DataFrame, recipe: Recipe, *, unpacked: bytes, digest: str, packed: int
) -> None:
    """Print what was written, in the numbers a reviewer checks a rebuild against.

    A blank cell counts towards nothing: it is a publisher who never named the gene,
    and counting it as a value would report a family and a category too many for every
    table with a second publisher in it.
    """
    families = len(set(table["animaltfdb_family"]) - {""})
    categories = len(set(table["animaltfdb_category"]) - {""})
    print(f"{recipe.file_name}: {len(table)} rows, {families} families, {categories} categories")
    print(f"  columns: {', '.join(table.columns)}")
    if recipe.unions_epifactors:
        counts = table["source"].value_counts()
        split = ", ".join(
            f"{counts.get(name, 0)} {name}" for name in (BOTH, ANIMALTFDB, EPIFACTORS)
        )
        print(f"  source:  {split}")
        for name in EPIFACTORS_COLUMNS:
            values = {
                value for cell in table[name] for value in cell.split(VALUE_SEPARATOR) if value
            }
            print(f"  {name}: {len(values)} distinct values over {(table[name] != '').sum()} genes")
    print(f"  sources: {', '.join(source.source for source in recipe.sources)}")
    print(f"  unpacked {len(unpacked)} bytes, sha256 {digest}")
    print(f"  packed   {packed} bytes ({len(unpacked) / packed:.1f}x)")


def build(
    table_name: str,
    genes_path: Path,
    summary_path: Path,
    data_dir: Path,
    *,
    epifactors_path: Path | None = None,
    hgnc_path: Path | None = None,
) -> None:
    """Build one cofactor table and its provenance rows, and report what was written."""
    recipe = RECIPES[table_name]
    table = assemble(
        recipe,
        genes_path,
        summary_path,
        epifactors_path=epifactors_path,
        hgnc_path=hgnc_path,
    )
    unpacked = render(table)
    data_dir.mkdir(parents=True, exist_ok=True)
    write_table(unpacked, data_dir / recipe.file_name)
    digest = hashlib.sha256(unpacked).hexdigest()
    write_metadata(
        data_dir / METADATA_FILE,
        METADATA_COLUMNS,
        [
            {
                "species": recipe.species,
                "ncbi_taxid": str(recipe.ncbi_taxid),
                "file": recipe.file_name,
                "sha256": digest,
            }
        ],
        keys=("species",),
    )
    write_metadata(
        data_dir / SOURCE_METADATA_FILE,
        SOURCE_METADATA_COLUMNS,
        [
            {
                "species": recipe.species,
                "source": source.source,
                "publisher": source.publisher,
                "version": source.version,
                "pubmed_id": str(source.pubmed_id),
                "source_url": source.source_url,
            }
            for source in recipe.sources
        ],
        keys=("species", "source"),
    )
    report(
        table,
        recipe,
        unpacked=unpacked,
        digest=digest,
        packed=(data_dir / recipe.file_name).stat().st_size,
    )


def main(argv: Sequence[str] | None = None) -> None:
    """Parse the command line and build the cofactor table it names."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("table", choices=sorted(RECIPES), help="which cofactor table to build")
    parser.add_argument("genes", type=Path, help="the publisher's per-species cofactor list")
    parser.add_argument("summary", type=Path, help="the publisher's cof_info_summary.tsv")
    parser.add_argument(
        "--epifactors",
        type=Path,
        default=None,
        help="EpiFactors' EpiGenes_main.csv, for the table that unions it in",
    )
    parser.add_argument(
        "--hgnc",
        type=Path,
        default=None,
        help="the pinned dated HGNC monthly archive, which supplies EpiFactors' gene id stems",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / DATA_SUBDIR,
        help="where the table and its two provenance tables are written",
    )
    arguments = parser.parse_args(argv)
    try:
        build(
            arguments.table,
            arguments.genes,
            arguments.summary,
            arguments.data_dir,
            epifactors_path=arguments.epifactors,
            hgnc_path=arguments.hgnc,
        )
    except CofactorSourceError as error:
        # The message already says what moved and what to do about it; a traceback on
        # top of it would only bury that.
        raise SystemExit(str(error)) from None


if __name__ == "__main__":
    main()
