#!/usr/bin/env python
r"""Rebuild one shipped **Cofactor table** from AnimalTFDB's own two files.

Run by hand, never at install time and never in CI — the publisher's files are
downloads CI cannot make, which is why the shipped tables are guarded by
``tests/test_tf_cofactor.py`` instead. This script lives outside the package for
the same reason: the wheel carries the table and not the tooling that made it.

Attribution
-----------
``animaltfdb4_mouse`` and ``animaltfdb4_worm`` — Shen *et al.*, "AnimalTFDB 4.0",
*Nucleic Acids Research* 51(D1):D39-D45, 2023 (PMID 36268869). Two files per
species, both from https://guolab.wchscu.cn/: the per-species cofactor list
(``.../Cof_list_final/<Genus_species>_Cof``), which is membership and family, and
the family summary (``.../cof_info_summary.tsv``), which is the category each
family sits in. The all-species bulk file is deliberately not read: the per-species
download is the one whose row count the summary reconciles against.

Cite the publisher whose table you use.

What it does
------------
Reads the publisher's two files, holds both headers to the exact columns that
release published, joins each gene's category onto its family, and writes a
gzipped TSV plus one row in each of the two provenance tables beside it. Three
properties are the point:

*It fails loudly on an unrecognised header.* A publisher who re-spells ``Family``
would otherwise have that column silently dropped, and nothing downstream would
say so. Every published name is listed here, so a re-spelling is one missing name
and one unexpected name, and both are named in the error.

*It fails loudly on a family with no category.* The publisher's own two files
spell some families differently, so the join goes through the hand-written
:data:`FAMILY_SPELLINGS` below. A family that survives that map and is still not in
the summary breaks the build rather than blanking a column, which is what a future
release renaming a family must do.

*It writes byte-stable output.* Rows keep the publisher's own order, cells are
written with no quoting, the line terminator is ``\n``, and gzip is given
``mtime=0``, so re-running on unchanged inputs produces no diff.

Usage
-----
``python scripts/build_tf_cofactor.py <table> <the species cofactor list>
<cof_info_summary.tsv>``, where ``<table>`` names one of :data:`RECIPES`.

Add ``--data-dir`` to write somewhere other than ``src/genome/data/tf_cofactor``.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
from collections.abc import Mapping, Sequence
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

#: The value of the uniform ``source`` column on every row this script writes. One of
#: the closed vocabulary ``genome.tf.cofactor`` validates on read.
ANIMALTFDB = "animaltfdb"

#: The header the per-species cofactor list publishes, in order.
GENE_LIST_COLUMNS = ("Species", "Symbol", "Ensembl", "Family", "Entrez_ID")

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
class Recipe:
    """How one species' AnimalTFDB download becomes a shipped **Cofactor table**.

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
        Who published the table.
    version : str
        The publisher's own release identifier.
    pubmed_id : int
        The PubMed id of the paper to cite.
    source_url : str
        Where the per-species cofactor list is downloaded from. The one URL the
        provenance row carries, because it is the one that names this species; the
        summary file is the same release and is named in ``ATTRIBUTION.md``.
    summary_url : str
        Where the family summary is downloaded from.
    """

    species: str
    ncbi_taxid: int
    animaltfdb_species: str
    publisher: str
    version: str
    pubmed_id: int
    source_url: str
    summary_url: str

    @property
    def slug(self) -> str:
        """Return the species slug this table's file is named by."""
        return species_slug(self.species)

    @property
    def file_name(self) -> str:
        """Return the table's file name."""
        return f"{self.slug}{COFACTOR_SUFFIX}"


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

#: Every table this script knows how to build, by the name given on the command line.
RECIPES: dict[str, Recipe] = {
    "animaltfdb4_mouse": ANIMALTFDB4_MOUSE,
    "animaltfdb4_worm": ANIMALTFDB4_WORM,
}


def read_source(path: Path, columns: tuple[str, ...]) -> pd.DataFrame:
    """Read one of the publisher's files, refusing it unless its header is the one expected.

    Every cell is text: the publisher's own missing-value spellings are not read as
    NaN here, since which of them a column uses is a fact about that column worth
    seeing rather than one pandas should guess at.
    """
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    frame.columns = pd.Index([str(name).strip() for name in frame.columns])
    check_header(tuple(frame.columns), columns, origin=str(path))
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


def build_table(
    genes: pd.DataFrame, summary: pd.DataFrame, recipe: Recipe, *, origin: str
) -> pd.DataFrame:
    """Return the shipped table: the uniform four then AnimalTFDB's two, in the publisher's own row order.

    The publisher's order is kept rather than sorted by **Gene id stem**. It is already
    deterministic, which is all byte-stability needs, and the file is ordered by family
    — sorting it scatters exactly the runs that compress.

    The family ships as the *gene list's* own spelling, never the summary's: the
    spelling map exists to find each gene's category and not to re-spell what the
    publisher classified it under.
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
    table = pd.DataFrame(
        {
            "gene_id_stem": rows["Ensembl"].str.strip(),
            "symbol": rows["Symbol"].str.strip(),
            # Every gene AnimalTFDB lists is one it accepts: it publishes no rejected
            # set, so the flag is supplied and reads `yes` on every row. It is kept
            # anyway — dropping it would make presence in the file the verdict, and a
            # future source could then not record a rejection without a format change.
            "is_cofactor": TRUE_CELL,
            "source": ANIMALTFDB,
            "animaltfdb_family": family,
            "animaltfdb_category": [category[map_family(name)] for name in family],
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


def build(table_name: str, genes_path: Path, summary_path: Path, data_dir: Path) -> None:
    """Build one cofactor table and its two provenance rows, and report what was written."""
    recipe = RECIPES[table_name]
    genes = read_source(genes_path, GENE_LIST_COLUMNS)
    summary = read_source(summary_path, SUMMARY_COLUMNS)
    table = build_table(genes, summary, recipe, origin=str(genes_path))
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
                "source": ANIMALTFDB,
                "publisher": recipe.publisher,
                "version": recipe.version,
                "pubmed_id": str(recipe.pubmed_id),
                "source_url": recipe.source_url,
            }
        ],
        keys=("species", "source"),
    )
    packed = (data_dir / recipe.file_name).stat().st_size
    families = table["animaltfdb_family"].nunique()
    categories = table["animaltfdb_category"].nunique()
    print(f"{recipe.file_name}: {len(table)} rows, {families} families, {categories} categories")
    print(f"  columns: {', '.join(table.columns)}")
    print(f"  unpacked {len(unpacked)} bytes, sha256 {digest}")
    print(f"  packed   {packed} bytes ({len(unpacked) / packed:.1f}x)")


def main(argv: Sequence[str] | None = None) -> None:
    """Parse the command line and build the cofactor table it names."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("table", choices=sorted(RECIPES), help="which cofactor table to build")
    parser.add_argument("genes", type=Path, help="the publisher's per-species cofactor list")
    parser.add_argument("summary", type=Path, help="the publisher's cof_info_summary.tsv")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / DATA_SUBDIR,
        help="where the table and its two provenance tables are written",
    )
    arguments = parser.parse_args(argv)
    try:
        build(arguments.table, arguments.genes, arguments.summary, arguments.data_dir)
    except CofactorSourceError as error:
        # The message already says what moved and what to do about it; a traceback on
        # top of it would only bury that.
        raise SystemExit(str(error)) from None


if __name__ == "__main__":
    main()
