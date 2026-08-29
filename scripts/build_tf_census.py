#!/usr/bin/env python
r"""Rebuild one shipped **TF gene table** from its publisher's own file.

Run by hand, never at install time and never in CI — the publishers' files are
downloads CI cannot make, which is why the shipped tables are guarded by
``tests/test_tf_census.py`` instead. This script lives outside the package for
the same reason: the wheel carries the census and not the tooling that made it.

Attribution
-----------
``lambert2018`` — Lambert *et al.*, "The Human Transcription Factors", *Cell*
172(4):650-665, 2018 (PMID 29425488). Database extract v_1.01, downloaded from
https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv.
Cite the paper when you use the census.

What it does
------------
Reads the publisher's own file, holds its header to the exact set of columns that
release published, keeps the ones the package ships, and writes a gzipped TSV plus
a row in the provenance table beside it. Two properties are the point:

*It fails loudly on an unrecognised header.* A publisher who re-spells ``DBD``
would otherwise have that column silently dropped from a re-generated census, and
nothing downstream would say so. Every published name is listed here, so a
re-spelling is one missing name and one unexpected name, and both are named in the
error.

*It writes byte-stable output.* Rows keep the publisher's own order, cells are
written with no quoting, the line terminator is ``\n``, and gzip is given
``mtime=0``, so re-running on an unchanged input produces no diff.

Usage
-----
``python scripts/build_tf_census.py lambert2018 <the publisher's file>``

Add ``--data-dir`` to write somewhere other than ``src/genome/data/tf_gene``.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

#: Where the shipped censuses live, relative to the repository root.
DATA_SUBDIR = Path("src") / "genome" / "data" / "tf_gene"

#: What one census file is called: the species slug, then this. Kept in step with
#: ``genome.tf.gene.census.CENSUS_SUFFIX``, which is what enumerates them.
CENSUS_SUFFIX = ".tf_gene_table.tsv.gz"

#: The provenance table beside the censuses.
METADATA_FILE = "census_metadata.tsv"

#: Its columns, in table order.
METADATA_COLUMNS = (
    "species",
    "ncbi_taxid",
    "file",
    "publisher",
    "version",
    "pubmed_id",
    "family_column",
    "source_url",
    "sha256",
)

#: The four columns every census carries under the same name, in the same place.
UNIFORM_COLUMNS = ("gene_id_stem", "symbol", "is_tf", "dbd_family")

#: How the uniform TF flag is spelled in a shipped census.
TRUE_CELL, FALSE_CELL = "yes", "no"

#: The publishers' own spelling of *nothing recorded here*, which ships as a blank
#: cell instead — blank is how every other table in this package spells unknown.
#: Deliberately one literal and not a list of plausible ones: blanking a spelling no
#: publisher uses would silently erase a real value that happens to look like it.
_ABSENT_CELLS = frozenset({"None"})


class CensusSourceError(ValueError):
    """The publisher's file is not the file this recipe was written against."""


@dataclass(frozen=True)
class Recipe:
    """How one publisher's census is turned into a shipped **TF gene table**.

    Attributes
    ----------
    species : str
        The species as the assembly metadata table spells it, e.g. ``"Homo sapiens"``.
        Its slug is the census's file name.
    ncbi_taxid : int
        NCBI taxonomy id for that species.
    publisher : str
        Who published the census.
    version : str
        The publisher's own release identifier.
    pubmed_id : int
        The PubMed id of the paper to cite.
    source_url : str
        Where the publisher's file is downloaded from.
    family_column : str
        The publisher's own name for the column the **DBD family** is taken from —
        recorded because the two vocabularies are deliberately not crosswalked
        (ADR-0014), so a reader has to know whose family names they are grouping by.
    separator : str
        The publisher's field separator.
    published : tuple of str
        Every column that release publishes, in order and with trailing whitespace
        stripped. The header is held to exactly this.
    shipped : tuple of tuple of str
        ``(published name, shipped name)`` for the columns that ship, in shipped
        order. The first four must be :data:`UNIFORM_COLUMNS`.
    """

    species: str
    ncbi_taxid: int
    publisher: str
    version: str
    pubmed_id: int
    source_url: str
    family_column: str
    separator: str
    published: tuple[str, ...]
    shipped: tuple[tuple[str, str], ...]

    @property
    def slug(self) -> str:
        """Return the species slug this census's file is named by."""
        return species_slug(self.species)

    @property
    def file_name(self) -> str:
        """Return the census's file name."""
        return f"{self.slug}{CENSUS_SUFFIX}"


def species_slug(species: str) -> str:
    """Return the file-name spelling of ``species``: lower case, underscores for runs of anything else."""
    kept = [character if character.isalnum() else " " for character in species.strip().lower()]
    return "_".join("".join(kept).split())


LAMBERT_2018 = Recipe(
    species="Homo sapiens",
    ncbi_taxid=9606,
    publisher="Lambert et al. 2018",
    version="v_1.01",
    pubmed_id=29425488,
    source_url="https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv",
    family_column="DBD",
    separator=",",
    # Twenty-nine columns. The first is the publisher's unnamed row-index column,
    # which pandas names ``Unnamed: 0``; three of the rest carry trailing whitespace
    # in the published file and are normalised out of it here rather than at read time.
    published=(
        "Unnamed: 0",
        "Ensembl ID",
        "HGNC symbol",
        "DBD",
        "Is TF?",
        "TF assessment",
        "Binding mode",
        "Motif status",
        "Final Notes",
        "Final Comments",
        "Interpro ID(s)",
        "EntrezGene ID",
        "EntrezGene Description",
        "PDB ID",
        "TF tested by HT-SELEX?",
        "TF tested by PBM?",
        "Conditional Binding Requirements",
        "Original Comments",
        "Vaquerizas 2009 classification",
        "CisBP considers it a TF?",
        "TFCat classification",
        "Is a GO TF?",
        "Initial assessment",
        "Curator 1",
        "Curator 2",
        "TFclass considers it a TF?",
        "Go Evidence",
        "Pfam Domains (By ENSP ID)",
        "Is C2H2 ZF(KRAB)?",
    ),
    # Thirteen of them ship. Dropped: the row index; the identifier blobs keyed to
    # 2018 Ensembl protein and transcript ids (``Pfam Domains (By ENSP ID)``,
    # ``Go Evidence``); the free-text curation notes (``Final Notes``,
    # ``Final Comments``, ``Original Comments``, ``Conditional Binding Requirements``,
    # ``Initial assessment``); the curator names; and the columns unusable without
    # them — the EntrezGene id and description, the PDB id, the two "tested by"
    # columns and the TFclass vote, each of which is a 2018 cross-reference into a
    # database that has moved on.
    shipped=(
        ("Ensembl ID", "gene_id_stem"),
        ("HGNC symbol", "symbol"),
        ("Is TF?", "is_tf"),
        ("DBD", "dbd_family"),
        ("TF assessment", "tf_assessment"),
        ("Binding mode", "binding_mode"),
        ("Motif status", "motif_status"),
        ("Interpro ID(s)", "interpro_ids"),
        ("Vaquerizas 2009 classification", "vaquerizas_2009_classification"),
        ("CisBP considers it a TF?", "cisbp_considers_it_a_tf"),
        ("TFCat classification", "tfcat_classification"),
        ("Is a GO TF?", "is_a_go_tf"),
        ("Is C2H2 ZF(KRAB)?", "is_c2h2_zf_krab"),
    ),
)

#: Every census this script knows how to build, by the name given on the command line.
RECIPES: dict[str, Recipe] = {"lambert2018": LAMBERT_2018}


def read_source(path: Path, recipe: Recipe) -> pd.DataFrame:
    """Read the publisher's file, refusing it unless its header is the one expected.

    Every cell is text: the publisher's own missing-value spellings are not read as
    NaN here, since which of them a column uses is a fact about that column worth
    seeing rather than one pandas should guess at.
    """
    frame = pd.read_csv(path, sep=recipe.separator, dtype=str, keep_default_na=False)
    frame.columns = pd.Index([str(name).strip() for name in frame.columns])
    check_header(tuple(frame.columns), recipe, origin=str(path))
    return frame


def check_header(seen: tuple[str, ...], recipe: Recipe, *, origin: str) -> None:
    """Hold ``seen`` to the columns this recipe was written against, naming every difference.

    Raises
    ------
    CensusSourceError
        If a published column is missing or an unexpected one appears. A publisher
        who re-spells a column shows up as both, and both are named — dropping the
        column silently is the outcome this exists to prevent.
    """
    missing = [name for name in recipe.published if name not in seen]
    unexpected = [name for name in seen if name not in recipe.published]
    if not missing and not unexpected:
        return
    parts = []
    if missing:
        parts.append(f"columns this recipe expects and {origin} does not carry: {missing}")
    if unexpected:
        parts.append(f"columns {origin} carries and this recipe does not name: {unexpected}")
    raise CensusSourceError(
        f"{origin} is not the file this recipe was written against — "
        + "; ".join(parts)
        + ". A publisher who re-spells a column would otherwise have it silently dropped, "
        "so update the recipe's 'published' and 'shipped' lists in this script to match "
        "the release you are building from, and say in the commit which columns moved."
    )


def build_table(frame: pd.DataFrame, recipe: Recipe) -> pd.DataFrame:
    """Return the shipped table: the kept columns, renamed and cleaned, in the publisher's own row order.

    The publisher's order is kept rather than sorted by **Gene id stem**. It is
    already deterministic, which is all byte-stability needs, and Lambert's file is
    ordered by domain family — sorting it costs 8 KB of gzip because it scatters
    exactly the runs that compress.
    """
    shipped = pd.DataFrame(
        {name: frame[published].map(clean_cell) for published, name in recipe.shipped}
    )
    shipped["is_tf"] = shipped["is_tf"].str.lower()
    check_table(shipped, origin=recipe.file_name)
    return shipped


def clean_cell(value: str) -> str:
    """Return one cell's text with surrounding whitespace and the publisher's *absent* spelling gone."""
    stripped = value.strip()
    return "" if stripped in _ABSENT_CELLS else stripped


def check_table(table: pd.DataFrame, *, origin: str) -> None:
    """Hold the built table to what a shipped census promises, or refuse to write it.

    Raises
    ------
    CensusSourceError
        If the uniform four are not first, a **Gene id stem** is blank or repeated,
        the TF flag is spelled a way no census spells one, or a cell carries a tab or
        a newline — the shipped file is read by splitting on those, and quoting it
        would make a plain TSV that is not plainly readable.
    """
    columns = tuple(table.columns)
    if columns[: len(UNIFORM_COLUMNS)] != UNIFORM_COLUMNS:
        raise CensusSourceError(
            f"{origin} would lead with {columns[: len(UNIFORM_COLUMNS)]} where every census "
            f"leads with {UNIFORM_COLUMNS}. Fix the recipe's 'shipped' list."
        )
    stems = table["gene_id_stem"]
    if (stems == "").any():
        raise CensusSourceError(f"{origin} would carry a blank gene id stem — fix the recipe.")
    if stems.duplicated().any():
        repeated = sorted(stems[stems.duplicated()].unique())
        raise CensusSourceError(
            f"{origin} would name these gene id stems more than once: {repeated}. A census "
            f"answers one verdict per gene, so the publisher's file has to be reconciled first."
        )
    flags = set(table["is_tf"].unique())
    if not flags <= {TRUE_CELL, FALSE_CELL}:
        raise CensusSourceError(
            f"{origin} would spell its TF flag {sorted(flags)} where a census spells it "
            f"{TRUE_CELL!r} or {FALSE_CELL!r}. Map the publisher's spelling in the recipe."
        )
    for name in columns:
        if table[name].str.contains("[\t\r\n]", regex=True).any():
            raise CensusSourceError(
                f"{origin} would carry a tab or a newline inside the {name!r} column. A census "
                f"is read by splitting on those, so such a cell has to be cleaned first."
            )


def render(table: pd.DataFrame) -> bytes:
    r"""Return the shipped table as unpacked TSV bytes — no quoting, ``\n`` throughout."""
    lines = ["\t".join(table.columns)]
    lines.extend("\t".join(row) for row in table.itertuples(index=False, name=None))
    return ("\n".join(lines) + "\n").encode("utf-8")


def write_census(unpacked: bytes, path: Path) -> None:
    """Write ``unpacked`` to ``path`` as gzip with no timestamp, so two runs agree byte for byte."""
    path.write_bytes(gzip.compress(unpacked, compresslevel=9, mtime=0))


def write_metadata(path: Path, row: dict[str, str]) -> None:
    """Rewrite the provenance table with ``row`` in place of any row for the same species."""
    rows = [
        dict(existing) for existing in read_metadata(path) if existing["species"] != row["species"]
    ]
    rows.append(row)
    rows.sort(key=lambda existing: existing["species"])
    lines = ["\t".join(METADATA_COLUMNS)]
    lines.extend("\t".join(entry[name] for name in METADATA_COLUMNS) for entry in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_metadata(path: Path) -> list[dict[str, str]]:
    """Return the provenance table's rows as text, empty when no table is there yet."""
    if not path.is_file():
        return []
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:] if line]


def build(census: str, source: Path, data_dir: Path) -> None:
    """Build one census and its provenance row, and report what was written."""
    recipe = RECIPES[census]
    table = build_table(read_source(source, recipe), recipe)
    unpacked = render(table)
    data_dir.mkdir(parents=True, exist_ok=True)
    write_census(unpacked, data_dir / recipe.file_name)
    digest = hashlib.sha256(unpacked).hexdigest()
    write_metadata(
        data_dir / METADATA_FILE,
        {
            "species": recipe.species,
            "ncbi_taxid": str(recipe.ncbi_taxid),
            "file": recipe.file_name,
            "publisher": recipe.publisher,
            "version": recipe.version,
            "pubmed_id": str(recipe.pubmed_id),
            "family_column": recipe.family_column,
            "source_url": recipe.source_url,
            "sha256": digest,
        },
    )
    packed = (data_dir / recipe.file_name).stat().st_size
    positive = int((table["is_tf"] == TRUE_CELL).sum())
    print(f"{recipe.file_name}: {len(table)} rows, {positive} assessed positive")
    print(f"  columns: {', '.join(table.columns)}")
    print(f"  unpacked {len(unpacked)} bytes, sha256 {digest}")
    print(f"  packed   {packed} bytes ({len(unpacked) / packed:.1f}x)")


def main(argv: Sequence[str] | None = None) -> None:
    """Parse the command line and build the census it names."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("census", choices=sorted(RECIPES), help="which census to build")
    parser.add_argument("source", type=Path, help="the publisher's own file")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / DATA_SUBDIR,
        help="where the census and its provenance table are written",
    )
    arguments = parser.parse_args(argv)
    try:
        build(arguments.census, arguments.source, arguments.data_dir)
    except CensusSourceError as error:
        # The message already says which columns moved and what to do about it; a
        # traceback on top of it would only bury that.
        raise SystemExit(str(error)) from None


if __name__ == "__main__":
    main()
