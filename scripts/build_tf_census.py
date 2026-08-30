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

``animaltfdb4_mouse`` — Shen *et al.*, "AnimalTFDB 4.0", *Nucleic Acids Research*
51(D1):D39-D45, 2023 (PMID 36268869). Mouse TF list, downloaded from
https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/Mus_musculus_TF.

Cite the publisher whose census you use.

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

*It writes rows in the publisher's own order*, which is already deterministic and
is all byte-stability needs of a recipe. Everything after that — the rendering, the
compression and the provenance merge — is :mod:`genome.shipped_writer`'s.

Usage
-----
``python scripts/build_tf_census.py <census> <the publisher's file>``, where
``<census>`` names one of :data:`RECIPES`.

Add ``--data-dir`` to write somewhere other than ``src/genome/data/tf_gene``.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from dataclasses import dataclass, replace
from pathlib import Path

import pandas as pd

# The format is the package's and this script writes through it: the columns, the file
# name and every rule a census is held to are declared once, beside the reader, and the
# writer holds what is about to be written to exactly them.
from genome.shipped_writer import merge_rows, shipped_name, write_table
from genome.tf.gene import (
    CENSUS_FORMAT,
    CENSUS_METADATA_FORMAT,
    CENSUS_SUBDIR,
    TRUE_CELL,
    species_slug,
)

#: Where the shipped censuses live, relative to the repository root — the package's own
#: resource directory, under the source tree it is packaged from.
DATA_SUBDIR = Path("src") / "genome" / CENSUS_SUBDIR

#: The publishers' own spelling of *nothing recorded here*, which ships as a blank
#: cell instead — blank is how every other table in this package spells unknown.
#: Deliberately one literal and not a list of plausible ones: blanking a spelling no
#: publisher uses would silently erase a real value that happens to look like it.
_ABSENT_CELLS = frozenset({"None"})


class CensusSourceError(ValueError):
    """The census asked for cannot be built.

    Either the publisher's file is not the file the recipe was written against, the
    recipe is not the shape a recipe has, or what the recipe would write is not a file
    a census reader would accept back.
    """


#: What repairs a refused build: never *re-run this*, which is what the reader's own
#: refusals say, because re-running is what has just failed.
_REPAIR = (
    "fix the recipe in scripts/build_tf_census.py, or reconcile the publisher's file "
    "before building from it"
)

#: A census as this script writes one, and the provenance table beside it: the package's
#: own declarations of both formats, refusing in this script's words. Everything about
#: what a census *is* comes from :mod:`genome.tf.gene.census`; only the error class and
#: the repair are the build's.
WRITTEN_CENSUS = replace(CENSUS_FORMAT, error=CensusSourceError, repair=_REPAIR)
WRITTEN_METADATA = replace(
    CENSUS_METADATA_FORMAT,
    error=CensusSourceError,
    repair="fix that file's header, or move it aside and let this script write it afresh",
)


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
    shipped : tuple of tuple of (str or None) and str
        ``(published name, shipped name)`` for the columns that ship, in shipped
        order. The first four must be the uniform columns every census leads with,
        which ``WRITTEN_CENSUS.columns`` declares. A published name of
        ``None`` means this script supplies the column because the publisher has no
        such column, which is legal for the TF flag and nothing else: a publisher
        that lists none but the genes it accepts publishes no rejected set, so every
        row it does publish is a ``yes``, and writing ``no`` rows for the genes it
        left out would fabricate a verdict nobody reached (ADR-0014).
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
    shipped: tuple[tuple[str | None, str], ...]

    def __post_init__(self) -> None:
        """Refuse a recipe that supplies any column but the TF flag out of thin air.

        Raises
        ------
        CensusSourceError
            If ``shipped`` sources anything else from ``None``. Every other column
            carries the publisher's own content, and a supplied one would carry this
            script's — a whole column of invented cells that nothing downstream could
            tell from published ones.
        """
        supplied = [name for published, name in self.shipped if published is None]
        if supplied not in ([], ["is_tf"]):
            raise CensusSourceError(
                f"the {self.species} recipe supplies the columns {supplied} rather than reading "
                f"them from the publisher's file. Only the TF flag may be supplied, and only for "
                f"a publisher that lists none but the genes it accepts. Name each of those "
                f"columns' published spelling in this script's 'shipped' list."
            )

    @property
    def slug(self) -> str:
        """Return the species slug this census's file is named by."""
        return species_slug(self.species)

    @property
    def file_name(self) -> str:
        """Return the census's file name, as the package's own declaration spells it."""
        return shipped_name(CENSUS_FORMAT, slug=self.slug)


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

ANIMALTFDB4_MOUSE = Recipe(
    species="Mus musculus",
    ncbi_taxid=10090,
    publisher="AnimalTFDB",
    version="4.0",
    pubmed_id=36268869,
    source_url=(
        "https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/Mus_musculus_TF"
    ),
    family_column="Family",
    separator="\t",
    # Six columns, none of them carrying trailing whitespace and none of them a TF
    # flag: AnimalTFDB lists the genes it judges transcription factors and no others.
    published=("Species", "Symbol", "Ensembl", "Family", "Protein", "Entrez_ID"),
    # Four ship, which is every judgement AnimalTFDB makes about a gene — that it is a
    # transcription factor, and of which family — and nothing else. Dropped: ``Species``,
    # the same word on all 1,611 rows and already the provenance row's; and the two
    # identifier cross-references, ``Protein`` (a semicolon-joined blob of Ensembl 105
    # protein ids, the class Lambert's ``Pfam Domains (By ENSP ID)`` was dropped for) and
    # ``Entrez_ID`` (NCBI's own id, already ``NA`` for 109 of the genes). The asymmetry
    # against Lambert's thirteen is the two publishers' and not this script's.
    shipped=(
        ("Ensembl", "gene_id_stem"),
        ("Symbol", "symbol"),
        # Supplied, because AnimalTFDB publishes no flag and no rejected set. Every gene
        # in the file is one it accepts, and the genes it left out are ones it says
        # nothing about — so mouse has no assessed-negative rows rather than fabricated
        # ones, and that absence is the census's shape rather than a defect (ADR-0014).
        (None, "is_tf"),
        ("Family", "dbd_family"),
    ),
)

#: Every census this script knows how to build, by the name given on the command line.
RECIPES: dict[str, Recipe] = {
    "lambert2018": LAMBERT_2018,
    "animaltfdb4_mouse": ANIMALTFDB4_MOUSE,
}


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

    A column the recipe sources from ``None`` is supplied as ``yes`` on every row,
    which the recipe holds to the TF flag alone. What a census must be — the uniform
    columns first, one row per **Gene id stem**, a flag spelled the one way — is
    checked as it is written, against the package's own declaration of the format.
    """
    shipped = pd.DataFrame(
        {
            name: (
                frame[published].map(clean_cell)
                if published is not None
                else pd.Series(TRUE_CELL, index=frame.index, dtype=object)
            )
            for published, name in recipe.shipped
        }
    )
    shipped["is_tf"] = shipped["is_tf"].str.lower()
    return shipped


def clean_cell(value: str) -> str:
    """Return one cell's text with surrounding whitespace and the publisher's *absent* spelling gone."""
    stripped = value.strip()
    return "" if stripped in _ABSENT_CELLS else stripped


def build(census: str, source: Path, data_dir: Path) -> None:
    """Build one census and its provenance row, and report what was written."""
    recipe = RECIPES[census]
    table = build_table(read_source(source, recipe), recipe)
    written = write_table(
        WRITTEN_CENSUS,
        data_dir,
        tuple(table.columns),
        table.itertuples(index=False, name=None),
        slug=recipe.slug,
    )
    merge_rows(
        WRITTEN_METADATA,
        data_dir,
        [
            {
                "species": recipe.species,
                "ncbi_taxid": str(recipe.ncbi_taxid),
                "file": recipe.file_name,
                "publisher": recipe.publisher,
                "version": recipe.version,
                "pubmed_id": str(recipe.pubmed_id),
                "family_column": recipe.family_column,
                "source_url": recipe.source_url,
                "sha256": written.sha256,
            }
        ],
    )
    positive = int((table["is_tf"] == TRUE_CELL).sum())
    print(f"{recipe.file_name}: {len(table)} rows, {positive} assessed positive")
    print(f"  columns: {', '.join(table.columns)}")
    print(f"  unpacked {len(written.unpacked)} bytes, sha256 {written.sha256}")
    print(f"  packed   {written.packed} bytes ({len(written.unpacked) / written.packed:.1f}x)")


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
