"""Tests for the Motif link tables shipped inside the package.

The shipped files answer here; no fixture stands in for them. That is the point of this
module rather than an accident of it: building a link table needs a JASPAR download CI
cannot make, so a guard over what ships is the only check that can exist. **No test here
proves a shipped table still matches JASPAR, and none can** — the counts below are pinned
so that drift becomes a loud failure instead of a table that quietly stopped agreeing
with the release it names, which is the most that is available (ADR-0015).

The tables are read as the files they are — gzipped TSV, no quoting — with :mod:`gzip`
and :mod:`csv` rather than through any module of this package, because reading them
without importing it is the property being defended.

Every test that is about one table iterates over *every* shipped table. A second species
or a third release is a file dropped into ``data/tf_link/``, and dropping one in must not
mean rewriting these.
"""

from __future__ import annotations

import csv
import gzip
import re
from importlib.resources import files
from importlib.resources.abc import Traversable

import pytest

from genome.tf.gene import TFGeneTable, census_species, species_slug, tf_gene_table
from genome.tf.motif.jaspar import JASPAR_RELEASES

#: Where the shipped link tables live inside the package. Kept in step with
#: ``scripts/build_tf_links.py``, which writes them; the module that reads them is not
#: built yet, and until it is these tests are the only reader.
LINK_SUBDIR = "data/tf_link"

#: What one link table is called: the species slug, the release under this prefix, then
#: the suffix — ``homo_sapiens.jaspar2026.motif_link_table.tsv.gz``.
RELEASE_PREFIX = "jaspar"
LINK_SUFFIX = ".motif_link_table.tsv.gz"

#: The alias table beside them, which is **plain**: three curated rows, where the tables
#: are bulk. Small metadata tables ship plain throughout this package and bulk tables
#: ship gzipped, and this directory is where both conventions meet.
ALIAS_FILE = "motif_name_alias.tsv"

#: Every link table's columns, in table order. Identical across tables so that two of
#: them concatenate into one frame that still says what each row came from.
LINK_COLUMNS = (
    "release",
    "species",
    "gene_id_stem",
    "symbol",
    "motif_id",
    "motif_name",
    "role",
    "partners",
    "motif_tax_ids",
    "is_cross_species",
    "total_information_content",
    "rank",
)

#: The alias table's columns, in table order.
ALIAS_COLUMNS = ("species", "motif_name_part", "gene_id_stem", "census_symbol", "note")

#: The two **Role**s, and the only two.
MONOMER, COMPLEX = "monomer", "complex"

#: How a flag column is spelled, which is how a census spells its own.
TRUE_FLAG, FALSE_FLAG = "yes", "no"

#: What separates one value from the next inside a cell.
VALUE_SEPARATOR = ";"

#: How a **Motif name** spells the join between the genes a complex names.
NAME_SEPARATOR = "::"

#: What each shipped table must contain: (genes, links, cross-species links, monomer
#: links). **Pinned because nothing else can be**: regenerating a table needs a download
#: CI has no network for, so these four numbers are what stands in for the check that
#: cannot exist. Keyed by ``(species slug, release)`` and checked to cover every table
#: that ships, so adding one without pinning it fails here rather than passing quietly.
#: Mouse is overwhelmingly cross-species — 732 of 896 on 2026 — which is ADR-0013's
#: coverage argument made concrete rather than an anomaly.
_PINNED: dict[tuple[str, str], tuple[int, int, int, int]] = {
    ("homo_sapiens", "2024"): (745, 946, 161, 794),
    ("homo_sapiens", "2026"): (876, 1085, 162, 931),
    ("mus_musculus", "2024"): (653, 851, 690, 702),
    ("mus_musculus", "2026"): (693, 896, 732, 745),
}

#: Genes whose only motifs are complexes, on the 2026 **Release**. Named because the
#: whole point of **Role** is that these are linked rather than reported motif-less: a
#: rule that merged ``FOS::JUN`` into a JUN motif would lose them silently.
_COMPLEX_ONLY_2026 = {"homo_sapiens": ("AHR", "DDIT3", "TAL1", "TLX1")}

#: The oncogenic fusion that names no gene. It stays unlinked by design, and asserting a
#: gene for it would be inventing one.
_FUSION = "EWSR1-FLI1"

#: A gene Lambert assessed and turned down, and the profile named for it. Unlinked
#: because only assessed-positive genes receive links — not because anything is missing.
_ASSESSED_NEGATIVE = ("homo_sapiens", "ENSG00000175387", "SMAD2")

#: The mouse profile that looks like a rename and is not one, with the part of its
#: **Motif name** an alias would have to carry. Unlinked because AnimalTFDB never
#: assessed the gene it names — see the test below, which records the two accessions.
_DIFFERENT_GENE = ("MA0611.3", "DUX")

#: The published span of a matrix's total **Information content**, in bits, with room on
#: both sides. A value outside it is a matrix read wrong rather than an unusual motif.
_IC_RANGE = (1.0, 60.0)


def _directory() -> Traversable:
    """Return the shipped link-table directory as a traversable resource."""
    return files("genome").joinpath(LINK_SUBDIR)


def _shipped_tables() -> tuple[tuple[str, str], ...]:
    """Return every shipped table as ``(species slug, release)``, sorted.

    Found by enumerating the directory and reading the two keys out of the file name, so
    a new release or a new species is a file drop rather than an edit here.
    """
    directory = _directory()
    if not directory.is_dir():
        return ()
    found = []
    for entry in directory.iterdir():
        if not entry.name.endswith(LINK_SUFFIX):
            continue
        slug, _, release = entry.name[: -len(LINK_SUFFIX)].rpartition(".")
        found.append((slug, release.removeprefix(RELEASE_PREFIX)))
    return tuple(sorted(found))


def _text(slug: str, release: str) -> str:
    """Return one shipped table unpacked, which is the one place the gzip is undone here."""
    name = f"{slug}.{RELEASE_PREFIX}{release}{LINK_SUFFIX}"
    return gzip.decompress(_directory().joinpath(name).read_bytes()).decode("utf-8")


def _read(slug: str, release: str) -> list[dict[str, str]]:
    """Return one shipped table's rows, read as the unquoted TSV it is."""
    lines = _text(slug, release).splitlines()
    return list(csv.DictReader(lines, delimiter="\t", quoting=csv.QUOTE_NONE))


def _header(slug: str, release: str) -> tuple[str, ...]:
    """Return one shipped table's header line, split on tabs."""
    return tuple(_text(slug, release).splitlines()[0].split("\t"))


def _aliases() -> list[dict[str, str]]:
    """Return the shipped alias table's rows."""
    text = _directory().joinpath(ALIAS_FILE).read_text(encoding="utf-8")
    return list(csv.DictReader(text.splitlines(), delimiter="\t", quoting=csv.QUOTE_NONE))


def _census(slug: str) -> TFGeneTable:
    """Return the census shipped for ``slug``, failing by name when none ships."""
    census = tf_gene_table(slug)
    assert census is not None, f"no census ships for {slug}, so its links answer for nothing"
    return census


def _by_stem(slug: str) -> dict[str, tuple[str, bool]]:
    """Return each **Gene id stem** the census assessed, with its symbol and its verdict."""
    frame = _census(slug).frame()
    return {
        str(stem): (str(symbol or ""), bool(flag))
        for stem, symbol, flag in zip(
            frame["gene_id_stem"], frame["symbol"], frame["is_tf"], strict=True
        )
    }


def _specificity_key(row: dict[str, str]) -> tuple[int, int, float, str]:
    """Return one shipped row's place under **Attribution specificity**, lowest first."""
    return (
        0 if row["role"] == MONOMER else 1,
        1 if row["is_cross_species"] == TRUE_FLAG else 0,
        -float(row["total_information_content"]),
        row["motif_id"],
    )


def _grouped(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    """Return ``rows`` grouped by **Gene id stem**, in file order within each gene."""
    genes: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        genes.setdefault(row["gene_id_stem"], []).append(row)
    return genes


_TABLES = _shipped_tables()


# ---------------------------------------------------------------------------------------
# What ships, and that all of it is pinned
# ---------------------------------------------------------------------------------------


def test_the_package_ships_link_tables_at_all() -> None:
    # The guard under every parametrized test below: with no files at all each of those
    # would collect zero cases and pass, which is the silent zero this module exists for.
    assert _TABLES


def test_every_shipped_table_has_its_counts_pinned() -> None:
    # A table added without pinned counts would be guarded by nothing at all — and pinned
    # counts are the only guard there can be, since regenerating one needs a download.
    assert sorted(_PINNED) == sorted(_TABLES)


def test_a_table_ships_for_every_census_and_every_release_the_package_prepares() -> None:
    # The linking layer covers what the package can load. Covering only the default
    # release would be a hole discovered at runtime, by whoever asked for the other one.
    assert sorted(_TABLES) == sorted(
        (slug, release) for slug in census_species() for release in JASPAR_RELEASES
    )


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_every_table_is_named_for_a_census_and_a_release(slug: str, release: str) -> None:
    # Two keys name one table, so both are in the file name.
    assert slug in census_species()
    assert release in JASPAR_RELEASES
    assert species_slug(slug) == slug


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_every_table_ships_gzipped_with_no_timestamp(slug: str, release: str) -> None:
    # Bulk data ships gzipped, as the censuses beside it do. The gzip header carries
    # ``mtime=0``, which is what lets two runs of the generator agree byte for byte —
    # without it every rebuild would diff whether or not the links changed.
    name = f"{slug}.{RELEASE_PREFIX}{release}{LINK_SUFFIX}"
    raw = _directory().joinpath(name).read_bytes()

    assert name.endswith(".tsv.gz")
    assert raw.startswith(b"\x1f\x8b")
    assert raw[4:8] == b"\x00\x00\x00\x00"
    assert _text(slug, release).endswith("\n")


def test_the_alias_table_ships_plain_beside_the_gzipped_tables() -> None:
    # The convention, and this directory is where its two halves meet: bulk tables
    # gzipped, small metadata tables plain — as `census_metadata.tsv` and the assembly
    # and annotation metadata tables already are. Three hand-curated rows are worth more
    # as a readable diff than as the couple of hundred bytes gzip would save.
    raw = _directory().joinpath(ALIAS_FILE).read_bytes()

    assert ALIAS_FILE.endswith(".tsv")
    assert not raw.startswith(b"\x1f\x8b")
    assert raw.endswith(b"\n")


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_every_table_carries_the_same_columns_in_the_same_order(slug: str, release: str) -> None:
    # Identical headers are what lets two tables concatenate into one frame.
    assert _header(slug, release) == LINK_COLUMNS


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_every_row_names_the_release_and_species_its_file_does(slug: str, release: str) -> None:
    # Carried on every row rather than left to the file name, so a concatenated frame
    # still says where each row came from.
    census = _census(slug)

    for row in _read(slug, release):
        assert row["release"] == release
        assert row["species"] == census.species


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_the_gene_link_and_species_counts_are_pinned(slug: str, release: str) -> None:
    rows = _read(slug, release)
    genes, links, cross, monomer = _PINNED[(slug, release)]

    assert len(rows) == links
    assert len({row["gene_id_stem"] for row in rows}) == genes
    assert sum(1 for row in rows if row["is_cross_species"] == TRUE_FLAG) == cross
    assert sum(1 for row in rows if row["role"] == MONOMER) == monomer


# ---------------------------------------------------------------------------------------
# Every link answers for a gene its census judged a transcription factor
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_every_linked_gene_is_in_its_census_and_assessed_positive(slug: str, release: str) -> None:
    # Only assessed-positive genes receive links: a gene the census assessed and turned
    # down is unlinked however a profile is named, and a gene it never assessed at all is
    # not this package's to link.
    assessed = _by_stem(slug)

    for row in _read(slug, release):
        stem = row["gene_id_stem"]
        assert stem in assessed, f"{stem} is linked and {slug} does not assess it"
        assert assessed[stem][1], f"{stem} is linked and {slug} did not judge it a TF"


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_every_link_carries_the_censuss_own_symbol_for_its_gene(slug: str, release: str) -> None:
    # The symbol on the row is the census's, never JASPAR's — they differ exactly where
    # the alias table says they do, and the row would be unreadable if it said which.
    assessed = _by_stem(slug)

    for row in _read(slug, release):
        assert row["symbol"] == assessed[row["gene_id_stem"]][0]


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_no_gene_names_one_motif_twice(slug: str, release: str) -> None:
    pairs = [(row["gene_id_stem"], row["motif_id"]) for row in _read(slug, release)]

    assert len(set(pairs)) == len(pairs)


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_every_motif_id_is_a_versioned_accession(slug: str, release: str) -> None:
    # A **Motif id** is versioned; a bare base id addresses whichever version a release
    # ships, which is not the same claim.
    for row in _read(slug, release):
        assert re.fullmatch(r"[A-Z]+\d+\.\d+", row["motif_id"]), row["motif_id"]


# ---------------------------------------------------------------------------------------
# Role: what the matrix is a motif of
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_every_role_is_one_of_the_two_declared(slug: str, release: str) -> None:
    assert {row["role"] for row in _read(slug, release)} <= {MONOMER, COMPLEX}


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_a_complex_names_at_least_one_partner_and_a_monomer_names_none(
    slug: str, release: str
) -> None:
    # The whole reason **Role** exists: a heterodimer matrix must never be read as a
    # monomer's, and a row that said `complex` with nothing beside it would be exactly
    # that reading with a label on it.
    for row in _read(slug, release):
        partners = [part for part in row["partners"].split(VALUE_SEPARATOR) if part]
        if row["role"] == COMPLEX:
            assert partners, f"{row['motif_id']} is a complex naming no partner"
        else:
            assert not partners, f"{row['motif_id']} is a monomer naming {partners}"


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_the_partners_are_the_other_genes_the_motif_name_names(slug: str, release: str) -> None:
    for row in _read(slug, release):
        parts = row["motif_name"].upper().split(NAME_SEPARATOR)
        partners = [part for part in row["partners"].split(VALUE_SEPARATOR) if part]

        assert len(partners) == len(parts) - 1
        assert set(partners) <= set(parts)
        assert (row["role"] == MONOMER) == (len(parts) == 1)


@pytest.mark.parametrize("slug", sorted(_COMPLEX_ONLY_2026))
def test_a_gene_whose_only_motifs_are_complexes_is_still_linked(slug: str) -> None:
    # AHR, DDIT3, TAL1 and TLX1 have no monomer matrix at all in the 2026 release. A join
    # that read `FOS::JUN` as a JUN motif would report them motif-less, which is the
    # failure **Role** was coined to make impossible — so they are named here.
    rows = _read(slug, "2026")
    by_symbol: dict[str, set[str]] = {}
    for row in rows:
        by_symbol.setdefault(row["symbol"], set()).add(row["role"])

    for symbol in _COMPLEX_ONLY_2026[slug]:
        assert symbol in by_symbol, f"{symbol} is reported motif-less"
        assert by_symbol[symbol] == {COMPLEX}


# ---------------------------------------------------------------------------------------
# Attribution specificity: the order one gene's links come back in
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_ranks_are_dense_from_one_within_a_gene(slug: str, release: str) -> None:
    for stem, rows in _grouped(_read(slug, release)).items():
        ranks = sorted(int(row["rank"]) for row in rows)

        assert ranks == list(range(1, len(rows) + 1)), f"{stem} ranks {ranks}"


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_ranks_encode_attribution_specificity(slug: str, release: str) -> None:
    # Four keys, so the order is total and stable: two machines and two releases mean the
    # same thing by "the motif for this factor". Re-sorting the shipped rows on the
    # shipped columns has to reproduce the shipped rank — which is also what makes the
    # rule checkable by a reader who never imports this package.
    for stem, rows in _grouped(_read(slug, release)).items():
        expected = sorted(rows, key=_specificity_key)

        assert [row["rank"] for row in expected] == [str(n) for n in range(1, len(rows) + 1)], stem


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_no_quality_score_is_shipped(slug: str, release: str) -> None:
    # JASPAR publishes none and this package invents none; the ordering states what a
    # matrix is attributable to and explicitly not which motif is better.
    for column in _header(slug, release):
        assert not any(word in column for word in ("score", "quality", "confidence", "weight"))


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_the_total_information_content_is_a_number_in_the_published_range(
    slug: str, release: str
) -> None:
    low, high = _IC_RANGE

    for row in _read(slug, release):
        value = float(row["total_information_content"])

        assert low < value < high, f"{row['motif_id']} carries {value} bits"


# ---------------------------------------------------------------------------------------
# Cross-species links are kept and marked (ADR-0013)
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_a_link_is_cross_species_exactly_when_the_profile_does_not_name_its_gene_species(
    slug: str, release: str
) -> None:
    # The flag is the only thing a caller whose question needs a species-matched profile
    # can filter on, so it has to be derivable from the row rather than trusted. A
    # profile the dump records no species for at all is marked cross-species: the row
    # cannot claim a match it has no evidence for.
    taxid = str(_census(slug).provenance.ncbi_taxid)

    for row in _read(slug, release):
        tax_ids = [part for part in row["motif_tax_ids"].split(VALUE_SEPARATOR) if part]
        assert row["is_cross_species"] in (TRUE_FLAG, FALSE_FLAG)
        assert (row["is_cross_species"] == TRUE_FLAG) == (taxid not in tax_ids)


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_every_tax_id_is_a_number_and_they_are_written_in_order(slug: str, release: str) -> None:
    for row in _read(slug, release):
        tax_ids = [part for part in row["motif_tax_ids"].split(VALUE_SEPARATOR) if part]

        assert all(part.isdigit() for part in tax_ids), row["motif_id"]
        assert tax_ids == sorted(tax_ids, key=int)
        assert len(set(tax_ids)) == len(tax_ids)


def test_cross_species_profiles_are_kept_rather_than_excluded() -> None:
    # Excluding them is the tidy reading a species column invites, and it costs mouse
    # most of its coverage — 732 of 896 links on 2026 (ADR-0013).
    rows = _read("mus_musculus", "2026")

    assert sum(1 for row in rows if row["is_cross_species"] == TRUE_FLAG) > len(rows) // 2


# ---------------------------------------------------------------------------------------
# The residuals: the alias table, and what stays unlinked on purpose
# ---------------------------------------------------------------------------------------


def test_the_alias_table_carries_the_columns_it_declares() -> None:
    text = _directory().joinpath(ALIAS_FILE).read_text(encoding="utf-8")

    assert tuple(text.splitlines()[0].split("\t")) == ALIAS_COLUMNS
    assert _aliases()


def test_every_alias_resolves_to_a_gene_present_in_its_census_and_assessed_positive() -> None:
    for alias in _aliases():
        slug = species_slug(alias["species"])
        assessed = _by_stem(slug)
        stem = alias["gene_id_stem"]

        assert slug in census_species(), f"{alias['species']} has no census"
        assert stem in assessed, f"{stem} is aliased and {slug} does not assess it"
        assert assessed[stem][1], f"{stem} is aliased and {slug} did not judge it a TF"


def test_every_alias_agrees_with_its_census_about_the_genes_symbol() -> None:
    # The recorded symbol is checked rather than trusted, so a row that went stale is an
    # error rather than a comment nobody reads.
    for alias in _aliases():
        assessed = _by_stem(species_slug(alias["species"]))

        assert assessed[alias["gene_id_stem"]][0] == alias["census_symbol"]


def test_the_alias_table_is_keyed_on_gene_id_and_not_on_symbol() -> None:
    # A symbol is the thing that moved — JASPAR's `SCAND3` is Lambert's `ZBED9` — so an
    # alias keyed on one would be keyed on the moving part. Each row names a gene id, and
    # its JASPAR-side name is one the census does not spell that way: an alias that
    # shadowed a symbol the census does carry would silently re-point a real gene.
    for alias in _aliases():
        assessed = _by_stem(species_slug(alias["species"]))
        symbols = {symbol.upper() for symbol, positive in assessed.values() if positive and symbol}

        assert alias["gene_id_stem"].startswith("ENS")
        assert alias["motif_name_part"] == alias["motif_name_part"].upper()
        assert alias["motif_name_part"] not in symbols
        assert alias["note"]


def test_every_alias_is_used_by_at_least_one_shipped_table() -> None:
    # An alias nothing uses is a claim nobody checks. Two of the three name profiles that
    # arrive only with the 2026 release, which is why this is asked of the tables
    # together rather than of each one.
    linked = {
        (row["species"], row["motif_name"].upper())
        for slug, release in _TABLES
        for row in _read(slug, release)
    }

    for alias in _aliases():
        assert any(
            species == alias["species"] and alias["motif_name_part"] in name.split(NAME_SEPARATOR)
            for species, name in linked
        ), f"{alias['motif_name_part']} is aliased and nothing uses it"


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_a_profile_that_names_no_gene_stays_unlinked(slug: str, release: str) -> None:
    # `EWSR1-FLI1` is an oncogenic fusion. It names no gene, and asserting one for it
    # would be inventing it — so it is absent by design rather than by omission.
    assert not [row for row in _read(slug, release) if row["motif_name"] == _FUSION]


def test_the_profile_whose_symbol_only_looks_like_a_rename_stays_unlinked() -> None:
    # JASPAR's `MA0611.3 Dux` is UniProt A1JVI8 — GeneID 664783, **MGI:3703875**.
    # AnimalTFDB's `Duxf3` is ENSMUSG00000075046 — GeneID 74399, **MGI:1921649**. Two MGI
    # accessions are two genes, so an alias joining them would assert an identity MGI
    # denies; the only thing linking them is a secondary EntrezGene xref Ensembl carries
    # on the `Duxf3` model, which is the Dux macrosatellite collapsing onto that locus.
    # So this profile names a gene the census never assessed and is unlinked for the same
    # structural reason `EWSR1-FLI1` is — a correct answer, not a gap. Pinned here
    # because the symbols look alike, and guessing at symbol history got `SCAND3` wrong
    # twice: an alias row for this is a regression, not a fix.
    motif_id, name_part = _DIFFERENT_GENE

    for _, release in [table for table in _TABLES if table[0] == "mus_musculus"]:
        assert not [row for row in _read("mus_musculus", release) if row["motif_id"] == motif_id]
    assert not [alias for alias in _aliases() if alias["motif_name_part"] == name_part]


def test_a_gene_its_census_turned_down_receives_no_link() -> None:
    # Lambert assessed SMAD2 and judged it not a sequence-specific TF, and JASPAR ships a
    # profile named for it. It is unlinked because only assessed-positive genes receive
    # links — not because anything failed to match, which is why it is not aliased.
    slug, stem, symbol = _ASSESSED_NEGATIVE
    assessed = _by_stem(slug)

    assert assessed[stem] == (symbol, False)
    for _, release in [table for table in _TABLES if table[0] == slug]:
        assert not [row for row in _read(slug, release) if row["gene_id_stem"] == stem]
