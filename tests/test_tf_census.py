"""Tests for genome.tf.gene.census — the TF gene tables shipped inside the package.

The shipped files answer here; no fixture stands in for them. That is the point of
this module rather than an accident of it: building a census needs a download CI
cannot make, so a guard over what ships is the only check that can exist. Nothing
touches the network, and nothing writes into the package's own data directory — the
malformed cases are handed to the reader as text rather than laid down as files.

Every test that is about one census iterates over *every* shipped census. A second
species is a file dropped into ``data/tf_gene/``, and dropping one in must not mean
rewriting these.
"""

from __future__ import annotations

import gzip
import hashlib
from importlib.resources import files

import pytest

from genome.metadata import assembly_table
from genome.tf.gene import (
    CENSUS_METADATA_RESOURCE,
    CENSUS_SUBDIR,
    CENSUS_SUFFIX,
    FALSE_CELL,
    TRUE_CELL,
    UNIFORM_COLUMNS,
    CensusProvenance,
    TFGeneTable,
    TFGeneTableError,
    census_metadata,
    census_species,
    species_slug,
    tf_gene_table,
)
from genome.tf.gene.census import _read_census, _read_metadata

#: What each shipped census must contain: (genes assessed, genes judged TFs). Pinned
#: because the biology is the shipped data's problem — a regeneration that silently
#: changed either number is exactly the drift there is no other way to catch. Keyed by
#: the species slug a census's file is named by, and checked to cover every census
#: that ships, so adding one without pinning it fails here rather than passing quietly.
_PINNED: dict[str, tuple[int, int]] = {"homo_sapiens": (2765, 1639)}

#: The species a census ships for that the assembly metadata table does not name. A
#: census keyed to a spelling no assembly uses could never be reached, since the
#: species is read off the assembly and never passed in.
_UNREACHABLE_SPECIES: tuple[str, ...] = ()

#: Every species the assembly metadata table spells, in its own spelling.
_TABLE_SPECIES = {record.species for record in assembly_table() if record.species}


def _shipped(slug: str) -> TFGeneTable:
    """Return the census shipped for ``slug``, failing by name when none ships."""
    census = tf_gene_table(slug)
    assert census is not None, f"no {slug}{CENSUS_SUFFIX} ships in the package"
    return census


def _provenance(**overrides: object) -> CensusProvenance:
    """Return one well-formed provenance record, with ``overrides`` applied."""
    fields: dict[str, object] = {
        "species": "Tiny beast",
        "ncbi_taxid": 1,
        "file": f"tiny_beast{CENSUS_SUFFIX}",
        "publisher": "Someone et al. 1999",
        "version": "v1",
        "pubmed_id": 2,
        "family_column": "DBD",
        "source_url": "https://example.org/census.csv",
        "sha256": "0" * 64,
    }
    fields.update(overrides)
    return CensusProvenance(**fields)  # type: ignore[arg-type]


def _census_text(*rows: str, header: str | None = None) -> str:
    """Return one census as text: a header of the uniform four, then ``rows``."""
    return "\n".join([header or "\t".join(UNIFORM_COLUMNS), *rows]) + "\n"


def _read(text: str) -> TFGeneTable:
    """Read ``text`` as a shipped census, under a name every message can name."""
    return _read_census(text, provenance=_provenance(), origin=f"tiny_beast{CENSUS_SUFFIX}")


# ---------------------------------------------------------------------------------------
# What ships: every census valid, and every one of them pinned
# ---------------------------------------------------------------------------------------


def test_the_package_ships_censuses_at_all() -> None:
    # The guard under every parametrized test below: with no files at all each of those
    # would collect zero cases and pass, which is the silent zero this module exists for.
    assert census_species()


def test_the_shipped_species_are_sorted_and_each_named_once() -> None:
    species = census_species()

    assert list(species) == sorted(species)
    assert len(set(species)) == len(species)


def test_every_shipped_census_has_its_counts_pinned() -> None:
    # A census added without a pinned count would be guarded by nothing at all.
    assert sorted(_PINNED) == sorted(census_species())


@pytest.mark.parametrize("slug", census_species())
def test_every_shipped_census_parses_and_names_its_species(slug: str) -> None:
    # Loading is what validates, so a census that cannot be trusted raises here rather
    # than answering. The assertions below say the same invariants again, so a loosened
    # reader is caught by a failure naming the offending file.
    census = _shipped(slug)

    assert species_slug(census.species) == slug
    assert census.provenance.file == f"{slug}{CENSUS_SUFFIX}"


@pytest.mark.parametrize("slug", census_species())
def test_every_shipped_census_leads_with_the_four_uniform_columns(slug: str) -> None:
    census = _shipped(slug)

    assert census.columns[: len(UNIFORM_COLUMNS)] == UNIFORM_COLUMNS
    assert len(set(census.columns)) == len(census.columns)


@pytest.mark.parametrize("slug", census_species())
def test_every_column_after_the_uniform_four_is_snake_case(slug: str) -> None:
    # The publisher's own name, respelled — never the published spelling with its
    # capitals, spaces, punctuation or the trailing whitespace three of Lambert's carry.
    for column in _shipped(slug).columns[len(UNIFORM_COLUMNS) :]:
        assert column == column.strip(), f"{slug}: {column!r} carries whitespace"
        assert column.replace("_", "").isalnum(), f"{slug}: {column!r} is not snake_case"
        assert column.islower() or column.replace("_", "").isdigit(), (
            f"{slug}: {column!r} is not lower case"
        )


@pytest.mark.parametrize("slug", census_species())
def test_every_gene_id_stem_is_unique_within_its_census(slug: str) -> None:
    stems = _shipped(slug).gene_id_stems

    assert all(stems), f"{slug}: a row carries no gene id stem"
    assert len(set(stems)) == len(stems), f"{slug}: a gene id stem is named twice"


@pytest.mark.parametrize("slug", census_species())
def test_the_row_and_assessed_positive_counts_are_pinned(slug: str) -> None:
    census = _shipped(slug)
    rows, positive = _PINNED[slug]

    assert len(census) == rows
    assert len(census.assessed_positive) == positive


@pytest.mark.parametrize("slug", census_species())
def test_a_rejected_gene_and_an_unassessed_one_stay_different_answers(slug: str) -> None:
    # The whole reason the rejected genes ship: a stem that is here and not in
    # assessed_positive was looked at and turned down, which is a verdict. A stem that
    # is in neither was never assessed, and only the first is the census speaking.
    census = _shipped(slug)

    assert set(census.assessed_positive) <= set(census.gene_id_stems)
    assert "no_such_gene" not in census.gene_id_stems


@pytest.mark.parametrize("slug", census_species())
def test_every_census_says_who_published_it_and_what_to_cite(slug: str) -> None:
    provenance = _shipped(slug).provenance

    assert provenance.publisher
    assert provenance.pubmed_id > 0
    assert provenance.ncbi_taxid > 0
    assert provenance.source_url.startswith("http")
    assert provenance.family_column
    assert provenance.publisher in provenance.attribution()
    assert str(provenance.pubmed_id) in provenance.attribution()


@pytest.mark.parametrize("slug", census_species())
def test_the_shipped_bytes_hash_to_what_the_provenance_pins(slug: str) -> None:
    # Over the unpacked TSV and not the gzip around it (ADR-0006), so recompressing the
    # file elsewhere does not break a match that ought to hold.
    resource = files("genome").joinpath(CENSUS_SUBDIR, f"{slug}{CENSUS_SUFFIX}")
    unpacked = gzip.decompress(resource.read_bytes())

    assert hashlib.sha256(unpacked).hexdigest() == _shipped(slug).provenance.sha256


def test_the_provenance_table_and_the_shipped_files_name_the_same_species() -> None:
    # A row with no file cannot be cited from anywhere; a file with no row cannot be
    # cited at all, which is the condition on redistributing a census here.
    assert sorted(species_slug(record.species) for record in census_metadata()) == sorted(
        census_species()
    )


@pytest.mark.parametrize("slug", census_species())
def test_every_censused_species_is_one_the_assembly_table_names(slug: str) -> None:
    # The species is read off the assembly and never passed in, so a census keyed to a
    # spelling no assembly uses would be unreachable.
    census = _shipped(slug)

    assert census.species in _TABLE_SPECIES or census.species in _UNREACHABLE_SPECIES


def test_the_species_spelling_the_assembly_table_uses_finds_its_census() -> None:
    assert tf_gene_table("Homo sapiens") is not None
    assert tf_gene_table("Homo sapiens") == tf_gene_table("homo_sapiens")


# ---------------------------------------------------------------------------------------
# The raw absence: one function's ``None``, and nobody else's
# ---------------------------------------------------------------------------------------


def test_a_species_no_census_ships_for_is_none() -> None:
    assert tf_gene_table("Caenorhabditis elegans") is None
    assert tf_gene_table("no such species") is None


def test_a_species_name_that_is_a_path_reaches_nothing() -> None:
    # Names are looked up among what ships rather than joined onto the resource
    # directory, so a name shaped like a path finds nothing instead of walking out of it.
    assert tf_gene_table("../assembly_metadata") is None
    assert tf_gene_table("") is None


# ---------------------------------------------------------------------------------------
# The species slug, which is the filename convention
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("species", "slug"),
    [
        ("Homo sapiens", "homo_sapiens"),
        ("homo_sapiens", "homo_sapiens"),
        ("  Mus   musculus  ", "mus_musculus"),
        ("Escherichia coli HT115", "escherichia_coli_ht115"),
        ("Saccharomyces cerevisiae", "saccharomyces_cerevisiae"),
    ],
)
def test_a_species_slugs_the_same_way_however_it_is_spelled(species: str, slug: str) -> None:
    assert species_slug(species) == slug


def test_slugging_a_slug_changes_nothing() -> None:
    for slug in census_species():
        assert species_slug(slug) == slug


# ---------------------------------------------------------------------------------------
# A file that is there and cannot be read is a packaging defect
# ---------------------------------------------------------------------------------------


def test_a_well_formed_census_reads_back_as_its_rows() -> None:
    census = _read(_census_text("g1\tA\tyes\tbHLH", "g2\tB\tno\t"))

    assert census.columns == UNIFORM_COLUMNS
    assert census.gene_id_stems == ("g1", "g2")
    assert census.assessed_positive == ("g1",)
    assert len(census) == 2


def test_a_blank_cell_reads_back_as_none() -> None:
    # The publisher recorded nothing there, which is the reading every other table in
    # this package gives a blank cell.
    census = _read(_census_text("g1\tA\tno\t"))

    assert census.rows == (("g1", "A", "no", None),)


def test_a_census_keeps_the_row_order_of_the_file() -> None:
    census = _read(_census_text("g3\tC\tyes\tX", "g1\tA\tyes\tY", "g2\tB\tyes\tZ"))

    assert census.gene_id_stems == ("g3", "g1", "g2")


@pytest.mark.parametrize(
    "header",
    [
        "gene_id\tsymbol\tis_tf\tdbd_family",
        "symbol\tgene_id_stem\tis_tf\tdbd_family",
        "gene_id_stem\tsymbol\tis_tf",
    ],
)
def test_a_header_that_does_not_lead_with_the_uniform_four_raises(header: str) -> None:
    with pytest.raises(TFGeneTableError, match="every census leads with"):
        _read(_census_text("g1\tA\tyes\tX", header=header))


def test_a_column_named_twice_raises() -> None:
    header = "\t".join([*UNIFORM_COLUMNS, "symbol"])

    with pytest.raises(TFGeneTableError, match="twice"):
        _read(_census_text("g1\tA\tyes\tX\tA", header=header))


def test_a_row_with_the_wrong_number_of_cells_raises() -> None:
    with pytest.raises(TFGeneTableError, match="line 3"):
        _read(_census_text("g1\tA\tyes\tX", "g2\tB\tyes"))


@pytest.mark.parametrize("flag", ["Yes", "YES", "true", "1", ""])
def test_a_tf_flag_nobody_spells_that_way_raises(flag: str) -> None:
    with pytest.raises(TFGeneTableError, match=TRUE_CELL):
        _read(_census_text(f"g1\tA\t{flag}\tX"))


def test_a_row_with_no_gene_id_stem_raises() -> None:
    with pytest.raises(TFGeneTableError, match="gene id stem"):
        _read(_census_text("\tA\tyes\tX"))


def test_the_same_gene_id_stem_twice_raises() -> None:
    # A census reaches one verdict per gene, so two rows for one stem would let a caller
    # read either — which is the disagreement between two analyses this feature is about.
    with pytest.raises(TFGeneTableError, match="g1"):
        _read(_census_text("g1\tA\tyes\tX", "g1\tB\tno\t"))


def test_a_census_with_a_header_and_no_genes_raises() -> None:
    # Absence is spelled by shipping no file. A file that assessed nothing would be a
    # second spelling of it, and the one that reads as *this species has no TFs*.
    with pytest.raises(TFGeneTableError, match="no genes"):
        _read(_census_text())


def test_an_empty_census_raises() -> None:
    with pytest.raises(TFGeneTableError, match="empty"):
        _read("")


def test_a_broken_census_is_a_bad_value_and_not_a_lookup() -> None:
    # A caller catching LookupError for the absences above it must not swallow a
    # packaging defect.
    assert issubclass(TFGeneTableError, ValueError)
    assert not issubclass(TFGeneTableError, LookupError)


# ---------------------------------------------------------------------------------------
# The provenance table beside the censuses
# ---------------------------------------------------------------------------------------


def test_the_provenance_table_reads_back_as_records() -> None:
    records = census_metadata()

    assert records
    assert all(isinstance(record, CensusProvenance) for record in records)


def test_a_provenance_header_that_is_not_the_tables_raises() -> None:
    with pytest.raises(TFGeneTableError, match=CENSUS_METADATA_RESOURCE):
        _read_metadata("species\tpublisher\n", origin=CENSUS_METADATA_RESOURCE)


def test_a_blank_provenance_cell_names_its_column() -> None:
    header = "\t".join(
        (
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
    )
    row = "Tiny beast\t1\ttiny.tsv.gz\t\tv1\t2\tDBD\thttps://example.org\tabc"

    with pytest.raises(TFGeneTableError, match="publisher"):
        _read_metadata(f"{header}\n{row}\n", origin="census_metadata.tsv")


def test_a_provenance_number_that_is_not_a_number_names_its_column() -> None:
    header = "\t".join(
        (
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
    )
    row = "Tiny beast\tmany\ttiny.tsv.gz\tSomeone\tv1\t2\tDBD\thttps://example.org\tabc"

    with pytest.raises(TFGeneTableError, match="ncbi_taxid"):
        _read_metadata(f"{header}\n{row}\n", origin="census_metadata.tsv")


# ---------------------------------------------------------------------------------------
# The values themselves
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize("slug", census_species())
def test_the_frame_is_the_census_with_the_flag_read_as_a_boolean(slug: str) -> None:
    census = _shipped(slug)
    frame = census.frame()

    assert list(frame.columns) == list(census.columns)
    assert len(frame) == len(census)
    assert frame["is_tf"].dtype == bool
    assert int(frame["is_tf"].sum()) == len(census.assessed_positive)


def test_the_frame_is_built_fresh_so_mutating_it_cannot_reach_the_census() -> None:
    census = _shipped(census_species()[0])
    frame = census.frame()
    frame.loc[0, "symbol"] = "MUTATED"

    assert census.frame().loc[0, "symbol"] != "MUTATED"


def test_a_census_is_frozen() -> None:
    census = _read(_census_text("g1\tA\tyes\tX"))

    with pytest.raises(AttributeError):
        census.species = "Other beast"  # type: ignore[misc]


def test_the_flag_spellings_are_the_two_a_census_uses() -> None:
    assert {TRUE_CELL, FALSE_CELL} == {"yes", "no"}
