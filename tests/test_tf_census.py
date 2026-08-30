"""Tests for genome.tf.gene.census — the TF gene tables shipped inside the package.

The shipped files answer here; no fixture stands in for them. That is the point of
this module rather than an accident of it: building a census needs a download CI
cannot make, so a guard over what ships is the only check that can exist. Nothing
touches the network, and nothing writes into the package's own data directory — the
malformed cases are handed to the reader as text rather than laid down as files.

Every test that is about one census runs against both shipped species — human, Lambert's
dual-classified case, and mouse, AnimalTFDB's all-positive one — since there are only two
and each is structurally distinct. A third species is a file dropped into
``data/tf_gene/``, and dropping one in must not mean rewriting these.
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
from tests._tables import table_text

#: What each shipped census must contain: (genes assessed, genes judged TFs, distinct
#: **DBD family** values). Pinned because the biology is the shipped data's problem — a
#: regeneration that silently changed any of the three is exactly the drift there is no
#: other way to catch. Keyed by the species slug a census's file is named by, and checked
#: to cover every census that ships, so adding one without pinning it fails here rather
#: than passing quietly. The two family counts sit side by side and are never compared:
#: 75 values under Lambert's ``DBD``, 72 under AnimalTFDB's ``Family``, two vocabularies
#: deliberately not crosswalked (ADR-0014). AnimalTFDB lists none but the genes it
#: accepts, which is why mouse's second number is its first.
_PINNED: dict[str, tuple[int, int, int]] = {
    "homo_sapiens": (2765, 1639, 75),
    "mus_musculus": (1611, 1611, 72),
}

#: Where the **DBD family** sits in every census, which the uniform four fix.
_DBD_FAMILY = UNIFORM_COLUMNS.index("dbd_family")

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


def _families(slug: str) -> set[str]:
    """Return the distinct **DBD family** values one shipped census classifies genes under."""
    families = (row[_DBD_FAMILY] for row in _shipped(slug).rows)
    return {family for family in families if family}


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
    return table_text(UNIFORM_COLUMNS, *rows, header=header)


def _read(text: str) -> TFGeneTable:
    """Read ``text`` as a shipped census, under a name every message can name."""
    return _read_census(text, provenance=_provenance(), origin=f"tiny_beast{CENSUS_SUFFIX}")


# ---------------------------------------------------------------------------------------
# What ships: every census valid, and every one of them pinned
# ---------------------------------------------------------------------------------------


def test_the_package_ships_sorted_deduplicated_censuses_with_every_one_pinned() -> None:
    # The guard under every parametrized test below: with no files at all each of those
    # would collect zero cases and pass, which is the silent zero this module exists for.
    species = census_species()

    assert species
    assert list(species) == sorted(species)
    assert len(set(species)) == len(species)
    assert sorted(_PINNED) == sorted(species)


def test_a_census_that_lists_none_but_the_genes_it_accepts_is_read_as_a_census() -> None:
    # AnimalTFDB publishes no rejected set, so mouse is all yes. That is the publisher
    # saying nothing at all about the genes it left out, never this package writing a
    # `no` for them — so an all-positive census is a census's shape and not a defect.
    census = _read(_census_text("g1\tA\tyes\tX", "g2\tB\tyes\tY"))

    assert census.assessed_positive == census.gene_id_stems


def test_the_two_family_vocabularies_are_the_publishers_and_are_not_crosswalked() -> None:
    # Uniform in position and deliberately not in content (ADR-0014). Lambert spells the
    # ARID family `ARID/BRIGHT` and AnimalTFDB spells its own `ARID`, and nothing here
    # asserts they are one family: inventing an equivalence nobody has checked is worse
    # than shipping two vocabularies that each say who spelled them.
    human, mouse = _families("homo_sapiens"), _families("mus_musculus")

    assert "ARID/BRIGHT" in human
    assert "ARID/BRIGHT" not in mouse
    assert "ARID" in mouse
    assert "ARID" not in human
    assert _shipped("homo_sapiens").provenance.family_column == "DBD"
    assert _shipped("mus_musculus").provenance.family_column == "Family"


@pytest.mark.parametrize("slug", census_species())
def test_every_shipped_census_parses_with_well_formed_columns_and_unique_gene_ids(
    slug: str,
) -> None:
    # Loading is what validates, so a census that cannot be trusted raises here rather
    # than answering: this and the test below say the same invariants a loosened reader
    # would otherwise let slip past, each naming the offending file.
    census = _shipped(slug)

    assert species_slug(census.species) == slug
    assert census.provenance.file == f"{slug}{CENSUS_SUFFIX}"
    assert census.columns[: len(UNIFORM_COLUMNS)] == UNIFORM_COLUMNS
    assert len(set(census.columns)) == len(census.columns)

    # The publisher's own name, respelled — never the published spelling with its
    # capitals, spaces, punctuation or the trailing whitespace three of Lambert's carry.
    for column in census.columns[len(UNIFORM_COLUMNS) :]:
        assert column == column.strip(), f"{slug}: {column!r} carries whitespace"
        assert column.replace("_", "").isalnum(), f"{slug}: {column!r} is not snake_case"
        assert column.islower() or column.replace("_", "").isdigit(), (
            f"{slug}: {column!r} is not lower case"
        )

    stems = census.gene_id_stems
    assert all(stems), f"{slug}: a row carries no gene id stem"
    assert len(set(stems)) == len(stems), f"{slug}: a gene id stem is named twice"


@pytest.mark.parametrize("slug", census_species())
def test_the_pinned_counts_hold_and_a_rejected_gene_differs_from_an_unassessed_one(
    slug: str,
) -> None:
    census = _shipped(slug)
    rows, positive, families = _PINNED[slug]

    assert len(census) == rows
    assert len(census.assessed_positive) == positive
    assert len(_families(slug)) == families

    # The whole reason the rejected genes ship: a stem that is here and not in
    # assessed_positive was looked at and turned down, which is a verdict. A stem that
    # is in neither was never assessed, and only the first is the census speaking. The
    # containment is deliberately not strict: a publisher that lists none but the genes
    # it accepts has no rejected set to ship, and inventing one would be the fabrication
    # this distinction exists to prevent.
    assert set(census.assessed_positive) <= set(census.gene_id_stems)
    assert "no_such_gene" not in census.gene_id_stems


@pytest.mark.parametrize("slug", census_species())
def test_a_censuss_provenance_names_its_publisher_hashes_correctly_and_is_reachable(
    slug: str,
) -> None:
    census = _shipped(slug)
    provenance = census.provenance

    assert provenance.publisher
    assert provenance.pubmed_id > 0
    assert provenance.ncbi_taxid > 0
    assert provenance.source_url.startswith("http")
    assert provenance.family_column
    assert provenance.publisher in provenance.attribution()
    assert str(provenance.pubmed_id) in provenance.attribution()

    # Over the unpacked TSV and not the gzip around it (ADR-0006), so recompressing the
    # file elsewhere does not break a match that ought to hold.
    resource = files("genome").joinpath(CENSUS_SUBDIR, f"{slug}{CENSUS_SUFFIX}")
    unpacked = gzip.decompress(resource.read_bytes())
    assert hashlib.sha256(unpacked).hexdigest() == provenance.sha256

    # The species is read off the assembly and never passed in, so a census keyed to a
    # spelling no assembly uses would be unreachable, and the assembly table's own
    # spelling has to be the one that finds it. Equal and not identical: the two
    # spellings are two entries in the lookup's cache, and a census is frozen, so value
    # equality is the whole of what a caller can tell about them.
    assert census.species in _TABLE_SPECIES or census.species in _UNREACHABLE_SPECIES
    assert tf_gene_table(census.species) == census


@pytest.mark.parametrize("slug", census_species())
def test_the_frame_is_the_census_with_the_flag_read_as_a_boolean(slug: str) -> None:
    census = _shipped(slug)
    frame = census.frame()

    assert list(frame.columns) == list(census.columns)
    assert len(frame) == len(census)
    assert frame["is_tf"].dtype == bool
    assert int(frame["is_tf"].sum()) == len(census.assessed_positive)


def test_the_provenance_table_and_the_shipped_files_name_the_same_species() -> None:
    # A row with no file cannot be cited from anywhere; a file with no row cannot be
    # cited at all, which is the condition on redistributing a census here.
    assert sorted(species_slug(record.species) for record in census_metadata()) == sorted(
        census_species()
    )


# ---------------------------------------------------------------------------------------
# The raw absence: one function's ``None``, and nobody else's
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "species",
    [
        # One species the assembly metadata table names and no census answers for,
        # spelled as that table spells it, plus one nobody has ever registered.
        "Caenorhabditis elegans",
        "no such species",
    ],
)
def test_a_species_no_census_ships_for_is_none(species: str) -> None:
    assert tf_gene_table(species) is None


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
        ("  Mus   musculus  ", "mus_musculus"),
        ("Escherichia coli HT115", "escherichia_coli_ht115"),
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


def test_a_well_formed_census_reads_back_with_blanks_as_none_and_row_order_kept() -> None:
    census = _read(_census_text("g1\tA\tyes\tbHLH", "g2\tB\tno\t"))
    assert census.columns == UNIFORM_COLUMNS
    assert census.gene_id_stems == ("g1", "g2")
    assert census.assessed_positive == ("g1",)
    assert len(census) == 2

    # The publisher recorded nothing there, which is the reading every other table in
    # this package gives a blank cell.
    blank = _read(_census_text("g1\tA\tno\t"))
    assert blank.rows == (("g1", "A", "no", None),)

    ordered = _read(_census_text("g3\tC\tyes\tX", "g1\tA\tyes\tY", "g2\tB\tyes\tZ"))
    assert ordered.gene_id_stems == ("g3", "g1", "g2")


def test_a_header_that_does_not_lead_with_the_uniform_four_raises() -> None:
    with pytest.raises(TFGeneTableError, match="every census leads with"):
        _read(_census_text("g1\tA\tyes\tX", header="gene_id_stem\tsymbol\tis_tf"))


def test_a_malformed_censuss_shape_raises_and_names_the_defect() -> None:
    header = "\t".join([*UNIFORM_COLUMNS, "symbol"])
    with pytest.raises(TFGeneTableError, match="twice"):
        _read(_census_text("g1\tA\tyes\tX\tA", header=header))

    with pytest.raises(TFGeneTableError, match="line 3"):
        _read(_census_text("g1\tA\tyes\tX", "g2\tB\tyes"))


@pytest.mark.parametrize("flag", ["YES", ""])
def test_a_tf_flag_nobody_spells_that_way_raises(flag: str) -> None:
    with pytest.raises(TFGeneTableError, match=TRUE_CELL):
        _read(_census_text(f"g1\tA\t{flag}\tX"))


def test_a_missing_duplicated_or_absent_gene_identity_raises_as_a_bad_value() -> None:
    with pytest.raises(TFGeneTableError, match="gene id stem"):
        _read(_census_text("\tA\tyes\tX"))

    # A census reaches one verdict per gene, so two rows for one stem would let a caller
    # read either — which is the disagreement between two analyses this feature is about.
    with pytest.raises(TFGeneTableError, match="g1"):
        _read(_census_text("g1\tA\tyes\tX", "g1\tB\tno\t"))

    # Absence is spelled by shipping no file. A file that assessed nothing would be a
    # second spelling of it, and the one that reads as *this species has no TFs*.
    with pytest.raises(TFGeneTableError, match="no genes"):
        _read(_census_text())

    with pytest.raises(TFGeneTableError, match="empty"):
        _read("")

    # A caller catching LookupError for the absences above must not swallow a packaging
    # defect.
    assert issubclass(TFGeneTableError, ValueError)
    assert not issubclass(TFGeneTableError, LookupError)


# ---------------------------------------------------------------------------------------
# The provenance table beside the censuses
# ---------------------------------------------------------------------------------------


def test_the_provenance_table_reads_back_as_records_and_rejects_malformed_rows() -> None:
    records = census_metadata()
    assert records
    assert all(isinstance(record, CensusProvenance) for record in records)

    with pytest.raises(TFGeneTableError, match=CENSUS_METADATA_RESOURCE):
        _read_metadata("species\tpublisher\n", origin=CENSUS_METADATA_RESOURCE)

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
    blank_publisher = "Tiny beast\t1\ttiny.tsv.gz\t\tv1\t2\tDBD\thttps://example.org\tabc"
    with pytest.raises(TFGeneTableError, match="publisher"):
        _read_metadata(f"{header}\n{blank_publisher}\n", origin="census_metadata.tsv")

    bad_number = "Tiny beast\tmany\ttiny.tsv.gz\tSomeone\tv1\t2\tDBD\thttps://example.org\tabc"
    with pytest.raises(TFGeneTableError, match="ncbi_taxid"):
        _read_metadata(f"{header}\n{bad_number}\n", origin="census_metadata.tsv")


# ---------------------------------------------------------------------------------------
# The values themselves
# ---------------------------------------------------------------------------------------


def test_a_census_is_frozen_with_a_freshly_built_frame_and_fixed_flag_spellings() -> None:
    census = _shipped(census_species()[0])
    frame = census.frame()
    frame.loc[0, "symbol"] = "MUTATED"
    assert census.frame().loc[0, "symbol"] != "MUTATED"

    frozen = _read(_census_text("g1\tA\tyes\tX"))
    with pytest.raises(AttributeError):
        frozen.species = "Other beast"  # type: ignore[misc]

    assert {TRUE_CELL, FALSE_CELL} == {"yes", "no"}
