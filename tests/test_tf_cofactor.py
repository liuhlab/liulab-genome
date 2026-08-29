"""Tests for genome.tf.cofactor — the cofactor tables shipped inside the package.

The shipped files answer here; no fixture stands in for them. That is the point of this
module rather than an accident of it: building a cofactor table needs downloads CI
cannot make, so a guard over what ships is the only check that can exist. Nothing
touches the network, and nothing writes into the package's own data directory — the
malformed cases are handed to the public parse entry point as text rather than laid down
as files, which is what that function is public for.

Every test that is about one shipped table iterates over *every* shipped table. A second
species is a file dropped into ``data/tf_cofactor/``, and dropping one in must not mean
rewriting these.
"""

from __future__ import annotations

import gzip
import hashlib
from collections.abc import Mapping, Sequence
from importlib.resources import files

import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.metadata import assembly_table
from genome.tf import VALUE_SEPARATOR
from genome.tf.cofactor import (
    ANIMALTFDB,
    BOTH,
    CITED_SOURCES,
    COFACTOR_SUBDIR,
    COFACTOR_SUFFIX,
    EPIFACTORS,
    FALSE_CELL,
    HGNC,
    SOURCES,
    TRUE_CELL,
    UNIFORM_COLUMNS,
    CofactorProvenance,
    CofactorSource,
    CofactorTable,
    CofactorTableError,
    cofactor_metadata,
    cofactor_species,
    cofactor_table,
    parse_cofactor_table,
    species_slug,
)

#: What each shipped table must contain, keyed by the species slug its file is named by
#: and checked to cover every table that ships, so adding one without pinning it fails
#: here rather than passing quietly. Pinned because the biology is the shipped data's
#: problem — a regeneration that silently changed a count is exactly the drift there is
#: no other way to catch. The three per-source counts are keyed by the closed vocabulary
#: itself, so a species listed by two publishers is an entry with different numbers
#: rather than an entry of a different shape. The family counts are AnimalTFDB's own
#: gene-list spellings and are never compared with another publisher's (ADR-0014).
_PINNED: dict[str, dict[str, int]] = {
    "caenorhabditis_elegans": {
        "rows": 317,
        "families": 57,
        "categories": 6,
        ANIMALTFDB: 317,
        EPIFACTORS: 0,
        BOTH: 0,
    },
    "homo_sapiens": {
        "rows": 1466,
        "families": 85,
        "categories": 6,
        ANIMALTFDB: 670,
        EPIFACTORS: 442,
        BOTH: 354,
    },
    "mus_musculus": {
        "rows": 970,
        "families": 84,
        "categories": 6,
        ANIMALTFDB: 970,
        EPIFACTORS: 0,
        BOTH: 0,
    },
}

#: AnimalTFDB's own two columns, which every shipped table carries today.
_FAMILY, _CATEGORY = "animaltfdb_family", "animaltfdb_category"

#: EpiFactors' own four columns, which only the human table carries.
_FUNCTION = "epifactors_function"
_TARGET = "epifactors_target"
_MODIFICATION = "epifactors_modification"
_COMPLEX = "epifactors_complex_name"
_EPIFACTORS_COLUMNS = (_FUNCTION, _TARGET, _MODIFICATION, _COMPLEX)

#: How many *atomic* values each EpiFactors vocabulary holds once a multi-valued cell is
#: split — the number that matters to a caller grouping by one of them, and not the
#: number of distinct cells, which counts every combination separately. Pinned because a
#: split on the wrong character shows up here first: splitting the complex names on a
#: bare comma rather than on this package's separator invents complexes called ``2`` and
#: ``4`` out of ``COMPASS-like MLL1,2`` and ``COMPASS-like MLL3,4``.
_EPIFACTORS_VOCABULARIES = {_FUNCTION: 19, _TARGET: 7, _MODIFICATION: 24, _COMPLEX: 70}

#: The five genes EpiFactors gives two rows each — not duplicates, but a histone
#: annotation on one row and an RNA one on the other — by the stem HGNC gives them.
#: They ship as one row apiece with their cells unioned; ``ALKBH1`` is the clearest,
#: carrying a function from each of its two rows.
_DOUBLE_ROWED = {
    "ENSG00000100601": "ALKBH1",
    "ENSG00000204389": "HSPA1A",
    "ENSG00000204388": "HSPA1B",
    "ENSG00000135372": "NAT10",
    "ENSG00000011304": "PTBP1",
}

#: Three genes EpiFactors still spells by a name HGNC has retired, and the spelling HGNC
#: approves today. The executable record that the join goes through the HGNC id: matching
#: on ``ACINU``, ``ARNTL`` or ``C11orf30`` would key these genes wrongly or drop them,
#: and 31 of EpiFactors' 801 rows are in that position.
_RENAMED_BY_HGNC = {
    "ENSG00000100813": "ACIN1",
    "ENSG00000133794": "BMAL1",
    "ENSG00000158636": "EMSY",
}

#: How many human stems are both a Lambert-positive **TF gene** and a **Transcription
#: cofactor**, and how many arrive from each publisher's side (28 from both, so the two
#: do not add up to the total). ADR-0016's central number: the two lists overlap, and a
#: caller who unions the two answers double-counts these genes.
_DUAL_CLASSIFIED, _DUAL_ANIMALTFDB, _DUAL_EPIFACTORS = 151, 57, 122

#: Three of them, named so that a failure says which genes moved rather than only how
#: many. TBP is a cofactor on the AnimalTFDB side alone; KMT2A and DNMT1 on both sides.
_DUAL_GENES = {
    "ENSG00000112592": "TBP",
    "ENSG00000118058": "KMT2A",
    "ENSG00000130816": "DNMT1",
}

#: The species a cofactor table ships for that the assembly metadata table does not name.
#: A table keyed to a spelling no assembly uses could never be reached, since the species
#: is read off the assembly and never passed in.
_UNREACHABLE_SPECIES: tuple[str, ...] = ()

#: Every species the assembly metadata table spells, in its own spelling.
_TABLE_SPECIES = {record.species for record in assembly_table() if record.species}

#: A cell of a well-formed table: no tab and no newline, and blank is legal.
_cell = st.text(alphabet="abAB01 -_/", max_size=6)

#: A gene id stem: never blank, and unique within a table.
_stem = st.text(alphabet="AB01", min_size=1, max_size=6)


def _shipped(slug: str) -> CofactorTable:
    """Return the cofactor table shipped for ``slug``, failing by name when none ships."""
    table = cofactor_table(slug)
    assert table is not None, f"no {slug}{COFACTOR_SUFFIX} ships in the package"
    return table


def _column(slug: str, column: str) -> tuple[str | None, ...]:
    """Return one shipped table's cells under ``column``, in row order."""
    table = _shipped(slug)
    index = table.columns.index(column)
    return tuple(row[index] for row in table.rows)


def _values(slug: str, column: str) -> set[str]:
    """Return the distinct non-blank values one shipped table records under ``column``."""
    return {cell for cell in _column(slug, column) if cell}


def _atomic(slug: str, column: str) -> set[str]:
    """Return one column's vocabulary once every multi-valued cell is split."""
    return {
        value for cell in _column(slug, column) if cell for value in cell.split(VALUE_SEPARATOR)
    }


def _by_stem(slug: str) -> dict[str, dict[str, str | None]]:
    """Return one shipped table keyed by **Gene id stem**, each row a column-to-cell map."""
    table = _shipped(slug)
    return {row[0] or "": dict(zip(table.columns, row, strict=True)) for row in table.rows}


def _filled(row: Mapping[str, str | None], columns: Sequence[str]) -> bool:
    """Return whether one publisher's group of columns says anything on this row."""
    return any(row[column] for column in columns)


def _source(**overrides: object) -> CofactorSource:
    """Return one well-formed source record, with ``overrides`` applied."""
    fields: dict[str, object] = {
        "species": "Tiny beast",
        "source": ANIMALTFDB,
        "publisher": "Someone et al. 1999",
        "version": "v1",
        "pubmed_id": 2,
        "source_url": "https://example.org/beast",
    }
    fields.update(overrides)
    return CofactorSource(**fields)  # type: ignore[arg-type]


def _provenance(**overrides: object) -> CofactorProvenance:
    """Return one well-formed provenance record, with ``overrides`` applied."""
    fields: dict[str, object] = {
        "species": "Tiny beast",
        "ncbi_taxid": 1,
        "file": f"tiny_beast{COFACTOR_SUFFIX}",
        "sha256": "0" * 64,
        "sources": (_source(),),
    }
    fields.update(overrides)
    return CofactorProvenance(**fields)  # type: ignore[arg-type]


def _table_text(*rows: str, header: str | None = None) -> str:
    """Return one cofactor table as text: a header of the uniform four, then ``rows``."""
    return "\n".join([header or "\t".join(UNIFORM_COLUMNS), *rows]) + "\n"


def _read(text: str) -> CofactorTable:
    """Read ``text`` as a shipped table, under a name every message can name."""
    return parse_cofactor_table(text, provenance=_provenance())


# ---------------------------------------------------------------------------------------
# What ships: every table valid, and every one of them pinned
# ---------------------------------------------------------------------------------------


def test_the_package_ships_cofactor_tables_at_all() -> None:
    # The guard under every parametrized test below: with no files at all each of those
    # would collect zero cases and pass, which is the silent zero this module exists for.
    assert cofactor_species()


def test_the_shipped_species_are_sorted_and_each_named_once() -> None:
    species = cofactor_species()

    assert list(species) == sorted(species)
    assert len(set(species)) == len(species)


def test_every_shipped_table_has_its_counts_pinned() -> None:
    # A table added without a pinned count would be guarded by nothing at all.
    assert sorted(_PINNED) == sorted(cofactor_species())


@pytest.mark.parametrize("slug", cofactor_species())
def test_every_shipped_table_parses_and_names_its_species(slug: str) -> None:
    # Loading is what validates, so a table that cannot be trusted raises here rather
    # than answering. The assertions below say the same invariants again, so a loosened
    # reader is caught by a failure naming the offending file.
    table = _shipped(slug)

    assert species_slug(table.species) == slug
    assert table.provenance.file == f"{slug}{COFACTOR_SUFFIX}"


@pytest.mark.parametrize("slug", cofactor_species())
def test_every_shipped_table_leads_with_the_four_uniform_columns(slug: str) -> None:
    table = _shipped(slug)

    assert table.columns[: len(UNIFORM_COLUMNS)] == UNIFORM_COLUMNS
    assert len(set(table.columns)) == len(table.columns)


@pytest.mark.parametrize("slug", cofactor_species())
def test_every_column_after_the_uniform_four_is_namespaced_snake_case(slug: str) -> None:
    # The publisher's own column under a snake_case name prefixed by which publisher's
    # vocabulary it is, so two publishers' columns never collide and a reader grouping by
    # one of them can see whose values they are (ADR-0014).
    for column in _shipped(slug).columns[len(UNIFORM_COLUMNS) :]:
        assert column == column.strip(), f"{slug}: {column!r} carries whitespace"
        assert column.replace("_", "").isalnum(), f"{slug}: {column!r} is not snake_case"
        assert column.islower(), f"{slug}: {column!r} is not lower case"
        assert any(column.startswith(f"{source}_") for source in SOURCES), (
            f"{slug}: {column!r} says whose vocabulary it is nowhere in its name"
        )


@pytest.mark.parametrize("slug", cofactor_species())
def test_every_gene_id_stem_is_unique_within_its_table(slug: str) -> None:
    stems = _shipped(slug).gene_id_stems

    assert all(stems), f"{slug}: a row carries no gene id stem"
    assert len(set(stems)) == len(stems), f"{slug}: a gene id stem is named twice"


@pytest.mark.parametrize("slug", cofactor_species())
def test_the_row_family_and_category_counts_are_pinned(slug: str) -> None:
    table = _shipped(slug)
    pinned = _PINNED[slug]

    assert len(table) == pinned["rows"]
    assert len(_values(slug, _FAMILY)) == pinned["families"]
    assert len(_values(slug, _CATEGORY)) == pinned["categories"]


@pytest.mark.parametrize("slug", cofactor_species())
def test_the_per_source_counts_are_pinned_and_the_vocabulary_is_closed(slug: str) -> None:
    sources = _column(slug, "source")

    assert set(sources) <= set(SOURCES), f"{slug}: a source outside the vocabulary"
    for source in SOURCES:
        assert sources.count(source) == _PINNED[slug][source]


@pytest.mark.parametrize("slug", cofactor_species())
def test_every_row_that_carries_a_family_carries_its_category(slug: str) -> None:
    # The category is joined onto the family as the table is built, and a family that
    # reached no category fails that build — so a blank category beside a filled family
    # would mean the join silently blanked a column instead.
    for family, category in zip(_column(slug, _FAMILY), _column(slug, _CATEGORY), strict=True):
        assert bool(family) == bool(category), f"{slug}: {family!r} sits in no category"


@pytest.mark.parametrize("slug", cofactor_species())
def test_every_listed_gene_is_one_the_publisher_calls_a_cofactor(slug: str) -> None:
    # No publisher shipping here releases a rejected set, so the flag reads yes on every
    # row. The column is kept anyway: presence in the file is not the verdict, and a
    # source that did record a rejection would ship `no` rows into this same format.
    table = _shipped(slug)

    assert set(table.cofactor_stems) <= set(table.gene_id_stems)
    assert table.cofactor_stems == table.gene_id_stems


@pytest.mark.parametrize("slug", cofactor_species())
def test_every_table_says_who_published_it_and_what_to_cite(slug: str) -> None:
    provenance = _shipped(slug).provenance

    assert provenance.ncbi_taxid > 0
    assert provenance.sources
    for source in provenance.sources:
        assert source.publisher
        assert source.pubmed_id > 0
        assert source.source in CITED_SOURCES
        assert source.source_url.startswith("http")
        assert source.publisher in provenance.attribution()
        assert str(source.pubmed_id) in provenance.attribution()


@pytest.mark.parametrize("slug", cofactor_species())
def test_one_attribution_line_is_rendered_per_species(slug: str) -> None:
    # One line however many publishers contributed, so the CLI, a notebook and an error
    # message all say it the same way.
    attribution = _shipped(slug).provenance.attribution()

    assert attribution
    assert "\n" not in attribution


@pytest.mark.parametrize("slug", cofactor_species())
def test_the_shipped_bytes_hash_to_what_the_provenance_pins(slug: str) -> None:
    # Over the unpacked TSV and not the gzip around it (ADR-0006), so recompressing the
    # file elsewhere does not break a match that ought to hold.
    resource = files("genome").joinpath(COFACTOR_SUBDIR, f"{slug}{COFACTOR_SUFFIX}")
    unpacked = gzip.decompress(resource.read_bytes())

    assert hashlib.sha256(unpacked).hexdigest() == _shipped(slug).provenance.sha256


@pytest.mark.parametrize("slug", cofactor_species())
def test_every_species_with_a_table_is_one_the_assembly_table_names(slug: str) -> None:
    # The species is read off the assembly and never passed in, so a table keyed to a
    # spelling no assembly uses would be unreachable.
    table = _shipped(slug)

    assert table.species in _TABLE_SPECIES or table.species in _UNREACHABLE_SPECIES


@pytest.mark.parametrize("slug", cofactor_species())
def test_the_species_spelling_the_assembly_table_uses_finds_its_table(slug: str) -> None:
    # An assembly carries its species in the metadata table's spelling and never a slug,
    # so that spelling has to reach the table the slug names. Equal and not identical:
    # the two spellings are two entries in the lookup's cache, and a table is frozen, so
    # value equality is the whole of what a caller can tell about them.
    table = _shipped(slug)

    assert cofactor_table(table.species) == table


def test_the_worm_table_ships_although_no_publisher_censused_worm_transcription_factors() -> None:
    # The asymmetry stated in ATTRIBUTION.md, pinned so that removing the worm table to
    # "match" the censuses is a failure rather than a tidy-up: absence is the publishers'
    # shape, and AnimalTFDB assessed worm cofactors while nobody has published its TFs.
    from genome.tf.gene import tf_gene_table

    assert cofactor_table("Caenorhabditis elegans") is not None
    assert tf_gene_table("Caenorhabditis elegans") is None


# ---------------------------------------------------------------------------------------
# The raw absence: one function's ``None``, and nobody else's
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "species",
    [
        # Two species the assembly metadata table names and no cofactor table answers
        # for, spelled as that table spells them, plus one nobody has ever registered.
        "Saccharomyces cerevisiae",
        "Escherichia coli HT115",
        "no such species",
    ],
)
def test_a_species_no_table_ships_for_is_none(species: str) -> None:
    assert cofactor_table(species) is None


def test_a_species_name_that_is_a_path_reaches_nothing() -> None:
    # Names are looked up among what ships rather than joined onto the resource
    # directory, so a name shaped like a path finds nothing instead of walking out of it.
    assert cofactor_table("../assembly_metadata") is None
    assert cofactor_table("../tf_gene/mus_musculus") is None
    assert cofactor_table("") is None


def test_the_slugger_is_the_one_the_censuses_use() -> None:
    # Imported rather than written a third time, so the two shipped-data directories
    # cannot drift into naming one species two ways.
    from genome.tf.gene import species_slug as census_species_slug

    assert species_slug is census_species_slug
    assert species_slug("Caenorhabditis elegans") == "caenorhabditis_elegans"


# ---------------------------------------------------------------------------------------
# A file that is there and cannot be read is a packaging defect
# ---------------------------------------------------------------------------------------


def test_a_well_formed_table_reads_back_as_its_rows() -> None:
    table = _read(_table_text("g1\tA\tyes\tanimaltfdb", "g2\tB\tno\tepifactors"))

    assert table.columns == UNIFORM_COLUMNS
    assert table.gene_id_stems == ("g1", "g2")
    assert table.cofactor_stems == ("g1",)
    assert len(table) == 2


def test_a_blank_cell_reads_back_as_none() -> None:
    # The publisher recorded nothing there, which is the reading every other table in
    # this package gives a blank cell.
    header = "\t".join([*UNIFORM_COLUMNS, _FAMILY])
    table = _read(_table_text("g1\tA\tyes\tanimaltfdb\t", header=header))

    assert table.rows == (("g1", "A", "yes", "animaltfdb", None),)


def test_a_table_keeps_the_row_order_of_the_file() -> None:
    table = _read(
        _table_text("g3\tC\tyes\tboth", "g1\tA\tyes\tanimaltfdb", "g2\tB\tyes\tepifactors")
    )

    assert table.gene_id_stems == ("g3", "g1", "g2")


def test_the_origin_defaults_to_the_file_the_provenance_names() -> None:
    with pytest.raises(CofactorTableError, match=f"tiny_beast{COFACTOR_SUFFIX}"):
        _read(_table_text())


@pytest.mark.parametrize(
    "header",
    [
        "gene_id\tsymbol\tis_cofactor\tsource",
        "symbol\tgene_id_stem\tis_cofactor\tsource",
        "gene_id_stem\tsymbol\tsource\tis_cofactor",
        "gene_id_stem\tsymbol\tis_cofactor",
    ],
)
def test_a_header_that_does_not_lead_with_the_uniform_four_raises(header: str) -> None:
    with pytest.raises(CofactorTableError, match="every cofactor table leads with"):
        _read(_table_text("g1\tA\tyes\tanimaltfdb", header=header))


def test_a_column_named_twice_raises() -> None:
    header = "\t".join([*UNIFORM_COLUMNS, "symbol"])

    with pytest.raises(CofactorTableError, match="twice"):
        _read(_table_text("g1\tA\tyes\tanimaltfdb\tA", header=header))


def test_a_row_with_the_wrong_number_of_cells_raises() -> None:
    with pytest.raises(CofactorTableError, match="line 3"):
        _read(_table_text("g1\tA\tyes\tanimaltfdb", "g2\tB\tyes"))


def test_a_cell_carrying_a_tab_raises() -> None:
    # A cofactor table is a plain TSV with no quoting, so a tab inside a cell reaches the
    # reader as a row with one cell too many rather than as something to parse around.
    with pytest.raises(CofactorTableError, match="carrying a tab"):
        _read(_table_text("g1\tA\tB\tyes\tanimaltfdb"))


@pytest.mark.parametrize("source", ["ANIMALTFDB", "animalTFDB", "lambert2018", "", "both;epi"])
def test_a_source_outside_the_closed_vocabulary_raises(source: str) -> None:
    with pytest.raises(CofactorTableError, match="the vocabulary is"):
        _read(_table_text(f"g1\tA\tyes\t{source}"))


@pytest.mark.parametrize("flag", ["Yes", "YES", "true", "1", ""])
def test_a_cofactor_flag_nobody_spells_that_way_raises(flag: str) -> None:
    with pytest.raises(CofactorTableError, match=TRUE_CELL):
        _read(_table_text(f"g1\tA\t{flag}\tanimaltfdb"))


def test_a_row_with_no_gene_id_stem_raises() -> None:
    with pytest.raises(CofactorTableError, match="gene id stem"):
        _read(_table_text("\tA\tyes\tanimaltfdb"))


def test_the_same_gene_id_stem_twice_raises() -> None:
    # One row per gene, so two rows for one stem would let a caller read either.
    with pytest.raises(CofactorTableError, match="g1"):
        _read(_table_text("g1\tA\tyes\tanimaltfdb", "g1\tB\tyes\tepifactors"))


def test_a_table_with_a_header_and_no_genes_raises() -> None:
    # Absence is spelled by shipping no file. A file listing nothing would be a second
    # spelling of it, and the one that reads as *this species has no cofactors*.
    with pytest.raises(CofactorTableError, match="no genes"):
        _read(_table_text())


def test_an_empty_table_raises() -> None:
    with pytest.raises(CofactorTableError, match="empty"):
        _read("")


def test_every_refusal_names_the_file_and_the_repair() -> None:
    # A packaging defect and not a caller error, so the message has to say which file and
    # what regenerates it — nothing else is available to whoever hits one.
    with pytest.raises(CofactorTableError) as raised:
        _read(_table_text("g1\tA\tyes\tnobody"))

    assert f"tiny_beast{COFACTOR_SUFFIX}" in str(raised.value)
    assert "scripts/build_tf_cofactor.py" in str(raised.value)


def test_a_broken_table_is_a_bad_value_and_not_a_lookup() -> None:
    # A caller catching LookupError for the absences above it must not swallow a
    # packaging defect.
    assert issubclass(CofactorTableError, ValueError)
    assert not issubclass(CofactorTableError, LookupError)


# ---------------------------------------------------------------------------------------
# The reader, over generated tables
# ---------------------------------------------------------------------------------------


@given(
    st.lists(
        st.tuples(_stem, _cell, st.sampled_from(SOURCES), _cell),
        min_size=1,
        max_size=8,
        unique_by=lambda row: row[0],
    )
)
def test_a_well_formed_table_round_trips_to_its_rows(
    rows: list[tuple[str, str, str, str]],
) -> None:
    # A blank cell reads back as None and nothing else changes: the round trip is the
    # whole contract between the generator and the reader.
    header = "\t".join([*UNIFORM_COLUMNS, _FAMILY])
    text = _table_text(
        *(
            f"{stem}\t{symbol}\t{TRUE_CELL}\t{source}\t{family}"
            for stem, symbol, source, family in rows
        ),
        header=header,
    )

    table = parse_cofactor_table(text, provenance=_provenance())

    assert table.columns == (*UNIFORM_COLUMNS, _FAMILY)
    assert table.gene_id_stems == tuple(stem for stem, _, _, _ in rows)
    assert table.rows == tuple(
        (stem, symbol or None, TRUE_CELL, source, family or None)
        for stem, symbol, source, family in rows
    )


@pytest.mark.parametrize("flag", [TRUE_CELL, FALSE_CELL])
def test_only_the_rows_flagged_a_cofactor_come_back_as_cofactor_stems(flag: str) -> None:
    table = _read(_table_text(f"g1\tA\t{flag}\t{ANIMALTFDB}"))

    assert table.cofactor_stems == (("g1",) if flag == TRUE_CELL else ())


# ---------------------------------------------------------------------------------------
# The two provenance tables beside the data
# ---------------------------------------------------------------------------------------


def test_the_provenance_tables_read_back_as_records() -> None:
    records = cofactor_metadata()

    assert records
    assert all(isinstance(record, CofactorProvenance) for record in records)
    assert all(
        isinstance(source, CofactorSource) for record in records for source in record.sources
    )


def test_the_provenance_tables_and_the_shipped_files_name_the_same_species() -> None:
    # A row with no file cannot be cited from anywhere; a file with no row cannot be
    # cited at all, which is the condition on redistributing a table here.
    assert sorted(species_slug(record.species) for record in cofactor_metadata()) == sorted(
        cofactor_species()
    )


def test_the_source_table_is_ragged_and_gives_every_species_at_least_one_publisher() -> None:
    # Two tables and not one: a species built from three publishers gets three rows here,
    # where one row joining them positionally inside a cell is the shape that breaks
    # quietly. Every species has at least one, since a table nobody can cite may not ship.
    for record in cofactor_metadata():
        assert record.sources
        assert len({source.source for source in record.sources}) == len(record.sources)
        assert all(source.species == record.species for source in record.sources)


def test_neither_provenance_table_carries_a_licence_or_a_family_column() -> None:
    # No licence column and no licence prose, matching the census metadata table; and no
    # family-column field either, since the shipped column names are already namespaced
    # and so already say whose vocabulary they are.
    for name in ("cofactor_metadata.tsv", "cofactor_source_metadata.tsv"):
        header = (
            files("genome")
            .joinpath(COFACTOR_SUBDIR, name)
            .read_text(encoding="utf-8")
            .splitlines()[0]
        )

        assert "licen" not in header.lower(), f"{name}: {header}"
        assert "family_column" not in header, f"{name}: {header}"


# ---------------------------------------------------------------------------------------
# The values themselves
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize("slug", cofactor_species())
def test_the_frame_is_the_table_with_the_flag_read_as_a_boolean(slug: str) -> None:
    table = _shipped(slug)
    frame = table.frame()

    assert list(frame.columns) == list(table.columns)
    assert len(frame) == len(table)
    assert frame["is_cofactor"].dtype == bool
    assert int(frame["is_cofactor"].sum()) == len(table.cofactor_stems)


def test_the_frame_is_built_fresh_so_mutating_it_cannot_reach_the_table() -> None:
    table = _shipped(cofactor_species()[0])
    frame = table.frame()
    frame.loc[0, "symbol"] = "MUTATED"

    assert table.frame().loc[0, "symbol"] != "MUTATED"


def test_a_cofactor_table_is_frozen() -> None:
    table = _read(_table_text("g1\tA\tyes\tanimaltfdb"))

    with pytest.raises(AttributeError):
        table.species = "Other beast"  # type: ignore[misc]


def test_the_flag_spellings_and_the_source_vocabulary_are_the_ones_a_table_uses() -> None:
    assert {TRUE_CELL, FALSE_CELL} == {"yes", "no"}
    assert SOURCES == ("animaltfdb", "epifactors", "both")


def test_the_two_source_vocabularies_differ_in_both_directions_on_purpose() -> None:
    # A row of a table says which publisher listed the gene; a provenance row says which
    # source a table owes a citation to. `both` belongs only to the first — it describes
    # a row and no publisher — and `hgnc` only to the second, since HGNC lists no gene
    # and only makes one readable. One vocabulary would have to admit both mistakes.
    assert CITED_SOURCES == ("animaltfdb", "epifactors", "hgnc")
    assert BOTH in SOURCES
    assert BOTH not in CITED_SOURCES
    assert HGNC in CITED_SOURCES
    assert HGNC not in SOURCES


def test_a_provenance_row_describing_no_publisher_raises() -> None:
    # `both` is the one value the two vocabularies disagree about most sharply: a
    # provenance row spelled that way names nobody to cite.
    row = {
        "species": "Tiny beast",
        "source": BOTH,
        "publisher": "Someone et al. 1999",
        "version": "v1",
        "pubmed_id": "2",
        "source_url": "https://example.org/beast",
    }

    with pytest.raises(CofactorTableError, match="a provenance row is spelled from"):
        CofactorSource.from_row(row, origin="cofactor_source_metadata.tsv")


def test_a_source_that_lists_no_gene_can_still_be_cited() -> None:
    # The stems of 442 human genes exist only because HGNC said so, so an identifier
    # crosswalk is a source with a row of its own rather than an implementation detail.
    row = {
        "species": "Tiny beast",
        "source": HGNC,
        "publisher": "HGNC",
        "version": "2026-08-07",
        "pubmed_id": "2",
        "source_url": "https://example.org/archive",
    }

    assert CofactorSource.from_row(row, origin="cofactor_source_metadata.tsv").source == HGNC


# ---------------------------------------------------------------------------------------
# Human: the union this package publishes, and nobody else's verdict
# ---------------------------------------------------------------------------------------


def test_the_human_table_carries_both_publishers_namespaced_columns() -> None:
    # A table built from two publishers is more columns and one more provenance row,
    # never a change of format: the uniform four still lead, and each publisher's group
    # follows under its own namespace with nothing crosswalked between them (ADR-0014).
    table = _shipped("homo_sapiens")

    assert table.columns == (
        *UNIFORM_COLUMNS,
        _FAMILY,
        _CATEGORY,
        _FUNCTION,
        _TARGET,
        _MODIFICATION,
        _COMPLEX,
    )


def test_membership_is_unioned_and_classification_is_not() -> None:
    # A gene either publisher lists is a row; its AnimalTFDB columns are filled only if
    # AnimalTFDB listed it and its EpiFactors columns only if EpiFactors did. `source`
    # reads `both` exactly when both groups say something, so a caller can tell an
    # agreement about membership from a blank one publisher never filled.
    for stem, row in _by_stem("homo_sapiens").items():
        animaltfdb = _filled(row, (_FAMILY, _CATEGORY))
        epifactors = _filled(row, _EPIFACTORS_COLUMNS)

        assert (row["source"] == BOTH) == (animaltfdb and epifactors), stem
        assert animaltfdb == (row["source"] in (ANIMALTFDB, BOTH)), stem
        # Not an equivalence in this direction, and deliberately not asserted as one:
        # HSFX3 is a gene EpiFactors lists and records nothing whatever about.
        if epifactors:
            assert row["source"] in (EPIFACTORS, BOTH), stem


def test_the_epifactors_vocabularies_are_pinned_as_atomic_values() -> None:
    # Counted after the split, which is the count a caller grouping by one of them sees.
    for column, pinned in _EPIFACTORS_VOCABULARIES.items():
        assert len(_atomic("homo_sapiens", column)) == pinned, column


def test_a_multi_valued_cell_splits_on_the_separator_the_whole_package_uses() -> None:
    # The same character `interpro_ids` and a Motif link's partners already use, so a
    # caller never has to remember which column uses which. A value that arrived with one
    # inside it would split into two here, which is what the build refuses to write.
    complexes = _atomic("homo_sapiens", _COMPLEX)

    assert VALUE_SEPARATOR == ";"
    assert "COMPASS-like MLL1,2" in complexes
    assert "COMPASS-like MLL3,4" in complexes
    assert not any(value.strip() != value or not value for value in complexes)


def test_the_five_double_rowed_genes_ship_as_one_row_with_their_cells_unioned() -> None:
    # EpiFactors gives each of these two rows, and they are not duplicates: one carries a
    # histone-modification annotation and the other an RNA-modification one. One row per
    # stem is required, so the cells are unioned and deduplicated within a cell as well
    # as across the two rows — the cost being that for these five the pairing between a
    # function and its own modification is gone, which ATTRIBUTION.md states.
    rows = _by_stem("homo_sapiens")

    for stem, symbol in _DOUBLE_ROWED.items():
        row = rows[stem]
        assert row["symbol"] == symbol
        for column in _EPIFACTORS_COLUMNS:
            values = (row[column] or "").split(VALUE_SEPARATOR)
            assert len(values) == len(set(values)), f"{symbol}: {column} repeats a value"

    alkbh1 = rows["ENSG00000100601"]
    assert (alkbh1[_FUNCTION] or "").split(VALUE_SEPARATOR) == [
        "Histone modification",
        "RNA modification",
        "DNA modification",
    ]
    assert (alkbh1[_TARGET] or "").split(VALUE_SEPARATOR) == ["histone", "DNA", "RNA"]


def test_the_human_symbol_is_hgncs_spelling_and_not_the_publishers_retired_one() -> None:
    # Why the build joins EpiFactors on its HGNC id and never on its symbol, made
    # executable: 31 of the publisher's 801 rows name a gene by a symbol HGNC has since
    # retired, so matching on the name would key those genes wrongly or drop them. These
    # three ship under the name HGNC approves, which a symbol match could not have found.
    rows = _by_stem("homo_sapiens")

    for stem, symbol in _RENAMED_BY_HGNC.items():
        assert rows[stem]["symbol"] == symbol


def test_the_human_table_cites_three_sources_including_the_one_that_lists_no_gene() -> None:
    # Two publishers of membership and one of identifiers. HGNC gets a row of its own
    # because the stems of 442 genes exist only because it said so.
    sources = {source.source: source for source in _shipped("homo_sapiens").provenance.sources}

    assert sorted(sources) == sorted(CITED_SOURCES)
    assert sources[ANIMALTFDB].publisher == "AnimalTFDB"
    assert sources[EPIFACTORS].publisher == "EpiFactors"
    assert sources[EPIFACTORS].version == "v2.0"
    # A dated monthly archive and never the rolling current file, so that the 442 stems
    # only HGNC can supply are reproducible. The date is the release identifier, and it
    # is the file the URL names rather than one constructed from a date.
    assert sources[HGNC].version == "2026-08-07"
    assert sources[HGNC].version in sources[HGNC].source_url


def test_the_genes_that_are_both_a_tf_and_a_cofactor_are_pinned() -> None:
    # ADR-0016's central number, and the one cost of publishing a union that a caller
    # feels: the TF gene list and the cofactor list overlap, so unioning the two answers
    # double-counts these genes. Read through the census's own public API rather than off
    # the file, so that a change to either shipped table is caught here.
    from genome.tf.gene import tf_gene_table

    census = tf_gene_table("Homo sapiens")
    assert census is not None
    positives = set(census.assessed_positive)
    rows = _by_stem("homo_sapiens")
    animaltfdb = {stem for stem, row in rows.items() if row["source"] in (ANIMALTFDB, BOTH)}
    epifactors = {stem for stem, row in rows.items() if row["source"] in (EPIFACTORS, BOTH)}

    assert len(positives & set(rows)) == _DUAL_CLASSIFIED
    assert len(positives & animaltfdb) == _DUAL_ANIMALTFDB
    assert len(positives & epifactors) == _DUAL_EPIFACTORS
    for stem, symbol in _DUAL_GENES.items():
        assert stem in positives, symbol
        assert stem in rows, symbol
        assert rows[stem]["symbol"] == symbol


def test_being_a_cofactor_never_removes_a_gene_from_the_census() -> None:
    # The two answers are independent verdicts about different questions, so a second
    # shipped table never suppresses one the census already reached.
    from genome.tf.gene import tf_gene_table

    census = tf_gene_table("Homo sapiens")
    assert census is not None

    assert "ENSG00000118058" in set(census.assessed_positive)
    assert "ENSG00000118058" in set(_shipped("homo_sapiens").cofactor_stems)
