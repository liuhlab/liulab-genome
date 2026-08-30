"""Tests for genome.shipped — the one reader every **Shipped table** in the wheel goes through.

Six modules used to bring their own loader along, and the six failures that recur across them
were checked in some of the copies and not others. This module is the guard on that not coming
back: the failures are parametrized over **every** declaration in the package, so a table added
without a declaration, or a check quietly dropped from one, fails here rather than passing
because nobody wrote the case again for the seventh table.

Nothing here touches the network and nothing writes into the package's own data directory. The
malformed cases are handed to a declaration as text, which is what :meth:`ShippedTable.parse` is
separate from the resource for; the shipped files themselves answer the read-back tests, since a
guard over what ships is the only check that can exist for data a download builds.
"""

from __future__ import annotations

import pytest

from genome.annotation import curated
from genome.annotation import metadata as annotation_module
from genome.assembly import metadata as assembly_module
from genome.assembly.metadata import MetadataRowError
from genome.homology import metadata as homology_module
from genome.shipped import FALSE_CELL, TRUE_CELL, ShippedTable, ShippedTableError, parse_cell
from genome.tf import link as link_module
from genome.tf.cofactor import table as cofactor_module
from genome.tf.gene import census as census_module
from genome.xref import metadata as xref_module

#: Every **Shipped table** this package reads, keyed by the name its own module gives it. The
#: whole population and not a sample: the point of one reader is that a rule reaches all of
#: them, which only a list of all of them can check. Adding a table means adding a line here,
#: and :func:`test_every_module_that_reads_a_shipped_table_declares_one` is what says so.
_TABLES: dict[str, ShippedTable] = {
    "assembly_metadata": assembly_module._ASSEMBLY_TABLE,
    "annotation_metadata": annotation_module._ANNOTATION_TABLE,
    "xref_metadata": xref_module._XREF_TABLE,
    "homology_metadata": homology_module._HOMOLOGY_METADATA,
    "census_metadata": census_module.CENSUS_METADATA_FORMAT,
    "census": census_module.CENSUS_FORMAT,
    "cofactor_metadata": cofactor_module.COFACTOR_METADATA_FORMAT,
    "cofactor_source_metadata": cofactor_module.COFACTOR_SOURCE_METADATA_FORMAT,
    "cofactor": cofactor_module.COFACTOR_FORMAT,
    "link": link_module.LINK_FORMAT,
}

#: The modules that own one, and the one module that reads shipped data and owns none.
_READERS = (
    assembly_module,
    annotation_module,
    xref_module,
    homology_module,
    census_module,
    cofactor_module,
    link_module,
)

#: What names one of the per-species or per-release tables, so a declaration whose resource is
#: a template can still be read. A plain path ignores them.
_KEYS = {"slug": "homo_sapiens", "release": "2026"}


def _cell(table: ShippedTable, column: str, number: int) -> str:
    """Return a cell that column will accept: a flag spelled ``yes``, a key spelled uniquely."""
    if column in table.flags:
        return TRUE_CELL
    if column in table.key:
        return f"key{number}"
    return "x"


def _text(table: ShippedTable, *, rows: int = 1, columns: tuple[str, ...] | None = None) -> str:
    """Return a table's text that nothing in ``table``'s declaration objects to."""
    header = table.columns if columns is None else columns
    lines = [
        "\t".join(header),
        *(
            "\t".join(_cell(table, column, number) for column in header)
            for number in range(1, rows + 1)
        ),
    ]
    return "\n".join(lines) + "\n"


#: The three tables whose header leads with uniform columns and whose rows are generated, so
#: they are the ones that declare a required column, a flag and — two of them — a key.
_GENERATED = ["census", "cofactor", "link"]

#: The four provenance tables, whose every column is required by its record's own field types.
_PROVENANCE = [
    "homology_metadata",
    "census_metadata",
    "cofactor_metadata",
    "cofactor_source_metadata",
]


# ---------------------------------------------------------------------------------------
# What ships: every declaration reads its own file back
# ---------------------------------------------------------------------------------------


def test_every_declared_table_reads_its_shipped_file_back_through_one_loader() -> None:
    # One loader, whether the bytes are gzipped or plain: the suffix decides, and no module
    # writes the decompression again. The three identical `gzip.decompress(...)` one-liners
    # this replaced are the reason the assertion crosses both kinds in one loop.
    plain, packed = 0, 0
    for name, table in _TABLES.items():
        read = table.read(**_KEYS)

        assert read.rows, f"{name}: the shipped file has no rows"
        if table.leading:
            assert read.columns[: len(table.columns)] == table.columns
        else:
            assert read.columns == table.columns
        if table.resource.endswith(".gz"):
            packed += 1
        else:
            plain += 1

    # Both halves of the convention are exercised — bulk gzipped, small metadata plain.
    assert plain
    assert packed


def test_every_module_that_reads_a_shipped_table_declares_one_and_gene_list_does_not() -> None:
    # The boundary the shared module's docstring draws, asserted rather than only stated: the
    # six that read a table declare one, and the curated gene lists — shipped JSON keyed by a
    # **Registered name**, nested rather than tabular — declare none, because they are not one.
    for module in _READERS:
        declared = [value for value in vars(module).values() if isinstance(value, ShippedTable)]
        assert declared, f"{module.__name__} reads shipped data and declares no table"
        assert all(table in _TABLES.values() for table in declared)

    assert not [value for value in vars(curated).values() if isinstance(value, ShippedTable)]


# ---------------------------------------------------------------------------------------
# The six failures, checked in one place and reaching every table that declares the rule
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize("name", list(_TABLES))
def test_an_empty_file_is_refused_by_every_table(name: str) -> None:
    table = _TABLES[name]

    with pytest.raises(table.error, match="empty") as raised:
        table.parse("", origin="broken.tsv")
    assert "broken.tsv" in str(raised.value)


@pytest.mark.parametrize("name", list(_TABLES))
def test_a_header_that_is_not_the_declared_columns_is_refused_by_every_table(name: str) -> None:
    # The gap this closes, not only the duplication it removes: the two pandas readers checked
    # no header at all, so a renamed or missing column in the assembly metadata table reached
    # the cell parser as a blank cell and was reported as one — or read as `None` where the
    # column was optional, which is the same silence with no message at all.
    table = _TABLES[name]
    renamed = ("renamed", *table.columns[1:])

    with pytest.raises(table.error) as raised:
        table.parse(_text(table, columns=renamed), origin="broken.tsv")

    message = str(raised.value)
    assert "broken.tsv" in message
    assert "columns" in message
    assert "blank" not in message


@pytest.mark.parametrize("name", list(_TABLES))
def test_a_row_with_the_wrong_cell_count_is_refused_by_every_table(name: str) -> None:
    table = _TABLES[name]
    short = _text(table).rsplit("\t", 1)[0] + "\n"

    with pytest.raises(table.error, match="line 2") as raised:
        table.parse(short, origin="broken.tsv")
    assert "cells" in str(raised.value)


@pytest.mark.parametrize("name", _GENERATED)
def test_a_blank_cell_in_a_required_column_is_refused(name: str) -> None:
    table = _TABLES[name]
    column = table.required[0]
    blanked = _text(table).replace(f"\n{_cell(table, column, 1)}", "\n", 1)

    with pytest.raises(table.error, match="line 2") as raised:
        table.parse(blanked, origin="broken.tsv")
    assert repr(column) in str(raised.value)


@pytest.mark.parametrize("name", _PROVENANCE)
def test_a_blank_cell_in_a_provenance_column_is_refused_with_that_tables_own_reason(
    name: str,
) -> None:
    # A provenance table's required columns are its record's own field types rather than a
    # list on the declaration, and the same refusal serves both: one message shape, one place.
    table = _TABLES[name]
    column = table.columns[0]

    with pytest.raises(table.error, match=repr(column)) as raised:
        table.record({column: ""}, {column: str}, origin="broken.tsv")

    message = str(raised.value)
    assert "broken.tsv" in message
    assert table.because in message


@pytest.mark.parametrize("name", _GENERATED)
def test_a_flag_spelled_a_way_no_table_spells_one_is_refused(name: str) -> None:
    table = _TABLES[name]
    column = table.flags[0]
    misspelled = _text(table).replace(f"\t{TRUE_CELL}", "\tTrue", 1)

    with pytest.raises(table.error, match="line 2") as raised:
        table.parse(misspelled, origin="broken.tsv")

    message = str(raised.value)
    assert repr(column) in message
    assert repr(TRUE_CELL) in message
    assert repr(FALSE_CELL) in message


@pytest.mark.parametrize("name", ["annotation_metadata", "xref_metadata"])
def test_a_curated_flag_column_takes_a_blank_as_no_and_refuses_a_spelling_nobody_uses(
    name: str,
) -> None:
    # The other flag rule, and the reason there are two: a curated table's flag column has a
    # real answer for a blank cell — *no* — where a generated table's does not.
    table = _TABLES[name]
    assert "default" in table.columns

    assert parse_cell("default", {}, {"default": bool}, table=table, origin="curated.tsv") is False
    with pytest.raises(table.error, match="is not a flag"):
        parse_cell("default", {"default": "y"}, {"default": bool}, table=table)


@pytest.mark.parametrize("name", ["census", "cofactor", "cofactor_source_metadata"])
def test_a_repeated_key_is_refused_where_a_table_declares_one(name: str) -> None:
    table = _TABLES[name]
    repeated = _text(table, rows=2).replace("key2", "key1")

    with pytest.raises(table.error, match="more than once") as raised:
        table.parse(repeated, origin="broken.tsv")
    assert "key1" in str(raised.value)


def test_a_link_table_declares_no_key_because_a_gene_has_many_motifs() -> None:
    # The refusal is declared per table, and this is the table that declines it: JASPAR has
    # several profiles for one factor, so a **Gene id stem** naming many rows is the format
    # working rather than a defect, and every shipped link table would be rejected by a rule
    # the other tables need. The shipped human table is what proves it is not hypothetical.
    assert _TABLES["link"].key == ()

    table = link_module.motif_link_table("Homo sapiens", "2026")
    assert table is not None
    stems = [link.gene_id_stem for link in table.links]
    assert len(set(stems)) < len(stems)


# ---------------------------------------------------------------------------------------
# The messages: each table's own noun, and its own repair
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize("name", list(_TABLES))
def test_every_refusal_names_the_file_the_table_is_called_and_what_repairs_it(name: str) -> None:
    # A shipped file that cannot be read is a packaging defect and not anything a caller did,
    # so naming the file and the repair is the whole of what a message can usefully carry.
    # Composed here from two halves the table declares, which is why they still read as prose.
    table = _TABLES[name]

    with pytest.raises(table.error) as raised:
        table.parse("", origin="broken.tsv")

    message = str(raised.value)
    assert "broken.tsv" in message
    assert table.noun in message
    # The repair opens this particular sentence, so its first letter is capitalized; the rest
    # of it is the declaration's own words either way.
    assert table.repair[1:] in message


def test_the_repair_each_table_names_is_pinned_word_for_word() -> None:
    # These read well because each names its own noun and its own repair, and the point of
    # composing them here is that they go on reading exactly as well. Pinned verbatim, since a
    # message that has become generic is the regression this refactor could have introduced.
    census = _raises(census_module.CENSUS_FORMAT, "")
    assert "Re-run scripts/build_tf_census.py" in census

    cofactor = _raises(cofactor_module.COFACTOR_FORMAT, "gene_id_stem\tsymbol\tis_cofactor\n")
    assert "re-run scripts/build_tf_cofactor.py for that species" in cofactor

    link = _raises(link_module.LINK_FORMAT, "nonsense\n")
    assert "re-run scripts/build_tf_links.py for that species and release" in link

    # Homology's is a reason rather than a command — the table is maintained by hand — and it
    # is the half a shared refusal cannot compose, so it is declared and pinned like the rest.
    blank = "\t".join(("116", "Homo sapiens", "Mus musculus", "Mus", "Ensembl", "1", "u", ""))
    with pytest.raises(homology_module.HomologyMetadataError) as raised:
        homology_module.read_metadata(
            "\t".join(homology_module.METADATA_COLUMNS) + f"\n{blank}\n", origin="broken.tsv"
        )
    assert (
        "Every provenance column is required: a set nobody can cite is one this package may "
        "not point anyone at" in str(raised.value)
    )


def _raises(table: ShippedTable, text: str) -> str:
    """Return the message ``table`` refuses ``text`` with, failing when it does not refuse."""
    with pytest.raises(table.error) as raised:
        table.parse(text, origin="broken.tsv")
    return str(raised.value)


# ---------------------------------------------------------------------------------------
# What moved, and the import paths that must still resolve
# ---------------------------------------------------------------------------------------


def test_the_names_that_moved_still_resolve_from_everywhere_they_used_to() -> None:
    # `species_slug` was imported from the curated-table module by nine modules and
    # re-exported by two more; `parse_cell` claimed in its own docstring to be the reader
    # every curated table shares. Both now live with the reader, and every path to them still
    # answers — with the same object, so a caller comparing them by identity is not surprised
    # either.
    from genome.assembly.metadata import parse_cell as from_metadata
    from genome.assembly.metadata import species_slug as slug_from_metadata
    from genome.shipped import species_slug as slug_from_shipped
    from genome.tf.cofactor import species_slug as slug_from_cofactor
    from genome.tf.gene import species_slug as slug_from_gene
    from genome.tf.gene.census import species_slug as slug_from_census

    assert (
        slug_from_metadata
        is slug_from_shipped
        is slug_from_gene
        is slug_from_cofactor
        is slug_from_census
    )
    assert from_metadata is parse_cell
    assert slug_from_metadata("Homo sapiens") == "homo_sapiens"


def test_every_tables_error_is_one_kind_of_defect_and_never_an_absence() -> None:
    # One base class, so a caller may catch a broken shipped file across all of them — and a
    # `ValueError` and never a `LookupError`, which is what absence is spelled with here, so
    # catching a species with no census cannot swallow a census that ships broken.
    for name, table in _TABLES.items():
        assert issubclass(table.error, ShippedTableError), name
        assert issubclass(table.error, ValueError), name
        assert not issubclass(table.error, LookupError), name

    assert issubclass(MetadataRowError, ShippedTableError)
