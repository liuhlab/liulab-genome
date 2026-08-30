"""Tests for genome.shipped_writer — the one writer the generators in ``scripts/`` go through.

Three build scripts used to carry the same renderer, the same gzip call and the same
provenance merge under the same names, and the promise that matters most travelled with
them: two runs write the same bytes. It is one implementation now, and this module is what
holds it to that promise **with no download at all** — the scripts themselves need a
publisher's file and CI has no network, so what can be checked here is the writing and not
the curation.

Nothing here writes into the package's own data directory: every case writes into ``tmp_path``
and the shipped files are only ever read. The last test is the one that crosses back — it
re-renders each shipped table's own rows and asserts the bytes are the bytes that ship, which
is what "the reader and the writer agree by construction" means when it is not merely stated.
"""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import pytest

from genome.shipped import ShippedTable, ShippedTableError
from genome.shipped_writer import WrittenTable, merge_rows, render, shipped_name, write_table
from genome.tf import LINK_FORMAT
from genome.tf.cofactor import (
    COFACTOR_FORMAT,
    COFACTOR_METADATA_FORMAT,
    COFACTOR_SOURCE_METADATA_FORMAT,
)
from genome.tf.gene import CENSUS_FORMAT, CENSUS_METADATA_FORMAT

#: The repository's own generators, which are what this module exists for. Read as text
#: rather than imported: ``scripts/`` is not a package, is not on the path a test runs
#: under, and is deliberately outside the wheel.
_SCRIPTS = sorted((Path(__file__).resolve().parents[1] / "scripts").glob("build_tf_*.py"))

#: A bulk table as one is declared: gzipped, one row per key, a flag column, a name with a
#: species in it. Small enough to write by hand and shaped like every table that ships.
_BULK = ShippedTable(
    resource="data/tiny/{slug}.tiny_table.tsv.gz",
    columns=("gene_id_stem", "symbol", "is_tiny"),
    noun="tiny table",
    repair="fix the recipe in scripts/build_tiny.py",
    unit="gene",
    key=("gene_id_stem",),
    required=("gene_id_stem",),
    flags=("is_tiny",),
)

#: A provenance table as one is declared: plain, keyed, and merged into rather than
#: rewritten, since a generator rebuilds one species and leaves the others alone.
_PROVENANCE = ShippedTable(
    resource="data/tiny/tiny_metadata.tsv",
    columns=("species", "file", "sha256"),
    noun="tiny provenance table",
    repair="fix the recipe in scripts/build_tiny.py",
    key=("species",),
)

#: Every shipped format and what names one of its files — the whole population, so a table
#: added without the writer being able to reproduce it fails here.
_SHIPPED = [
    ("census", CENSUS_FORMAT, {"slug": "homo_sapiens"}),
    ("census_metadata", CENSUS_METADATA_FORMAT, {}),
    ("cofactor", COFACTOR_FORMAT, {"slug": "mus_musculus"}),
    ("cofactor_metadata", COFACTOR_METADATA_FORMAT, {}),
    ("cofactor_source_metadata", COFACTOR_SOURCE_METADATA_FORMAT, {}),
    ("link", LINK_FORMAT, {"slug": "homo_sapiens", "release": "2026"}),
]

_ROWS = [("g1", "A", "yes"), ("g2", "B", "no")]


def _write(directory: Path, *rows: tuple[str, ...]) -> WrittenTable:
    """Write the tiny bulk table into ``directory``, defaulting to two ordinary rows."""
    return write_table(_BULK, directory, _BULK.columns, rows or _ROWS, slug="worm")


# ---------------------------------------------------------------------------------------
# Two runs write the same bytes
# ---------------------------------------------------------------------------------------


def test_two_renders_of_the_same_rows_are_one_byte_string() -> None:
    # The first half of the promise, and the only half that is a pure function: same
    # columns, same rows, same bytes — no dictionary order, no locale, no clock.
    assert render(_BULK.columns, _ROWS) == render(_BULK.columns, _ROWS)
    assert render(_BULK.columns, _ROWS) == b"gene_id_stem\tsymbol\tis_tiny\ng1\tA\tyes\ng2\tB\tno\n"


def test_two_writes_of_the_same_rows_are_one_file_byte_for_byte(tmp_path: Path) -> None:
    # The other half, which the compression is where it could be lost: gzip records the
    # time it ran unless it is told not to, and a rebuild would then diff every time.
    first = _write(tmp_path / "first")
    second = _write(tmp_path / "second")

    assert first.path.read_bytes() == second.path.read_bytes()
    assert first.unpacked == second.unpacked
    assert first.sha256 == second.sha256


def test_the_gzip_header_carries_a_zeroed_timestamp(tmp_path: Path) -> None:
    # Where the clock would have got in: bytes 4 to 8 of a gzip member are the modification
    # time, and this is the only place in the package that decides what goes there.
    raw = _write(tmp_path).path.read_bytes()

    assert raw.startswith(b"\x1f\x8b")
    assert raw[4:8] == b"\x00\x00\x00\x00"
    assert gzip.decompress(raw) == render(_BULK.columns, _ROWS)


def test_the_declaration_decides_the_name_and_whether_the_bytes_are_packed(
    tmp_path: Path,
) -> None:
    # Bulk gzipped and small metadata plain, decided by the name the declaration gives the
    # file — the same rule the reader unpacks by, so the two halves cannot pick differently.
    packed = _write(tmp_path)
    plain = merge_rows(_PROVENANCE, tmp_path, [{"species": "Worm", "file": "w", "sha256": "0"}])

    assert packed.path.name == shipped_name(_BULK, slug="worm") == "worm.tiny_table.tsv.gz"
    assert packed.path.read_bytes().startswith(b"\x1f\x8b")
    assert plain.path.name == "tiny_metadata.tsv"
    assert plain.path.read_bytes() == plain.unpacked
    assert plain.unpacked.endswith(b"\n")


def test_the_digest_is_over_the_unpacked_bytes(tmp_path: Path) -> None:
    # ADR-0006: what a provenance row pins is the TSV and never the gzip around it, so a
    # copy recompressed elsewhere still matches what the table says about itself.
    written = _write(tmp_path)

    assert written.sha256 == hashlib.sha256(written.unpacked).hexdigest()
    assert written.sha256 != hashlib.sha256(written.path.read_bytes()).hexdigest()
    assert written.packed == written.path.stat().st_size


# ---------------------------------------------------------------------------------------
# The merge: a row replaced by its key, and the file re-sorted
# ---------------------------------------------------------------------------------------


def test_a_merge_replaces_a_row_by_key_and_re_sorts_the_file(tmp_path: Path) -> None:
    # How a generator rebuilds one species and leaves the others alone. The key is the one
    # the declaration already names as the identity no two rows may repeat, so what makes a
    # row the same row is stated once and the merge cannot key on something else.
    rows = [
        {"species": "Mus musculus", "file": "mouse", "sha256": "1"},
        {"species": "Caenorhabditis elegans", "file": "worm", "sha256": "2"},
    ]
    merge_rows(_PROVENANCE, tmp_path, rows)

    merged = merge_rows(
        _PROVENANCE, tmp_path, [{"species": "Mus musculus", "file": "mouse", "sha256": "3"}]
    )

    assert merged.unpacked == (
        b"species\tfile\tsha256\nCaenorhabditis elegans\tworm\t2\nMus musculus\tmouse\t3\n"
    )
    # Sorted by key and not by arrival: worm was written second and reads first, and one
    # row per species, so the replaced row is gone rather than shadowed by a later one.
    assert merged.unpacked.count(b"Mus musculus") == 1


def test_a_merge_into_no_file_writes_the_one_row_and_its_header(tmp_path: Path) -> None:
    # The first run of a generator for a new species: there is nothing to merge into, which
    # is ordinary rather than an error, and the header comes from the declaration either way.
    merged = merge_rows(_PROVENANCE, tmp_path, [{"species": "Worm", "file": "w", "sha256": "0"}])

    assert merged.unpacked == b"species\tfile\tsha256\nWorm\tw\t0\n"


def test_a_merge_refuses_a_file_that_is_not_the_table_it_declares(tmp_path: Path) -> None:
    # Merging into a file whose header has moved would write one neither shape can read, so
    # the file already there is held to the declaration exactly as the reader holds it.
    (tmp_path / "tiny_metadata.tsv").write_text("species\tfile\nWorm\tw\n", encoding="utf-8")

    with pytest.raises(ShippedTableError) as raised:
        merge_rows(_PROVENANCE, tmp_path, [{"species": "Worm", "file": "w", "sha256": "0"}])

    message = str(raised.value)
    assert "tiny_metadata.tsv" in message
    assert _PROVENANCE.repair in message


def test_a_table_that_names_no_key_cannot_be_merged_into(tmp_path: Path) -> None:
    # A merge replaces a row *by* its key, so a table that declares none — a link table,
    # where one gene has many rows by design — says so rather than replacing arbitrarily.
    assert LINK_FORMAT.key == ()

    with pytest.raises(ValueError, match="write_table"):
        merge_rows(LINK_FORMAT, tmp_path, [{"release": "2026"}], slug="worm", release="2026")


def test_a_merged_row_missing_a_column_is_refused_rather_than_blanked(tmp_path: Path) -> None:
    # Every column is written as it stands, and a blank cell is a value rather than a column
    # nobody supplied — so a row that leaves one out is a defect in the generator.
    with pytest.raises(ValueError, match="sha256"):
        merge_rows(_PROVENANCE, tmp_path, [{"species": "Worm", "file": "w"}])


# ---------------------------------------------------------------------------------------
# What the reader would refuse never reaches disk
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("character", "named"), [("\t", "a tab"), ("\n", "a newline"), ("\r", "a carriage return")]
)
def test_a_cell_that_a_plain_tsv_cannot_carry_is_refused(
    tmp_path: Path, character: str, named: str
) -> None:
    # The file carries no quoting and is read by splitting on the tab and the newline, so
    # such a cell is a defect in what built it. The message names the file, the column, the
    # row and the repair, because fixing the recipe is the only thing anyone can do.
    with pytest.raises(ShippedTableError) as raised:
        _write(tmp_path, ("g1", f"A{character}B", "yes"))

    message = str(raised.value)
    assert named in message
    assert "'symbol'" in message
    assert "gene 1" in message
    assert _BULK.repair in message
    assert not list(tmp_path.iterdir())


def test_a_cell_that_is_not_text_is_refused(tmp_path: Path) -> None:
    # Every cell is written as it stands, so a value that has not been rendered to text is
    # named rather than reaching `join` as a TypeError with no file in it.
    with pytest.raises(ShippedTableError, match="not text"):
        write_table(_BULK, tmp_path, _BULK.columns, [("g1", 2, "yes")], slug="worm")  # type: ignore[list-item]


@pytest.mark.parametrize(
    ("what", "rows"),
    [
        ("a repeated key", [("g1", "A", "yes"), ("g1", "B", "yes")]),
        ("a blank required cell", [("", "A", "yes")]),
        ("a flag nobody spells that way", [("g1", "A", "true")]),
        ("a row of the wrong width", [("g1", "A")]),
    ],
    ids=["repeated-key", "blank-required", "misspelled-flag", "wrong-width"],
)
def test_what_the_reader_would_refuse_is_refused_before_it_reaches_disk(
    tmp_path: Path, what: str, rows: list[tuple[str, ...]]
) -> None:
    # The point of handing the writer the reader's own declaration: the file is held to
    # every rule it will be read under *before it exists*, so a generator cannot ship a
    # table that the package then refuses to load. None of these rules is written twice.
    with pytest.raises(ShippedTableError) as raised:
        _write(tmp_path, *rows)

    assert _BULK.repair in str(raised.value), what
    assert not list(tmp_path.iterdir()), what


def test_a_refusal_names_the_recipe_and_never_tells_the_builder_to_re_run_the_build(
    tmp_path: Path,
) -> None:
    # A generator replaces the declaration's repair with its own, because *re-run this* —
    # which is what a shipped file's own refusals say — is the one thing that cannot help
    # somebody whose run has just failed.
    from dataclasses import replace

    building = replace(_BULK, repair="fix the recipe in scripts/build_tiny.py")

    with pytest.raises(ShippedTableError) as raised:
        write_table(
            building,
            tmp_path,
            building.columns,
            [("g1", "A", "no"), ("g1", "B", "no")],
            slug="worm",
        )

    assert "fix the recipe in scripts/build_tiny.py" in str(raised.value)


# ---------------------------------------------------------------------------------------
# The reader and the writer agree by construction
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(("name", "table", "keys"), _SHIPPED, ids=[row[0] for row in _SHIPPED])
def test_every_shipped_table_re_renders_to_the_bytes_it_ships(
    name: str, table: ShippedTable, keys: dict[str, str]
) -> None:
    # The crossing: read one shipped file through the declaration, hand its own rows back to
    # the writer, and the bytes are the bytes that ship. The compression is left out of it
    # deliberately — a digest is pinned over the unpacked TSV (ADR-0006) precisely because
    # another zlib may pack the same bytes differently, and this is the same distinction.
    read = table.read(**keys)
    rows = [tuple(cell or "" for cell in row) for row in read.rows]

    assert render(read.columns, rows) == table.text(**keys).encode("utf-8"), name


def test_no_build_script_writes_a_table_of_its_own() -> None:
    # The duplication this module ended, guarded rather than remembered: three scripts each
    # declared a renderer, a gzip call and a provenance merge. A fourth generator inherits
    # them — and inherits the promise about the timestamp, which it would otherwise learn
    # only from a reader who happened to notice.
    assert len(_SCRIPTS) == 3

    for script in _SCRIPTS:
        text = script.read_text(encoding="utf-8")
        assert "from genome.shipped_writer import" in text, script.name
        for spelled in ("gzip.compress", "compresslevel", "mtime=0", "def render("):
            assert spelled not in text, f"{script.name} spells {spelled!r} again"
