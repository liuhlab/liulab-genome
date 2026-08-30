"""Tests for genome.assembly.metadata — the curated assembly table.

The table ships inside the package, so these run against the real TSV: no fixture stands
in for it and nothing here touches the network. Its annotation half is the other
context's, and is tested in ``tests/annotation/test_metadata.py``.
"""

from __future__ import annotations

import io
import string
from dataclasses import asdict

import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.assembly.chimera import ChimeraNamingError, split_name
from genome.assembly.metadata import (
    METADATA_FIELDS,
    AssemblyMetadata,
    MetadataRowError,
    assembly_metadata,
    assembly_table,
    format_table_row,
    lookup_assembly,
)

#: Columns a shipped row is allowed to leave blank — the pinned source and digest fill
#: in over time, a reference UCSC has never carried has no name in that namespace, and
#: an assembly nobody has chosen an intron bound for carries neither cap column.
#: Every identifier column may be blank too; a curated row is expected to fill them.
_OPTIONAL_FIELDS = (
    "ucsc_name",
    "source_url",
    "sha256",
    "intron_length_cap",
    "intron_length_cap_rationale",
)

#: Every assembly the shipped table officially supports, read from the table itself: a
#: row added without a test to match is then covered by whichever kind it turns out to be
#: rather than quietly unguarded.
_SHIPPED_ASSEMBLIES: tuple[str, ...] = tuple(record.assembly_name for record in assembly_table())


def _is_chimera_row(assembly: str) -> bool:
    """Whether a shipped row names a chimera — what the row *is*, not which name it has.

    It splits into two or more parts and the table lists every one of them, which is the
    same test that separates ``ce11_ecHT115`` from a local key someone chose.
    """
    try:
        components = split_name(assembly)
    except ChimeraNamingError:
        return False
    return all(component in _SHIPPED_ASSEMBLIES for component in components)


#: The two kinds, because a blank ``sha256`` means opposite things in them: on a
#: downloaded row nobody has pinned it yet, and on a chimera's pinning one would be a
#: mistake — its bytes were derived here, not fetched from anywhere.
_SHIPPED_CHIMERAS: tuple[str, ...] = tuple(filter(_is_chimera_row, _SHIPPED_ASSEMBLIES))
_SHIPPED_DOWNLOADS: tuple[str, ...] = tuple(
    assembly for assembly in _SHIPPED_ASSEMBLIES if not _is_chimera_row(assembly)
)

#: Every intron bound the shipped table registers, and the whole of it: a value edited or
#: a row filled in goes red here rather than reaching a consumer's aligner unremarked.
#: Each is a hand-set round number with a recorded reason, never a maximum computed from
#: an annotation — whose longest intron is a floor on what the organism does rather than
#: a ceiling on it (ADR-0010).
_REGISTERED_INTRON_CAPS = {"ce11": 50_000, "ecHT115": 1, "hg38": 1_000_000, "mm39": 1_000_000}


#: A name a caller might ask for — mostly one nothing lists, sometimes one the shipped
#: table does, so the total accessor is generated onto both of its answers.
_ASSEMBLY_NAMES = st.sampled_from(_SHIPPED_ASSEMBLIES) | st.text(
    alphabet=string.ascii_letters + string.digits + "._-", min_size=1
)


def test_lookup_and_total_accessor_answer_listed_and_unlisted_assemblies() -> None:
    assert lookup_assembly("hg38") == AssemblyMetadata(
        assembly_name="hg38",
        species="Homo sapiens",
        ucsc_name="hg38",
        ncbi_name="GRCh38",
        ncbi_assembly_id="GCF_000001405.40",
        ncbi_taxid=9606,
        source_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        sha256="5be01555d98347fdb3714dc84c6f77c9d8bc774adcf32c6f7a8fa06f5baf5e51",
        intron_length_cap=1_000_000,
        intron_length_cap_rationale="ENCODE's value, a convention rather than a measurement; it clips 8 annotated introns",
    )
    # The table is a cross-reference, not an allow-list: no row is not an error.
    assert lookup_assembly("no_such_assembly") is None

    # *unknown* is a record whose fields are unknown, not the absence of one, so nothing
    # downstream has to guard a missing row before reading a field.
    assert AssemblyMetadata.unknown("my_ref") == AssemblyMetadata(
        assembly_name="my_ref",
        species=None,
        ucsc_name=None,
        ncbi_name=None,
        ncbi_assembly_id=None,
        ncbi_taxid=None,
        source_url=None,
        sha256=None,
        intron_length_cap=None,
        intron_length_cap_rationale=None,
    )

    # The total accessor hands back the shipped row when the table lists one...
    assert assembly_metadata("hg38") == lookup_assembly("hg38")
    # ...and an unknown record — never an absence — when it does not.
    record = assembly_metadata("no_such_assembly")
    assert record == AssemblyMetadata.unknown("no_such_assembly")
    assert record.assembly_name == "no_such_assembly"
    assert all(getattr(record, field) is None for field in METADATA_FIELDS[1:])

    # The two questions are different and stay in different functions. *Is this listed?*
    # is what separates a chimera's derived name from a free-form local key on a machine
    # holding neither (ADR-0003), so it keeps its ``None``; making it total would read
    # every name as listed and resolve ``my_ref`` as a chimera of ``my`` and ``ref``.
    assert lookup_assembly("my_ref") is None
    assert assembly_metadata("my_ref") == AssemblyMetadata.unknown("my_ref")
    assert lookup_assembly("ce11") is not None


@given(name=_ASSEMBLY_NAMES)
def test_the_total_accessor_answers_every_name(name: str) -> None:
    listed = lookup_assembly(name)
    assert assembly_metadata(name) == (
        listed if listed is not None else AssemblyMetadata.unknown(name)
    )


# --- the table is a parameter, and the shipped one is only its default -------


#: A row no shipped table carries, to hand to a lookup in place of the shipped one.
_TINY = AssemblyMetadata(
    assembly_name="tiny",
    species="Testus minimus",
    ucsc_name="tiny",
    ncbi_name="TINY.1",
    ncbi_assembly_id="GCF_0.0",
    ncbi_taxid=1,
)


def test_a_table_parameter_overrides_the_shipped_default_without_mutating_it() -> None:
    # The table is a cross-reference, not an allow-list (ADR-0003), so a caller curating
    # their own rows is the ordinary case rather than a test's trick.
    assert lookup_assembly("tiny", table=[_TINY]) == _TINY
    assert lookup_assembly("hg38", table=[_TINY]) is None
    assert assembly_metadata("tiny", table=[_TINY]) == _TINY
    assert assembly_metadata("hg38", table=[_TINY]) == AssemblyMetadata.unknown("hg38")

    # An empty table lists nothing, and that is not an error either.
    assert lookup_assembly("hg38", table=[]) is None
    assert assembly_metadata("hg38", table=[]) == AssemblyMetadata.unknown("hg38")

    # The shipped table is what a lookup reads when no table is handed...
    assert lookup_assembly("hg38") in assembly_table()

    # ...and nothing is installed or cached: handing one call a table leaves it alone
    # for the next.
    assert lookup_assembly("tiny", table=[_TINY]) is not None
    assert lookup_assembly("tiny") is None
    assert lookup_assembly("hg38") is not None


def test_shipped_rows_parse_types_pin_provenance_and_bound_introns_where_backed() -> None:
    sac_cer3 = lookup_assembly("sacCer3")
    assert sac_cer3 is not None
    assert type(sac_cer3.ncbi_taxid) is int
    assert sac_cer3.ncbi_taxid == 559292
    assert sac_cer3.sha256 == "6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3"

    record = lookup_assembly("mm39")
    assert record is not None
    assert METADATA_FIELDS == (
        "assembly_name",
        "species",
        "ucsc_name",
        "ncbi_name",
        "ncbi_assembly_id",
        "ncbi_taxid",
        "source_url",
        "sha256",
        "intron_length_cap",
        "intron_length_cap_rationale",
    )
    identifiers = [f for f in METADATA_FIELDS if f not in _OPTIONAL_FIELDS]
    assert all(getattr(record, field) is not None for field in identifiers)

    # E. coli HT115 is a real, officially supported reference with no UCSC name at all:
    # the assembly id is a local key and UCSC is only the default source.
    echt115 = lookup_assembly("ecHT115")
    assert echt115 is not None
    assert echt115.ucsc_name is None
    assert echt115.species == "Escherichia coli HT115"
    assert echt115.ncbi_assembly_id == "GCF_004354945.1"

    # The blank cell is absence, not a value — looking up "" or NaN finds nothing.
    assert lookup_assembly("") is None
    assert lookup_assembly("nan") is None

    # Every officially supported assembly whose bytes were fetched says where they come
    # from and what they hash to. The URL's file name is the source's business, not the
    # assembly's — ce11 is pinned to WormBase, whose names carry the bioproject and
    # release rather than the assembly, and ecHT115 to NCBI, which spells it .fna.gz.
    assert _SHIPPED_DOWNLOADS
    for assembly in _SHIPPED_DOWNLOADS:
        record = lookup_assembly(assembly)
        assert record is not None
        assert record.source_url is not None
        assert record.source_url.endswith(".gz")
        assert record.sha256 is not None
        assert len(record.sha256) == 64

    # The other meaning of a blank cell: not *nobody got round to it* but *pinning this
    # would be a mistake*. A chimera's bytes are derived by a pure function from
    # components that are themselves pinned, so it is proven transitively; a digest here
    # would turn this package's own concatenation into a contract that fails on every
    # user's disk instead of in its test suite. Its identifiers are its components' too.
    assert _SHIPPED_CHIMERAS, "nothing is guarded if the table lists no chimera at all"
    for assembly in _SHIPPED_CHIMERAS:
        record = lookup_assembly(assembly)
        assert record is not None
        assert record.source_url is None
        assert record.sha256 is None
        assert record.species is None
        assert record.ncbi_taxid is None

    # Readable at all only because a blank identifier cell reads back as unknown: a blank
    # ncbi_taxid used to raise, which is what made this row impossible to ship. It exists
    # so that a machine holding neither component can still tell this name from a
    # free-form local key, by splitting it into components the table lists.
    assert lookup_assembly("ce11_ecHT115") == AssemblyMetadata(
        assembly_name="ce11_ecHT115",
        species=None,
        ucsc_name=None,
        ncbi_name=None,
        ncbi_assembly_id=None,
        ncbi_taxid=None,
    )
    assert split_name("ce11_ecHT115") == ("ce11", "ecHT115")

    # sacCer3 is the row to look at: nothing backs a number for yeast the way WormBase's
    # own pipeline backs the worm's, so its cell is blank deliberately rather than
    # pending. A chimera's is blank for a third reason again — its bound is the maximum
    # over its components, derived by whoever aligns against it from a component set this
    # row does not repeat. Blank reads back as ``None``, which is what leaves an assembly
    # nobody has characterised aligning exactly as it did before.
    caps = {
        assembly: assembly_metadata(assembly).intron_length_cap for assembly in _SHIPPED_ASSEMBLIES
    }

    assert caps == {
        assembly: _REGISTERED_INTRON_CAPS.get(assembly) for assembly in _SHIPPED_ASSEMBLIES
    }
    assert caps["sacCer3"] is None
    assert all(caps[assembly] is None for assembly in _SHIPPED_CHIMERAS)

    # The reason rides in the row beside the value rather than in a commit message, so a
    # later reader can tell a convention — ENCODE's mammalian million, run for a decade
    # and below the longest annotated intron in both mammals it was written for — from a
    # number read off an organism's own curators. A value with no reason beside it is the
    # failure this guards: it is unrevisitable, since nobody can tell what it would take
    # to change it.
    for assembly in _SHIPPED_ASSEMBLIES:
        record = assembly_metadata(assembly)
        carries_reason = record.intron_length_cap_rationale is not None
        assert carries_reason is (record.intron_length_cap is not None), assembly


def test_from_row_reads_blank_and_missing_optional_cells_as_none() -> None:
    # An unpinned checksum is unverified, not wrong — and never the string "nan". Every
    # downloaded row pins both optional columns today, and the chimera row's blanks say
    # something else again, so this path is exercised against a row stood up for it.
    record = AssemblyMetadata.from_row(
        {
            "assembly_name": "unpinned",
            "species": "Testus minimus",
            "ucsc_name": "unpinned",
            "ncbi_name": "TINY.1",
            "ncbi_assembly_id": "GCF_0.0",
            "ncbi_taxid": "1",
            "source_url": "https://example.invalid/unpinned.fa.gz",
            "sha256": None,
        }
    )
    assert record.source_url == "https://example.invalid/unpinned.fa.gz"
    assert record.sha256 is None

    # A hand-written row need not carry every column: what it does not say is unknown.
    assert AssemblyMetadata.from_row({"assembly_name": "sparse"}) == AssemblyMetadata.unknown(
        "sparse"
    )

    # And the same defaults apply when a record is built by hand rather than parsed.
    by_hand = AssemblyMetadata("tiny", "Testus minimus", "tiny", "TINY.1", "GCF_0.0", 1)
    assert (by_hand.source_url, by_hand.sha256) == (None, None)
    assert (by_hand.intron_length_cap, by_hand.intron_length_cap_rationale) == (None, None)

    # What `genome table-row` emits for an assembly the table does not list yet: the
    # name, the source and the digest, and blanks for everything only a person supplies.
    # Pasted into the table, it used to give species the string "nan" and raise on taxid.
    line = format_table_row(
        {
            "assembly_name": "newAsm",
            "source_url": "https://example.invalid/newAsm.fa.gz",
            "sha256": "0" * 64,
        }
    )
    emitted = lookup_assembly("newAsm", table=[AssemblyMetadata.from_row(_pasted(line))])
    assert emitted is not None
    assert emitted.species is None
    assert emitted.ucsc_name is None
    assert emitted.ncbi_name is None
    assert emitted.ncbi_assembly_id is None
    assert emitted.ncbi_taxid is None
    assert emitted.intron_length_cap is None
    assert emitted.sha256 == "0" * 64


def test_from_row_raises_naming_the_offending_cell_without_partial_state() -> None:
    with pytest.raises(MetadataRowError, match="ncbi_taxid") as raised:
        AssemblyMetadata.from_row({"assembly_name": "tiny", "ncbi_taxid": "many"})
    assert "int" in str(raised.value)
    assert "many" in str(raised.value)

    # The other half of what a blank cell means: an optional column reads back unknown,
    # and a column with no unknown is a malformed row rather than the string "nan".
    with pytest.raises(MetadataRowError, match="assembly_name") as blank:
        AssemblyMetadata.from_row({"assembly_name": None, "species": "Testus minimus"})
    assert "blank" in str(blank.value)

    # A row whose good columns come first still raises rather than handing back a record
    # carrying them: a caller gets a record or an error naming the cell, never both.
    row = {"assembly_name": "tiny", "species": "Testus minimus", "ncbi_taxid": "many"}
    with pytest.raises(MetadataRowError, match="ncbi_taxid") as malformed:
        AssemblyMetadata.from_row(row)
    assert "assembly_name" not in str(malformed.value)


def test_format_table_row_renders_columns_in_order_and_leaves_unknowns_blank() -> None:
    # sacCer3 pins a source and a digest and carries no intron bound, so the line it
    # renders to also says what an unfilled cell looks like in the middle of a full row.
    record = lookup_assembly("sacCer3")
    assert record is not None
    assert format_table_row(asdict(record)).split("\t") == [
        "sacCer3",
        "Saccharomyces cerevisiae",
        "sacCer3",
        "R64-1-1",
        "GCF_000146045.2",
        "559292",
        "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz",
        "6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3",
        "",
        "",
    ]

    row = format_table_row({"assembly_name": "newAsm", "sha256": "abc123"})
    assert row.split("\t") == ["newAsm", "", "", "", "", "", "", "abc123", "", ""]


# --- a rendered row reads back as itself -------------------------------------


#: Spellings the table's reader takes for a blank cell whatever a row meant by them, so
#: a cell holding one is blank by the reader's own rules rather than a counterexample.
_BLANK_SPELLINGS = frozenset(
    {"NA", "N/A", "n/a", "NULL", "null", "None", "NaN", "nan", "-NaN", "-nan"}
)

#: A filled cell: an identifier, a URL or a digest. Never empty, since an empty cell is
#: how the table spells unknown, and never carrying a tab, which would split the row.
_CELLS = st.text(alphabet=string.ascii_letters + string.digits + "._-:/", min_size=1).filter(
    lambda cell: cell not in _BLANK_SPELLINGS
)

#: A row with an arbitrary subset of the columns a row may leave blank left blank.
_RECORDS = st.builds(
    AssemblyMetadata,
    assembly_name=_CELLS,
    species=st.none() | _CELLS,
    ucsc_name=st.none() | _CELLS,
    ncbi_name=st.none() | _CELLS,
    ncbi_assembly_id=st.none() | _CELLS,
    ncbi_taxid=st.none() | st.integers(min_value=1),
    source_url=st.none() | _CELLS,
    sha256=st.none() | _CELLS,
    intron_length_cap=st.none() | st.integers(min_value=1),
    intron_length_cap_rationale=st.none() | _CELLS,
)


def _pasted(line: str) -> dict[str, object]:
    """Read a rendered line back the way the shipped table is read.

    Pasting the line into ``data/assembly_metadata.tsv`` is what the register-then-paste
    flow does, and reading it as text is what makes a blank cell arrive as the NaN pandas
    reads it as rather than as ``""`` — which is the whole of the bug. Only the reading
    is stood in for; the row is parsed by :meth:`AssemblyMetadata.from_row` itself.
    """
    text = "\t".join(METADATA_FIELDS) + "\n" + line + "\n"
    frame = pd.read_csv(io.StringIO(text), sep="\t", dtype=str)
    return dict(frame.iloc[0])


@given(record=_RECORDS)
def test_a_row_survives_rendering_and_parsing_unchanged(record: AssemblyMetadata) -> None:
    assert AssemblyMetadata.from_row(_pasted(format_table_row(asdict(record)))) == record


@given(record=_RECORDS)
def test_a_record_is_rebuilt_from_its_own_fields(record: AssemblyMetadata) -> None:
    # The other direction of the same seam: a record's fields are a row, so a caller
    # holding one can correct a cell and hand it straight back without rendering it.
    assert AssemblyMetadata.from_row(asdict(record)) == record


@given(name=_CELLS)
def test_an_unknown_record_survives_rendering_and_parsing_unchanged(name: str) -> None:
    # An unknown record is a row like any other, and this is the one that says so: it
    # renders to the line `genome table-row` emits for an assembly nobody has curated —
    # the name, and every other cell blank — and pasting that line back yields it again.
    record = AssemblyMetadata.unknown(name)
    assert AssemblyMetadata.from_row(_pasted(format_table_row(asdict(record)))) == record
