"""Tests for genome.metadata — the curated assembly and annotation tables.

Both tables ship inside the package, so these run against the real TSVs: no
fixture stands in for them and nothing here touches the network.
"""

from __future__ import annotations

import io
import string
from dataclasses import asdict

import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome import metadata
from genome.chimera import ChimeraNamingError, split_name
from genome.metadata import (
    ANNOTATION_FIELDS,
    METADATA_FIELDS,
    AnnotationMetadata,
    AssemblyMetadata,
    format_table_row,
    list_annotation_metadata,
    lookup_annotation,
    lookup_assembly,
)

#: Columns a shipped row is allowed to leave blank — the pinned source and digest fill
#: in over time, and a reference UCSC has never carried has no name in that namespace.
#: Every identifier column may be blank too; a curated row is expected to fill them.
_OPTIONAL_FIELDS = ("ucsc_name", "source_url", "sha256")

#: Every assembly the shipped table officially supports, read from the table itself: a
#: row added without a test to match is then covered by whichever kind it turns out to be
#: rather than quietly unguarded.
_SHIPPED_ASSEMBLIES: tuple[str, ...] = tuple(metadata._metadata_table()["assembly_name"])


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


def test_lookup_returns_the_row_for_a_listed_assembly() -> None:
    assert lookup_assembly("hg38") == AssemblyMetadata(
        assembly_name="hg38",
        species="Homo sapiens",
        ucsc_name="hg38",
        ncbi_name="GRCh38",
        ncbi_assembly_id="GCF_000001405.40",
        ncbi_taxid=9606,
        source_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        sha256="5be01555d98347fdb3714dc84c6f77c9d8bc774adcf32c6f7a8fa06f5baf5e51",
    )


def test_lookup_returns_none_for_an_unlisted_assembly() -> None:
    # The table is a cross-reference, not an allow-list: no row is not an error.
    assert lookup_assembly("no_such_assembly") is None


def test_taxid_is_parsed_as_a_python_int() -> None:
    record = lookup_assembly("sacCer3")
    assert record is not None
    assert type(record.ncbi_taxid) is int
    assert record.ncbi_taxid == 559292


def test_every_declared_field_is_filled_from_the_table() -> None:
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
    )
    identifiers = [f for f in METADATA_FIELDS if f not in _OPTIONAL_FIELDS]
    assert all(getattr(record, field) is not None for field in identifiers)


def test_every_shipped_downloaded_row_pins_a_source_and_a_checksum() -> None:
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


def test_a_shipped_chimera_pins_nothing() -> None:
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


def test_the_shipped_chimera_row_carries_a_name_and_nothing_else() -> None:
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


def test_an_assembly_ucsc_never_carried_is_looked_up_by_its_own_name() -> None:
    # E. coli HT115 is a real, officially supported reference with no UCSC name at
    # all: the assembly id is a local key and UCSC is only the default source.
    record = lookup_assembly("ecHT115")
    assert record is not None
    assert record.ucsc_name is None
    assert record.species == "Escherichia coli HT115"
    assert record.ncbi_assembly_id == "GCF_004354945.1"


def test_a_blank_ucsc_name_never_matches_a_lookup() -> None:
    # The blank cell is absence, not a value — looking up "" or NaN finds nothing.
    assert lookup_assembly("") is None
    assert lookup_assembly("nan") is None


def test_a_pinned_checksum_is_read_back_as_text() -> None:
    record = lookup_assembly("sacCer3")
    assert record is not None
    assert record.sha256 == "6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3"


def test_a_blank_optional_cell_reads_back_as_none(monkeypatch: pytest.MonkeyPatch) -> None:
    # An unpinned checksum is unverified, not wrong — and never the string "nan". Every
    # downloaded row pins both optional columns today, and the chimera row's blanks say
    # something else again, so this path is exercised against a table stood up for it.
    table = pd.DataFrame(
        [
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
        ],
        dtype=str,
    )
    monkeypatch.setattr(metadata, "_metadata_table", lambda: table)
    metadata.lookup_assembly.cache_clear()
    try:
        record = metadata.lookup_assembly("unpinned")
    finally:
        # The cache is module-global; leaving the stand-in's row in it would leak
        # into every later test.
        metadata.lookup_assembly.cache_clear()

    assert record is not None
    assert record.source_url == "https://example.invalid/unpinned.fa.gz"
    assert record.sha256 is None


def test_optional_fields_default_to_none_when_a_record_is_built_by_hand() -> None:
    record = AssemblyMetadata("tiny", "Testus minimus", "tiny", "TINY.1", "GCF_0.0", 1)
    assert (record.source_url, record.sha256) == (None, None)


def test_format_table_row_renders_the_columns_in_table_order() -> None:
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
    ]


def test_format_table_row_leaves_unknown_values_blank() -> None:
    row = format_table_row({"assembly_name": "newAsm", "sha256": "abc123"})
    assert row.split("\t") == ["newAsm", "", "", "", "", "", "", "abc123"]


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
)


def _pasted(line: str) -> pd.DataFrame:
    """Read a rendered line back the way the shipped table is read.

    Pasting the line into ``data/assembly_metadata.tsv`` is what the register-then-paste
    flow does, and reading it as text is what makes a blank cell arrive as the NaN pandas
    reads it as rather than as ``""`` — which is the whole of the bug.
    """
    text = "\t".join(METADATA_FIELDS) + "\n" + line + "\n"
    return pd.read_csv(io.StringIO(text), sep="\t", dtype=str)


def _parsed(table: pd.DataFrame) -> AssemblyMetadata:
    """Build a record from the first row of ``table``, column by declared column."""
    row = table.iloc[0]
    return AssemblyMetadata(
        **{name: metadata._parse_cell(name, row, metadata._FIELD_TYPES) for name in METADATA_FIELDS}
    )


@given(record=_RECORDS)
def test_a_row_survives_rendering_and_parsing_unchanged(record: AssemblyMetadata) -> None:
    assert _parsed(_pasted(format_table_row(asdict(record)))) == record


def test_a_row_with_every_identifier_blank_reads_back_as_unknown(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
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
    monkeypatch.setattr(metadata, "_metadata_table", lambda: _pasted(line))
    metadata.lookup_assembly.cache_clear()
    try:
        record = metadata.lookup_assembly("newAsm")
    finally:
        # The cache is module-global; leaving the stand-in's row in it would leak.
        metadata.lookup_assembly.cache_clear()

    assert record is not None
    assert record.species is None
    assert record.ucsc_name is None
    assert record.ncbi_name is None
    assert record.ncbi_assembly_id is None
    assert record.ncbi_taxid is None
    assert record.sha256 == "0" * 64


class TestAnnotationTable:
    """The shipped annotation table — keyed by assembly and registered name."""

    def test_lookup_returns_the_row_for_a_listed_annotation(self) -> None:
        assert lookup_annotation("sacCer3", "ensgene_v101") == AnnotationMetadata(
            assembly="sacCer3",
            name="ensgene_v101",
            provider="UCSC",
            version="ensGene.v101",
            url="https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz",
            sha256="d3f33fbf97deef26e2495f709f1c5bb2e2e1bf1ce71fb80758c2c9de42ad7026",
            default=True,
        )

    def test_lookup_returns_none_for_an_unlisted_annotation(self) -> None:
        # A cross-reference, not an allow-list: an unlisted GTF registers by path.
        assert lookup_annotation("sacCer3", "no_such_annotation") is None
        assert lookup_annotation("no_such_assembly", "gencode_v50") is None

    def test_an_annotation_belongs_to_exactly_one_assembly(self) -> None:
        # The key is the pair: hg38's GENCODE row is not hg19's.
        assert lookup_annotation("hg19", "gencode_v50") is None

    def test_the_table_carries_every_declared_column(self) -> None:
        assert ANNOTATION_FIELDS == (
            "assembly",
            "name",
            "provider",
            "version",
            "url",
            "sha256",
            "default",
        )
        record = lookup_annotation("hg38", "gencode_v50")
        assert record is not None
        assert all(getattr(record, field) is not None for field in ANNOTATION_FIELDS)

    def test_every_shipped_annotation_pins_a_url_and_a_checksum(self) -> None:
        for assembly in _SHIPPED_DOWNLOADS:
            listed = list_annotation_metadata(assembly)
            assert listed, f"{assembly} offers no annotation"
            for record in listed:
                assert record.url.startswith("https://")
                assert record.sha256 is not None
                assert len(record.sha256) == 64

    def test_every_shipped_assembly_names_exactly_one_default(self) -> None:
        for assembly in _SHIPPED_DOWNLOADS:
            defaults = [r.name for r in list_annotation_metadata(assembly) if r.default]
            assert len(defaults) == 1

    def test_a_shipped_chimera_offers_no_annotation_row(self) -> None:
        # A merged annotation is derived from its components' own and fetched from
        # nowhere, so it has nothing a row could pin: no source, no digest, and a name
        # that would be computed from the flags in the rows beside it.
        for assembly in _SHIPPED_CHIMERAS:
            assert list_annotation_metadata(assembly) == []

    def test_listing_an_assembly_returns_what_the_table_offers_it(self) -> None:
        assert [record.name for record in list_annotation_metadata("ce11")] == ["wormbase_ws298"]
        assert list_annotation_metadata("no_such_assembly") == []

    def test_the_flag_column_is_parsed_as_a_bool(self) -> None:
        record = lookup_annotation("mm39", "gencode_vM39")
        assert record is not None
        assert record.default is True

    def test_a_blank_flag_reads_as_false_and_a_blank_digest_as_none(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # Every shipped row pins a digest and is its assembly's default, so the blank
        # cells are exercised against a table stood up for them.
        table = pd.DataFrame(
            [
                {
                    "assembly": "tiny",
                    "name": "unpinned",
                    "provider": "Nobody",
                    "version": "0",
                    "url": "https://example.invalid/unpinned.gtf.gz",
                    "sha256": None,
                    "default": None,
                }
            ],
            dtype=str,
        )
        monkeypatch.setattr(metadata, "_annotation_table", lambda: table)

        record = lookup_annotation("tiny", "unpinned")

        assert record is not None
        assert record.sha256 is None
        assert record.default is False

    def test_a_flag_nobody_spells_that_way_says_so(self, monkeypatch: pytest.MonkeyPatch) -> None:
        table = pd.DataFrame(
            [
                {
                    "assembly": "tiny",
                    "name": "typo",
                    "provider": "Nobody",
                    "version": "0",
                    "url": "https://example.invalid/typo.gtf.gz",
                    "sha256": None,
                    "default": "y",
                }
            ],
            dtype=str,
        )
        monkeypatch.setattr(metadata, "_annotation_table", lambda: table)

        with pytest.raises(ValueError, match="is not a flag"):
            lookup_annotation("tiny", "typo")

    def test_optional_fields_default_when_a_record_is_built_by_hand(self) -> None:
        record = AnnotationMetadata(
            assembly="tiny",
            name="by_hand",
            provider="Nobody",
            version="0",
            url="https://example.invalid/by_hand.gtf.gz",
        )
        assert (record.sha256, record.default) == (None, False)
