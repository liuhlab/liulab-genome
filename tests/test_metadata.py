"""Tests for genome.metadata — the curated assembly cross-reference table.

The table ships inside the package, so these run against the real TSV: no
fixture stands in for it and nothing here touches the network.
"""

from __future__ import annotations

from dataclasses import asdict

from genome.metadata import (
    METADATA_FIELDS,
    AssemblyMetadata,
    format_table_row,
    lookup_assembly,
)

#: Columns a row carries beyond its six identifiers, both of which may be blank.
_OPTIONAL_FIELDS = ("source_url", "sha256")


def test_lookup_returns_the_row_for_a_listed_assembly() -> None:
    assert lookup_assembly("hg38") == AssemblyMetadata(
        assembly_name="hg38",
        species="Homo sapiens",
        ucsc_name="hg38",
        ncbi_name="GRCh38",
        ncbi_assembly_id="GCF_000001405.40",
        ncbi_taxid=9606,
        source_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        sha256=None,
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


def test_every_shipped_row_pins_a_source_url() -> None:
    # The checksum column fills in over time, but every officially supported assembly
    # says where its FASTA comes from.
    for assembly in ("hg38", "hg19", "mm39", "mm10", "sacCer3", "ce11"):
        record = lookup_assembly(assembly)
        assert record is not None
        assert record.source_url is not None
        assert record.source_url.endswith(f"/{assembly}.fa.gz")


def test_a_pinned_checksum_is_read_back_as_text() -> None:
    record = lookup_assembly("sacCer3")
    assert record is not None
    assert record.sha256 == "6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3"


def test_a_blank_optional_cell_reads_back_as_none() -> None:
    # An unpinned checksum is unverified, not wrong — and never the string "nan".
    record = lookup_assembly("ce11")
    assert record is not None
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
