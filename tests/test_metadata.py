"""Tests for genome.metadata — the curated assembly cross-reference table.

The table ships inside the package, so these run against the real TSV: no
fixture stands in for it and nothing here touches the network.
"""

from __future__ import annotations

from genome.metadata import METADATA_FIELDS, AssemblyMetadata, lookup_assembly


def test_lookup_returns_the_row_for_a_listed_assembly() -> None:
    assert lookup_assembly("hg38") == AssemblyMetadata(
        assembly_name="hg38",
        species="Homo sapiens",
        ucsc_name="hg38",
        ncbi_name="GRCh38",
        ncbi_assembly_id="GCF_000001405.40",
        ncbi_taxid=9606,
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
    )
    assert all(getattr(record, field) is not None for field in METADATA_FIELDS)
