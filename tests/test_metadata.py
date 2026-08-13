"""Tests for genome.metadata — the curated assembly and annotation tables.

Both tables ship inside the package, so these run against the real TSVs: no
fixture stands in for them and nothing here touches the network.
"""

from __future__ import annotations

from dataclasses import asdict

import pandas as pd
import pytest

from genome import metadata
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

#: Columns of the assembly table that may be blank — the pinned source and digest fill
#: in over time, and a reference UCSC has never carried has no name in that namespace.
_OPTIONAL_FIELDS = ("ucsc_name", "source_url", "sha256")

#: Every assembly the shipped table officially supports.
_SHIPPED_ASSEMBLIES = ("hg38", "hg19", "mm39", "mm10", "sacCer3", "ce11", "ecHT115")


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


def test_every_shipped_row_pins_a_source_and_a_checksum() -> None:
    # Every officially supported assembly says where its FASTA comes from and what
    # that FASTA hashes to. The URL's file name is the source's business, not the
    # assembly's — ce11 is pinned to WormBase, whose names carry the bioproject and
    # release rather than the assembly, and ecHT115 to NCBI, which spells it .fna.gz.
    for assembly in _SHIPPED_ASSEMBLIES:
        record = lookup_assembly(assembly)
        assert record is not None
        assert record.source_url is not None
        assert record.source_url.endswith(".gz")
        assert record.sha256 is not None
        assert len(record.sha256) == 64


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
    # An unpinned checksum is unverified, not wrong — and never the string "nan".
    # Every shipped row pins both optional columns today, so the blank-cell path is
    # exercised against a table stood up for it rather than by leaving a row unpinned.
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
        for assembly in _SHIPPED_ASSEMBLIES:
            listed = list_annotation_metadata(assembly)
            assert listed, f"{assembly} offers no annotation"
            for record in listed:
                assert record.url.startswith("https://")
                assert record.sha256 is not None
                assert len(record.sha256) == 64

    def test_every_shipped_assembly_names_exactly_one_default(self) -> None:
        for assembly in _SHIPPED_ASSEMBLIES:
            defaults = [r.name for r in list_annotation_metadata(assembly) if r.default]
            assert len(defaults) == 1

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
