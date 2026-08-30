"""Tests for genome.annotation.metadata — the curated annotation table.

Keyed by assembly and **Registered name**, and shipped inside the package, so these run
against the real TSV. Which assemblies the table lists at all is the assembly half's
question and is answered in ``tests/assembly/test_metadata.py``, whose reading of the
shipped rows this imports rather than repeats.
"""

from __future__ import annotations

from dataclasses import asdict

import pytest

from genome.annotation.metadata import (
    ANNOTATION_FIELDS,
    AnnotationMetadata,
    MetadataRowError,
    annotation_table,
    list_annotation_metadata,
    lookup_annotation,
)

from ..assembly.test_metadata import _SHIPPED_CHIMERAS, _SHIPPED_DOWNLOADS


class TestAnnotationTable:
    """The shipped annotation table — keyed by assembly and registered name."""

    def test_shipped_rows_fill_every_column_pin_provenance_and_flag_one_default(self) -> None:
        assert lookup_annotation("sacCer3", "ensgene_v101") == AnnotationMetadata(
            assembly="sacCer3",
            name="ensgene_v101",
            provider="UCSC",
            version="ensGene.v101",
            url="https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz",
            sha256="d3f33fbf97deef26e2495f709f1c5bb2e2e1bf1ce71fb80758c2c9de42ad7026",
            default=True,
        )
        # A cross-reference, not an allow-list: an unlisted GTF registers by path.
        assert lookup_annotation("sacCer3", "no_such_annotation") is None
        assert lookup_annotation("no_such_assembly", "gencode_v50") is None
        # The key is the pair: hg38's GENCODE row is not hg19's.
        assert lookup_annotation("hg19", "gencode_v50") is None
        # And listing an assembly returns exactly what the table offers it.
        assert [record.name for record in list_annotation_metadata("ce11")] == ["wormbase_ws298"]
        assert list_annotation_metadata("no_such_assembly") == []

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

        for assembly in _SHIPPED_DOWNLOADS:
            listed = list_annotation_metadata(assembly)
            assert listed, f"{assembly} offers no annotation"
            for row in listed:
                assert row.url.startswith("https://")
                assert row.sha256 is not None
                assert len(row.sha256) == 64
            defaults = [r.name for r in listed if r.default]
            assert len(defaults) == 1

        # A merged annotation is derived from its components' own and fetched from
        # nowhere, so it has nothing a row could pin: no source, no digest, and a name
        # that would be computed from the flags in the rows beside it.
        for assembly in _SHIPPED_CHIMERAS:
            assert list_annotation_metadata(assembly) == []

        mm39_default = lookup_annotation("mm39", "gencode_vM39")
        assert mm39_default is not None
        assert mm39_default.default is True

    def test_a_blank_flag_is_false_a_bad_flag_or_a_missing_url_raises(self) -> None:
        # Every shipped row pins a digest and is its assembly's default, so the blank
        # cells are exercised against a row stood up for them.
        record = AnnotationMetadata.from_row(
            {
                "assembly": "tiny",
                "name": "unpinned",
                "provider": "Nobody",
                "version": "0",
                "url": "https://example.invalid/unpinned.gtf.gz",
                "sha256": None,
                "default": None,
            }
        )
        assert record.sha256 is None
        assert record.default is False

        with pytest.raises(MetadataRowError, match="is not a flag"):
            AnnotationMetadata.from_row(
                {
                    "assembly": "tiny",
                    "name": "typo",
                    "provider": "Nobody",
                    "version": "0",
                    "url": "https://example.invalid/typo.gtf.gz",
                    "sha256": None,
                    "default": "y",
                }
            )

        # An annotation row that says where nothing is fetched from is malformed, not an
        # annotation whose source is unknown: nothing could register from it.
        with pytest.raises(MetadataRowError, match="url"):
            AnnotationMetadata.from_row(
                {"assembly": "tiny", "name": "sourceless", "provider": "Nobody", "version": "0"}
            )

    def test_a_record_is_rebuilt_from_its_own_fields_including_hand_built_defaults(self) -> None:
        record = lookup_annotation("sacCer3", "ensgene_v101")
        assert record is not None
        assert AnnotationMetadata.from_row(asdict(record)) == record

        by_hand = AnnotationMetadata(
            assembly="tiny",
            name="by_hand",
            provider="Nobody",
            version="0",
            url="https://example.invalid/by_hand.gtf.gz",
        )
        assert (by_hand.sha256, by_hand.default) == (None, False)

    def test_the_lookups_read_the_table_they_are_handed(self) -> None:
        offered = [
            AnnotationMetadata(
                assembly="tiny",
                name="mine",
                provider="Nobody",
                version="0",
                url="https://example.invalid/mine.gtf.gz",
                default=True,
            ),
            AnnotationMetadata(
                assembly="other",
                name="theirs",
                provider="Nobody",
                version="0",
                url="https://example.invalid/theirs.gtf.gz",
            ),
        ]

        # The shipped table is what a lookup reads when no table is handed.
        assert lookup_annotation("sacCer3", "ensgene_v101") in annotation_table()

        assert lookup_annotation("tiny", "mine", table=offered) == offered[0]
        assert lookup_annotation("tiny", "theirs", table=offered) is None
        assert list_annotation_metadata("tiny", table=offered) == [offered[0]]
        assert list_annotation_metadata("tiny", table=[]) == []
        # …and the shipped table is where it was.
        assert lookup_annotation("tiny", "mine") is None
        assert lookup_annotation("sacCer3", "ensgene_v101") is not None
