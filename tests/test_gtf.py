"""Tests for genome.io.gtf — registering GTF annotations and building gffutils DBs.

These need only ``gffutils`` (a pure-Python/SQLite default dependency), not the
native bioinformatics binaries, so the module is not gated on a tool skip.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from genome.io.gtf import annotation_dir, list_annotations, register_gtf

# A minimal but valid GTF: one gene with a transcript and an exon. Standard
# gene/transcript features are declared, so the default no-inference path applies.
_GTF = (
    "\n".join(
        [
            'chrI\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
            'chrI\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
            'chrI\ttest\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
        ]
    )
    + "\n"
)


def test_register_plain_gtf_copies_and_builds_db(tmp_path: Path) -> None:
    src = tmp_path / "ann.gtf"
    src.write_text(_GTF)
    assembly = tmp_path / "asm"

    ann = register_gtf(assembly, src, "WS298")

    assert ann.gtf == annotation_dir(assembly, "WS298") / "WS298.gtf"
    assert ann.gtf.read_text() == _GTF
    assert ann.db.is_file()
    assert list(list_annotations(assembly)) == ["WS298"]


def test_register_gzipped_gtf_is_decompressed(tmp_path: Path) -> None:
    src = tmp_path / "ann.gtf.gz"
    with gzip.open(src, "wt") as fh:
        fh.write(_GTF)
    assembly = tmp_path / "asm"

    ann = register_gtf(assembly, src, "WS298")

    # Stored as a plain .gtf with decompressed contents, and the db builds.
    assert ann.gtf.suffix == ".gtf"
    assert ann.gtf.read_text() == _GTF
    assert ann.db.is_file()
    assert list(list_annotations(assembly)) == ["WS298"]


def test_register_missing_source_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="GTF file not found"):
        register_gtf(tmp_path / "asm", tmp_path / "nope.gtf", "X")


def test_reregister_warns_and_skips_without_force(tmp_path: Path) -> None:
    src = tmp_path / "ann.gtf"
    src.write_text(_GTF)
    assembly = tmp_path / "asm"
    first = register_gtf(assembly, src, "WS298")
    db_mtime = first.db.stat().st_mtime_ns

    with pytest.warns(UserWarning, match="already registered"):
        second = register_gtf(assembly, src, "WS298")

    # Returns the existing annotation without rebuilding the database.
    assert second == first
    assert second.db.stat().st_mtime_ns == db_mtime


def test_reregister_with_force_rebuilds(tmp_path: Path) -> None:
    src = tmp_path / "ann.gtf"
    src.write_text(_GTF)
    assembly = tmp_path / "asm"
    register_gtf(assembly, src, "WS298")

    ann = register_gtf(assembly, src, "WS298", force=True)

    assert ann.db.is_file()
    assert list(list_annotations(assembly)) == ["WS298"]
