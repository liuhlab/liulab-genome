"""Tests for genome.io.fasta.

The end-to-end tests need the native binaries (samtools, faToTwoBit,
twoBitInfo); they run inside the pixi env and skip cleanly outside it. The
input-validation tests need no binaries and always run.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from genome.io.fasta import (
    GenomeFiles,
    faidx,
    fasta_to_2bit,
    prepare_fasta,
    twobit_to_chrom_sizes,
)

_REQUIRED = ("samtools", "faToTwoBit", "twoBitInfo")
_TOOLS_PRESENT = all(shutil.which(t) is not None for t in _REQUIRED)
_needs_tools = pytest.mark.skipif(
    not _TOOLS_PRESENT, reason="samtools/faToTwoBit/twoBitInfo not on PATH"
)

# A tiny two-sequence genome: chr1 is 16 bp, chr2 is 8 bp.
_CHR1 = "ACGTACGTACGTACGT"
_CHR2 = "TTTTGGGG"


@pytest.fixture
def fasta(tmp_path: Path) -> Path:
    path = tmp_path / "mini.fa"
    path.write_text(f">chr1\n{_CHR1}\n>chr2\n{_CHR2}\n")
    return path


def test_faidx_missing_input_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="not found"):
        faidx(tmp_path / "nope.fa")


def test_fasta_to_2bit_missing_input_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="not found"):
        fasta_to_2bit(tmp_path / "nope.fa")


def test_chrom_sizes_missing_input_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="not found"):
        twobit_to_chrom_sizes(tmp_path / "nope.2bit")


@_needs_tools
def test_faidx_creates_index(fasta: Path) -> None:
    fai = faidx(fasta)
    assert fai == fasta.with_name("mini.fa.fai")
    assert fai.is_file()
    # .fai columns: name, length, offset, linebases, linewidth
    names_lengths = {
        line.split("\t")[0]: int(line.split("\t")[1]) for line in fai.read_text().splitlines()
    }
    assert names_lengths == {"chr1": len(_CHR1), "chr2": len(_CHR2)}


@_needs_tools
def test_fasta_to_2bit_default_name(fasta: Path) -> None:
    twobit = fasta_to_2bit(fasta)
    assert twobit == fasta.with_name("mini.2bit")
    assert twobit.is_file()
    assert twobit.read_bytes()[:4] != b""  # non-empty binary output


@_needs_tools
def test_fasta_to_2bit_explicit_path(fasta: Path, tmp_path: Path) -> None:
    dest = tmp_path / "custom.2bit"
    twobit = fasta_to_2bit(fasta, dest)
    assert twobit == dest
    assert dest.is_file()


@_needs_tools
def test_prepare_fasta_end_to_end(fasta: Path) -> None:
    files = prepare_fasta(fasta)
    assert isinstance(files, GenomeFiles)
    assert files.fasta == fasta
    assert files.fai.is_file()
    assert files.twobit == fasta.with_name("mini.2bit")
    assert files.twobit.is_file()
    assert files.chrom_sizes == fasta.with_name("mini.chrom.sizes")
    assert files.chrom_sizes.is_file()

    sizes = {
        line.split("\t")[0]: int(line.split("\t")[1])
        for line in files.chrom_sizes.read_text().splitlines()
    }
    assert sizes == {"chr1": len(_CHR1), "chr2": len(_CHR2)}
