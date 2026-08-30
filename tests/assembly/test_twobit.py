"""Tests for genome.assembly.twobit.TwoBit.

The query tests need a real 2bit file, so they build one with the native tools
(``faToTwoBit`` etc.) and skip cleanly when those are absent. The open-error
tests need only ``py2bit`` (always present in the env) and always run.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from genome.assembly.fasta import prepare_fasta
from genome.assembly.twobit import TwoBit

_REQUIRED = ("samtools", "faToTwoBit", "twoBitInfo")
_TOOLS_PRESENT = all(shutil.which(t) is not None for t in _REQUIRED)
_needs_tools = pytest.mark.skipif(
    not _TOOLS_PRESENT, reason="samtools/faToTwoBit/twoBitInfo not on PATH"
)

# chrA: 20 bp with a soft-masked (lower-case) stretch at [8, 12); chrB: 8 bp.
_CHR_A = "ACGTACGTacgtACGTACGT"
_CHR_B = "TTTTGGGG"


@pytest.fixture
def twobit_path(tmp_path: Path) -> Path:
    fasta = tmp_path / "mini.fa"
    fasta.write_text(f">chrA\n{_CHR_A}\n>chrB\n{_CHR_B}\n")
    return prepare_fasta(fasta).twobit


# --- open errors (no native tools required) ---


def test_open_errors_name_the_problem(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="not found"):
        TwoBit(tmp_path / "nope.2bit")

    bad = tmp_path / "bad.2bit"
    bad.write_text("definitely not a 2bit file")
    with pytest.raises(RuntimeError, match="could not open"):
        TwoBit(bad)


# --- queries ---


@_needs_tools
def test_chroms_and_whole_chromosome_sequence(twobit_path: Path) -> None:
    with TwoBit(twobit_path) as tb:
        assert tb.chroms() == {"chrA": 20, "chrB": 8}
        assert tb.sequence("chrA") == _CHR_A
        assert tb.sequence("chrB") == _CHR_B


@_needs_tools
def test_sequence_is_zero_based_half_open(twobit_path: Path) -> None:
    with TwoBit(twobit_path) as tb:
        assert tb.sequence("chrA", 0, 8) == "ACGTACGT"
        assert tb.sequence("chrA", 4, 4) == ""  # empty half-open interval
        assert tb.sequence("chrA", 0, 20) == _CHR_A  # end == length is allowed


@_needs_tools
def test_masked_controls_soft_masking(twobit_path: Path) -> None:
    with TwoBit(twobit_path) as tb:
        assert tb.sequence("chrA", 8, 12) == "acgt"
    with TwoBit(twobit_path, masked=False) as tb:
        assert tb.sequence("chrA", 8, 12) == "ACGT"


@_needs_tools
def test_sequence_bounds_checked(twobit_path: Path) -> None:
    cases = [
        (("chrA", 0, 21), "exceeds sequence length"),  # boundary: one past the end
        (("chrA", -1, 5), "start must be >= 0"),
        (("chrA", 15, 10), "is past end"),
        (("missing", 0, 5), "unknown sequence"),
    ]
    with TwoBit(twobit_path) as tb:
        for args, message in cases:
            with pytest.raises(ValueError, match=message):
                tb.sequence(*args)


@_needs_tools
def test_nocheck_sequence_skips_bounds_check(twobit_path: Path) -> None:
    with TwoBit(twobit_path) as tb:
        # py2bit silently clamps an over-long end — nocheck does NOT raise.
        assert tb.nocheck_sequence("chrA", 0, 1000) == _CHR_A
        # ...but an unknown chromosome is still rejected.
        with pytest.raises(ValueError, match="unknown sequence"):
            tb.nocheck_sequence("missing")


# --- lifecycle ---


@_needs_tools
def test_close_is_idempotent_and_the_context_manager_closes_on_exit(twobit_path: Path) -> None:
    with TwoBit(twobit_path) as tb:
        assert repr(tb).endswith("open)")
        assert tb.sequence("chrA", 0, 4) == "ACGT"
    with pytest.raises(ValueError, match="closed"):
        tb.sequence("chrA", 0, 4)

    tb = TwoBit(twobit_path)
    tb.close()
    tb.close()  # second close is a no-op, not an error
    assert "closed" in repr(tb)
    with pytest.raises(ValueError, match="closed"):
        tb.sequence("chrA", 0, 4)
    with pytest.raises(ValueError, match="closed"):
        tb.chroms()
