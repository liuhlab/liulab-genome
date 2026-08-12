"""Tests for genome.Genome — sequence retrieval over a prepared assembly.

A real tiny assembly is built with the native tools, and the network download is
replaced by monkeypatching ``UCSCGenomeDownloader.fetch_genome`` to return those
prebuilt files. The whole module therefore skips when the native tools are absent.
"""

from __future__ import annotations

import shutil
from collections.abc import Iterator
from pathlib import Path

import pandas as pd
import pytest

import genome.genome as genome_mod
from genome import DNA, Genome, Region
from genome.io.fasta import prepare_fasta
from genome.metadata import METADATA_FIELDS, AssemblyMetadata

_REQUIRED = ("samtools", "faToTwoBit", "twoBitInfo")
_TOOLS_PRESENT = all(shutil.which(t) is not None for t in _REQUIRED)
pytestmark = pytest.mark.skipif(
    not _TOOLS_PRESENT, reason="samtools/faToTwoBit/twoBitInfo not on PATH"
)

# chrA: 20 bp with a soft-masked stretch at [8, 12); chrB: 8 bp.
_CHR_A = "ACGTACGTacgtACGTACGT"
_CHR_B = "TTTTGGGG"


@pytest.fixture
def prepared_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Prepare a tiny assembly on disk and return the cache dir to construct against.

    Skips the network/UCSC validation entirely: every ``Genome(...)`` built with
    this cache dir is handed back the prebuilt files, whatever it is named.
    """
    fasta = tmp_path / "tiny.fa"
    fasta.write_text(f">chrA\n{_CHR_A}\n>chrB\n{_CHR_B}\n")
    files = prepare_fasta(fasta)
    monkeypatch.setattr(
        genome_mod.UCSCGenomeDownloader,
        "fetch_genome",
        lambda self, **kwargs: files,
    )
    return tmp_path


@pytest.fixture
def genome(prepared_dir: Path) -> Iterator[Genome]:
    g = Genome("tiny", cache_dir=prepared_dir)
    yield g
    g.close()


def test_fetch_range_returns_dna(genome: Genome) -> None:
    result = genome.fetch_sequence("chrA:0-8")
    assert result == DNA("ACGTACGT")
    assert isinstance(result, DNA)


def test_fetch_preserves_soft_masking(genome: Genome) -> None:
    assert genome.fetch_sequence("chrA:8-12") == DNA("acgt")


def test_fetch_bare_chromosome_returns_whole_sequence(genome: Genome) -> None:
    assert genome.fetch_sequence("chrB") == DNA(_CHR_B)


def test_getitem_is_sugar_for_fetch_sequence(genome: Genome) -> None:
    assert genome["chrA:0-4"] == DNA("ACGT")


def test_fetch_accepts_region_object(genome: Genome) -> None:
    assert genome.fetch_sequence(Region("chrA", 0, 8)) == DNA("ACGTACGT")


def test_minus_strand_region_is_reverse_complemented(genome: Genome) -> None:
    # chrA[0:6] == "ACGTAC" — not a palindrome, so the reverse complement differs.
    plus = genome.fetch_sequence(Region("chrA", 0, 6, "+"))
    minus = genome.fetch_sequence(Region("chrA", 0, 6, "-"))
    assert plus == DNA("ACGTAC")
    assert minus == plus.reverse_complement() == DNA("GTACGT")


def test_end_equal_to_size_is_allowed(genome: Genome) -> None:
    assert genome.fetch_sequence("chrA:0-20") == DNA(_CHR_A)


def test_chrom_sizes_is_an_independent_series_copy(genome: Genome) -> None:
    sizes = genome.chrom_sizes
    assert isinstance(sizes, pd.Series)
    assert sizes.to_dict() == {"chrA": 20, "chrB": 8}
    sizes["chrA"] = 0  # mutate the returned copy...
    assert genome.chrom_sizes["chrA"] == 20  # ...the genome is unaffected


def test_chromosomes_preserve_reference_order(genome: Genome) -> None:
    assert genome.chromosomes == ["chrA", "chrB"]


def test_file_path_properties_mirror_genome_files(genome: Genome) -> None:
    assert genome.fasta_path == genome.files.fasta
    assert genome.twobit_path == genome.files.twobit
    assert genome.chrom_sizes_path == genome.files.chrom_sizes
    # chrom_sizes_path is the on-disk file, distinct from the in-memory sizes Series.
    assert isinstance(genome.chrom_sizes_path, Path)
    assert genome.chrom_sizes_path.is_file()
    assert isinstance(genome.chrom_sizes, pd.Series)


def test_repr(genome: Genome) -> None:
    assert repr(genome) == "Genome('tiny', 2 sequences)"


@pytest.mark.parametrize(
    ("region", "message"),
    [
        ("chrZ:0-5", "unknown chromosome"),
        ("chrA:0-21", "exceeds chrA length"),
        ("chrA:10-5", "is past end"),
        ("chrA:bad", "malformed region"),
    ],
)
def test_invalid_region_raises(genome: Genome, region: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        genome.fetch_sequence(region)


def test_context_manager_closes_handle(genome: Genome) -> None:
    with genome as g:
        assert g.fetch_sequence("chrA:0-4") == DNA("ACGT")
    with pytest.raises(ValueError, match="closed"):
        g.fetch_sequence("chrA:0-4")


_OVERRIDE = AssemblyMetadata(
    assembly_name="tinyAsm",
    species="Testus minimus",
    ucsc_name="tinyUcsc",
    ncbi_name="TINY.1",
    ncbi_assembly_id="GCF_000000000.0",
    ncbi_taxid=1,
)


def test_metadata_record_replaces_the_curated_row(prepared_dir: Path) -> None:
    # sacCer3 is listed in the shipped table; a record given here wins over every field of it.
    with Genome("sacCer3", cache_dir=prepared_dir, metadata=_OVERRIDE) as g:
        assert g.metadata == _OVERRIDE
        assert [getattr(g, field) for field in METADATA_FIELDS] == [
            "tinyAsm",
            "Testus minimus",
            "tinyUcsc",
            "TINY.1",
            "GCF_000000000.0",
            1,
        ]


def test_metadata_falls_back_to_the_curated_table(prepared_dir: Path) -> None:
    with Genome("sacCer3", cache_dir=prepared_dir) as g:
        assert g.assembly_name == "sacCer3"
        assert g.species == "Saccharomyces cerevisiae"
        assert g.ucsc_name == "sacCer3"
        assert g.ncbi_name == "R64-1-1"
        assert g.ncbi_assembly_id == "GCF_000146045.2"
        assert g.ncbi_taxid == 559292


def test_assembly_absent_from_the_table_still_constructs(prepared_dir: Path) -> None:
    # The table is a cross-reference, not an allow-list: every identifier is simply unknown.
    with Genome("tiny", cache_dir=prepared_dir) as g:
        assert g.metadata is None
        assert all(getattr(g, field) is None for field in METADATA_FIELDS)
        assert g.chromosomes == ["chrA", "chrB"]


def test_path_or_url_seeds_from_local_fasta(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # A local FASTA seed prepares the assembly end-to-end (copy + native tools)
    # without ever contacting UCSC.
    src = tmp_path / "src.fa"
    src.write_text(f">chrA\n{_CHR_A}\n>chrB\n{_CHR_B}\n")

    def boom(self: object, **_kwargs: object) -> None:
        raise AssertionError("UCSC fetch_genome must not run when path_or_url is given")

    monkeypatch.setattr(genome_mod.UCSCGenomeDownloader, "fetch_genome", boom)

    with Genome("tiny", path_or_url=src, cache_dir=tmp_path / "cache") as g:
        assert g.fetch_sequence("chrA:0-8") == DNA("ACGTACGT")
        assert g.chromosomes == ["chrA", "chrB"]
