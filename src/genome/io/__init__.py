"""File readers, writers, and downloaders — the I/O boundary.

Modules here are where the package talks to the outside world: the network
(:mod:`genome.io.download`) and on-disk genomic files plus the native binaries
that process them (:mod:`genome.io.fasta`). Keep real, side-effect-free logic
in ``core``/``features``; this layer only moves bytes.
"""

from genome.io.download import Downloader, UCSCGenomeDownloader
from genome.io.fasta import (
    GenomeFiles,
    faidx,
    fasta_to_2bit,
    prepare_fasta,
    twobit_to_chrom_sizes,
)

__all__ = [
    "Downloader",
    "GenomeFiles",
    "UCSCGenomeDownloader",
    "faidx",
    "fasta_to_2bit",
    "prepare_fasta",
    "twobit_to_chrom_sizes",
]
