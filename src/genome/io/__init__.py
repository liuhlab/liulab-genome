"""File readers, writers, and downloaders — the I/O boundary.

Modules here are where the package talks to the outside world: the network
(:mod:`genome.io.download`) and on-disk genomic files plus the native binaries
that process them (:mod:`genome.io.fasta`). Keep real, side-effect-free logic
in ``core``/``features``; this layer only moves bytes.
"""

from genome.io.download import (
    Downloader,
    UCSCGenomeDownloader,
    assembly_data_dir,
    liulab_data_dir,
)
from genome.io.fasta import (
    GenomeFiles,
    faidx,
    fasta_to_2bit,
    prepare_fasta,
    read_chrom_sizes,
    twobit_to_chrom_sizes,
)
from genome.io.gtf import (
    GtfAnnotation,
    annotation_dir,
    list_annotations,
    register_gtf,
)
from genome.io.twobit import TwoBit

__all__ = [
    "Downloader",
    "GenomeFiles",
    "GtfAnnotation",
    "TwoBit",
    "UCSCGenomeDownloader",
    "annotation_dir",
    "assembly_data_dir",
    "faidx",
    "fasta_to_2bit",
    "list_annotations",
    "liulab_data_dir",
    "prepare_fasta",
    "read_chrom_sizes",
    "register_gtf",
    "twobit_to_chrom_sizes",
]
