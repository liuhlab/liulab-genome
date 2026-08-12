"""File readers, writers, and downloaders — the I/O boundary.

Modules here are where the package talks to the outside world: the network
(:mod:`genome.io.download`) and on-disk genomic files plus the native binaries
that process them (:mod:`genome.io.fasta`). Keep real, side-effect-free logic
in ``core``/``features``; this layer only moves bytes.
"""

from genome.io.completion import (
    CompletionRecord,
    FileDisagreement,
    RegistrationError,
    RegistrationMismatchError,
    UnfinishedRegistrationError,
    build_record,
    check_registration,
    disagreements,
    read_record,
    write_record,
)
from genome.io.download import (
    Downloader,
    UCSCGenomeDownloader,
    assembly_data_dir,
    assembly_table_row,
    liulab_data_dir,
    register_assembly,
    verify_assembly,
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
    fetch_annotation,
    list_annotations,
    register_annotation,
    register_gtf,
)
from genome.io.twobit import TwoBit
from genome.io.utils import ChecksumMismatchError, sha256_file

__all__ = [
    "ChecksumMismatchError",
    "CompletionRecord",
    "Downloader",
    "FileDisagreement",
    "GenomeFiles",
    "GtfAnnotation",
    "RegistrationError",
    "RegistrationMismatchError",
    "TwoBit",
    "UCSCGenomeDownloader",
    "UnfinishedRegistrationError",
    "annotation_dir",
    "assembly_data_dir",
    "assembly_table_row",
    "build_record",
    "check_registration",
    "disagreements",
    "faidx",
    "fasta_to_2bit",
    "fetch_annotation",
    "list_annotations",
    "liulab_data_dir",
    "prepare_fasta",
    "read_chrom_sizes",
    "read_record",
    "register_annotation",
    "register_assembly",
    "register_gtf",
    "sha256_file",
    "twobit_to_chrom_sizes",
    "verify_assembly",
    "write_record",
]
