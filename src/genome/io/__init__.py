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
    RegisteredAssembly,
    UCSCGenomeDownloader,
    VerifiedAssembly,
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

# Four names the registry absorbed are deliberately absent — `list_annotations`,
# `list_broken_annotations`, `default_annotation` and `fetch_annotation`. Each is now
# reached through `AnnotationRegistry`, nothing in the package calls them by name, and
# what a caller wants from an assembly's annotations is the registry rather than the
# four steps it takes internally. They stay importable from `genome.io.gtf` for the
# tests that hold each rule on its own.
from genome.io.gtf import (
    AnnotationNotRegisteredError,
    AnnotationRegistry,
    AnnotationStatus,
    AnnotationStatusRow,
    BrokenAnnotation,
    ChromosomeMismatchError,
    GtfAnnotation,
    RegisteredAnnotation,
    annotation_dir,
    annotation_status,
    chromosome_check_summary,
    register_annotation,
    register_annotation_by_path,
    register_gtf,
)
from genome.io.twobit import TwoBit
from genome.io.utils import ChecksumMismatchError, sha256_file

__all__ = [
    "AnnotationNotRegisteredError",
    "AnnotationRegistry",
    "AnnotationStatus",
    "AnnotationStatusRow",
    "BrokenAnnotation",
    "ChecksumMismatchError",
    "ChromosomeMismatchError",
    "CompletionRecord",
    "Downloader",
    "FileDisagreement",
    "GenomeFiles",
    "GtfAnnotation",
    "RegisteredAnnotation",
    "RegisteredAssembly",
    "RegistrationError",
    "RegistrationMismatchError",
    "TwoBit",
    "UCSCGenomeDownloader",
    "UnfinishedRegistrationError",
    "VerifiedAssembly",
    "annotation_dir",
    "annotation_status",
    "assembly_data_dir",
    "assembly_table_row",
    "build_record",
    "check_registration",
    "chromosome_check_summary",
    "disagreements",
    "faidx",
    "fasta_to_2bit",
    "liulab_data_dir",
    "prepare_fasta",
    "read_chrom_sizes",
    "read_record",
    "register_annotation",
    "register_annotation_by_path",
    "register_assembly",
    "register_gtf",
    "sha256_file",
    "twobit_to_chrom_sizes",
    "verify_assembly",
    "write_record",
]
