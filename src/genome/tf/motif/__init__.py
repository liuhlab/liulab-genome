"""TF binding motifs — the matrices themselves, and finding where they occur."""

from genome.tf.motif.background import (
    BACKGROUND_FLOOR,
    BACKGROUND_MODES,
    UNIFORM_BACKGROUND,
    BackgroundArg,
)
from genome.tf.motif.compare import (
    COMPARISON_VARIABLES,
    FRAME_COLUMNS,
    MotifComparison,
    RaggedComparisonError,
)
from genome.tf.motif.jaspar import (
    DEFAULT_RELEASE,
    DEFAULT_TAX_GROUP,
    JASPAR_RELEASES,
    JASPAR_TAX_GROUPS,
    JasparDatabase,
    JasparReleaseError,
    TransfacError,
    parse_transfac,
)
from genome.tf.motif.mixin import REGION_HIT_PROVENANCE, MotifScanMixin
from genome.tf.motif.motif import (
    BASES,
    DEFAULT_SEQUENCE_NAME,
    DEFAULT_THRESHOLD,
    MIN_MOTIF_LENGTH,
    AmbiguousBaseIdError,
    AmbiguousMotifNameError,
    Motif,
    MotifNotFoundError,
    MotifSet,
)
from genome.tf.motif.parquet import HIT_PROVENANCE_KEY, hit_count, provenance_of, read_hits
from genome.tf.motif.scan import (
    HIT_COLUMNS,
    HIT_DTYPES,
    HIT_PROVENANCE,
    FastaFormatError,
)
from genome.tf.motif.thresholds import threshold_cache_dir

# Worker resolution lives at the package root — it is lab-HPC infrastructure and names
# nothing in this context. Re-exported here so `from genome.tf.motif import
# resolve_workers` keeps answering, since a scan is what most callers ask it about.
from genome.workers import DEFAULT_WORKERS, SLURM_CPU_VARS, resolve_workers

__all__ = [
    "BACKGROUND_FLOOR",
    "BACKGROUND_MODES",
    "BASES",
    "COMPARISON_VARIABLES",
    "DEFAULT_RELEASE",
    "DEFAULT_SEQUENCE_NAME",
    "DEFAULT_TAX_GROUP",
    "DEFAULT_THRESHOLD",
    "DEFAULT_WORKERS",
    "FRAME_COLUMNS",
    "HIT_COLUMNS",
    "HIT_DTYPES",
    "HIT_PROVENANCE",
    "HIT_PROVENANCE_KEY",
    "JASPAR_RELEASES",
    "JASPAR_TAX_GROUPS",
    "MIN_MOTIF_LENGTH",
    "REGION_HIT_PROVENANCE",
    "SLURM_CPU_VARS",
    "UNIFORM_BACKGROUND",
    "AmbiguousBaseIdError",
    "AmbiguousMotifNameError",
    "BackgroundArg",
    "FastaFormatError",
    "JasparDatabase",
    "JasparReleaseError",
    "Motif",
    "MotifComparison",
    "MotifNotFoundError",
    "MotifScanMixin",
    "MotifSet",
    "RaggedComparisonError",
    "TransfacError",
    "hit_count",
    "parse_transfac",
    "provenance_of",
    "read_hits",
    "resolve_workers",
    "threshold_cache_dir",
]
