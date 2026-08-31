"""liulab-genome: handling genomic files (metadata, processing, feature extraction)."""

from importlib.metadata import PackageNotFoundError, version

from genome.annotation import (
    AnnotationMetadata,
    AnnotationRegistry,
    CuratedGeneListError,
    GeneCategoryNotDeclaredError,
    GeneList,
    GeneListAssemblyMismatchError,
    NoGeneCategoriesError,
)
from genome.assembly.chimera_build import AmbiguousDefaultAnnotationError
from genome.assembly.genome import Genome
from genome.assembly.metadata import AssemblyMetadata
from genome.external import ToolNotFoundError
from genome.region import Region
from genome.seq import DNA, RNA, Protein
from genome.shipped import MetadataRowError, ShippedTableError
from genome.tf.cofactor import NoCofactorTableError, TFCofactorList
from genome.tf.gene import NoTFCensusError, TFGeneList
from genome.tf.species import UnknownSpeciesError

try:
    __version__ = version("liulab-genome")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"

__all__ = [
    "DNA",
    "RNA",
    "AmbiguousDefaultAnnotationError",
    "AnnotationMetadata",
    "AnnotationRegistry",
    "AssemblyMetadata",
    "CuratedGeneListError",
    "GeneCategoryNotDeclaredError",
    "GeneList",
    "GeneListAssemblyMismatchError",
    "Genome",
    "MetadataRowError",
    "NoCofactorTableError",
    "NoGeneCategoriesError",
    "NoTFCensusError",
    "Protein",
    "Region",
    "ShippedTableError",
    "TFCofactorList",
    "TFGeneList",
    "ToolNotFoundError",
    "UnknownSpeciesError",
    "__version__",
]
