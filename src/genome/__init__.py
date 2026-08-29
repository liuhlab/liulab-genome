"""liulab-genome: handling genomic files (metadata, processing, feature extraction)."""

from importlib.metadata import PackageNotFoundError, version

from genome.gene_list import GeneCategoryNotDeclaredError, NoGeneCategoriesError
from genome.genome import Genome
from genome.io.gtf import AnnotationRegistry, NoTFCensusError, UnknownSpeciesError
from genome.io.results import GeneList, TFGeneList
from genome.metadata import AnnotationMetadata, AssemblyMetadata
from genome.region import Region
from genome.seq import DNA, RNA, Protein

try:
    __version__ = version("liulab-genome")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"

__all__ = [
    "DNA",
    "RNA",
    "AnnotationMetadata",
    "AnnotationRegistry",
    "AssemblyMetadata",
    "GeneCategoryNotDeclaredError",
    "GeneList",
    "Genome",
    "NoGeneCategoriesError",
    "NoTFCensusError",
    "Protein",
    "Region",
    "TFGeneList",
    "UnknownSpeciesError",
    "__version__",
]
