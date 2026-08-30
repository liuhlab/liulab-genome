"""liulab-genome: handling genomic files (metadata, processing, feature extraction)."""

from importlib.metadata import PackageNotFoundError, version

from genome.gene_list import GeneCategoryNotDeclaredError, NoGeneCategoriesError
from genome.genome import Genome
from genome.io.gtf import AnnotationRegistry, GeneList
from genome.metadata import AnnotationMetadata, AssemblyMetadata
from genome.region import Region
from genome.seq import DNA, RNA, Protein
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
    "AnnotationMetadata",
    "AnnotationRegistry",
    "AssemblyMetadata",
    "GeneCategoryNotDeclaredError",
    "GeneList",
    "Genome",
    "NoCofactorTableError",
    "NoGeneCategoriesError",
    "NoTFCensusError",
    "Protein",
    "Region",
    "TFCofactorList",
    "TFGeneList",
    "UnknownSpeciesError",
    "__version__",
]
