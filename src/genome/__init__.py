"""liulab-genome: handling genomic files (metadata, processing, feature extraction)."""

from importlib.metadata import PackageNotFoundError, version

from genome.seq import DNA, RNA, Protein

try:
    __version__ = version("liulab-genome")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"

__all__ = ["DNA", "RNA", "Protein", "__version__"]
