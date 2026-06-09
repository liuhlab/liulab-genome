"""liulab-genome: handling genomic files (metadata, processing, feature extraction)."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("liulab-genome")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
