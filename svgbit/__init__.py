from importlib.metadata import version, PackageNotFoundError

from .core.STDataset import STDataset
from .run import run

from .core.io import load_10X, load_anndata_h5, load_table
from .core.combinations import find_combinations
from . import filters, normalizers, plot

try:
    __version__ = version("svgbit")
except PackageNotFoundError:
    __version__ = "unknown version"

__all__ = [
    "STDataset",
    "run",
    "load_10X",
    "load_anndata_h5",
    "load_table",
    "filters",
    "normalizers",
    "plot",
    "find_combinations",
]
