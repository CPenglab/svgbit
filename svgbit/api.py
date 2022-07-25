from __future__ import annotations

from multiprocessing import cpu_count
from pathlib import Path
from typing import Optional, Union

from .core.STDataset import STDataset
from .core.plot import _svg_heatmap


def run(
        dataset: STDataset,
        k: int = 6,
        n_svgs: int = 1000,
        n_svg_clusters: int = 8,
        cores: int = cpu_count(),
) -> STDataset:
    """
    Run all part of svgbit within one function.

    Parameters
    ==========
    dataset : STDataset
        A STDataset for running svgbit.

    k : int, default 6
        Number of nearest neighbors for KNN network.

    n_svgs : int, default 1000
        Number of SVGs to find clusters.

    n_svg_clusters : int, default 8
        Number of SVG clusters to find.

    cores : int
        Number of threads to run svgbit. Use all available cpus by default.

    Returns
    =======
    dataset : STDataset
        A STDataset with all evaluates done.

    """
    dataset.acquire_weight(k=k)
    dataset.acquire_hotspot(cores=cores)
    dataset.acquire_density(cores=cores)
    dataset.find_clusters(n_svgs=n_svgs, n_svg_clusters=n_svg_clusters)

    return dataset


def svg_heatmap(
    dataset: STDataset,
    save_path: Union[str, Path],
    he_image: Optional[Union[str, Path]] = None,
) -> None:
    """
    Draw SVG distribution heatmap.

    Patameters
    ==========
    dataset : STDataset
        A STDataset with hotspot and SVG cluster estimation finished.

    save_path : str or pathlib.Path
        Heatmap save path.

    he_image : str or pathlib.Path, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.
    """
    _svg_heatmap(
        dataset.hotspot_df,
        dataset.coordinate_df,
        dataset.svg_cluster,
        save_path,
        he_image,
    )
