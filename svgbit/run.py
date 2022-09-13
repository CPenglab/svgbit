from __future__ import annotations

from multiprocessing import cpu_count

from .core.STDataset import STDataset


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
