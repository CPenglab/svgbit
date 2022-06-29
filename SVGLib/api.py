from multiprocessing import cpu_count

from .core.STDataset import STDataset


def run(
        dataset: STDataset,
        k: int = 6,
        n_genes: int = 1000,
        n_gene_clusters: int = 8,
        cores: int = cpu_count(),
) -> STDataset:
    """
    Run all part of SVGLib within one function.

    Parameters
    ==========
    dataset : STDataset
        A STDataset for running SVGLib.

    k: int, default 6
        Number of nearest neighbors for KNN network.

    n_genes: int, default 1000
        Number of genes to find clusters.

    n_gene_clusters: int, default 8
        Number of gene clusters to find.

    cores: int
        Number of threads to run SVGLib. Use all available cpus by default.

    Returns
    =======
    dataset : STDataset
        A STDataset with all evaluates done.

    """
    dataset.acquire_weight(k=k)
    dataset.acquire_hotspot(cores=cores)
    dataset.acquire_density(cores=cores)
    dataset.find_clusters(n_genes=n_genes, n_gene_clusters=n_gene_clusters)

    return dataset
