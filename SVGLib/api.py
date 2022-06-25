from multiprocessing import cpu_count

from .core.STDataset import STDataset


def run(
    dataset: STDataset,
    k: int = 6,
    n_genes: int = 1000,
    n_gene_clusters: int = 8,
    cores: int = cpu_count,
) -> STDataset:
    dataset.acquire_weight(k=k)
    dataset.acquire_hotspot(cores=cores)
    dataset.acquire_density(cores=cores)
    dataset.find_clusters(n_genes=n_genes, n_gene_clusters=n_gene_clusters)

    return dataset
