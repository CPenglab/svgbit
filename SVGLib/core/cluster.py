import pandas as pd
from scipy.cluster import hierarchy as sch


def cluster(
    hotspot_df: pd.DataFrame,
    AI_series: pd.Series,
    n_genes: int = 1000,
    n_gene_clusters: int = 8,
) -> pd.Series:
    selected_genes = AI_series.sort_values(ascending=False)[:n_genes].index
    hotspot_set = hotspot_df[selected_genes]
    gene_distmat = sch.distance.pdist(hotspot_set.T, metric="jaccard")
    Z_gene = sch.linkage(gene_distmat, method="ward")
    gene_result = pd.Series(
        sch.fcluster(Z_gene, t=n_gene_clusters, criterion="maxclust"),
        index=hotspot_df.index,
    ).sort_values()
    return gene_result
