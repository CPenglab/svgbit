from collections import Counter
from typing import Optional, Union

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.mixture import BayesianGaussianMixture

from .STDataset import STDataset


def _find_combinations(
    hotspot_df: pd.DataFrame,
    coordinate_df: pd.DataFrame,
    spot_type: pd.DataFrame,
    svg_cluster: pd.Series,
    center_spots: Union[int, list],
    selected_genes: Optional[list] = None,
    use_neighbor: bool = True,
) -> pd.DataFrame:
    """
    Find gene pairs in certain SVG cluster.

    Parameters
    ==========
    hotspot_df : pd.DataFrame
        A hotspot DataFrame generated by svgbit.

    coordinate_df : pd.DataFrame
        A pd.DataFrame for coordinate files.

    spot_type : pd.DataFrame
        A pd.DataFrame for assigned type_df.

    center_spots : int or list
        If a ``int`` is given, find gene pairs for this SVG cluster. If a
        ``list`` of spots is given, find gene pairs within those spots.
        ``selected_genes`` should be given if ``center_cluster`` is a ``list``.

    selected_genes : list or None, default None
        If a ``list`` of genes is given, find gene pairs within given genes.

    use_neighbor : bool, default True
        Whether to find gene pairs with neighbor SVG clusters.

    Returns
    =======
    gene_pairs_df : pd.DataFrame
        A pd.DataFrame for gene pairs in center_cluster with weights.
    """
    try:
        center_spots = list(center_spots)
    except TypeError:
        try:
            center_spots = int(center_spots)
        except TypeError:
            err = (f"center_spots must be a number or a list-like object, "
                   "not '{type(center_spots)}'")
            raise TypeError(err)

    if use_neighbor and (not isinstance(center_spots, list)):
        # find neighbor clusters
        certain_series = spot_type[spot_type["spot_type"] != "uncertain"]
        certain_series = certain_series["type_1"]
        nbrs = NearestNeighbors(n_neighbors=7).fit(coordinate_df)
        neighbor_clusters_ = {}
        for present_cluster in set(spot_type["type_1"]):
            used_spots = []
            neighbor_clusters_[present_cluster] = []
            present_spots = certain_series[certain_series ==
                                           present_cluster].index
            distances, indices = nbrs.kneighbors(
                coordinate_df.reindex(index=present_spots))
            for record in indices:
                for neighbor_spot in record[1:]:
                    if neighbor_spot in used_spots:
                        continue
                    used_spots.append(neighbor_spot)
                    neighbor_cluster = spot_type.iloc[neighbor_spot, 1]
                    if neighbor_cluster == present_cluster:
                        continue
                    neighbor_clusters_[present_cluster].append(
                        neighbor_cluster)

        neighbor_clusters = {}
        for i in neighbor_clusters_:
            counter = Counter(neighbor_clusters_[i])
            coverage_ = len(certain_series[certain_series == i]) * 0.03
            neighbor_clusters[i] = [
                j[0] for j in counter.most_common()[:3] if j[1] >= coverage_
            ]
            neighbor_clusters[i].append(i)

        selected_spots = [
            j for i in neighbor_clusters[center_spots]
            for j in spot_type[spot_type["type_1"] == i].index
        ]
        selected_genes = [
            j for i in neighbor_clusters[center_spots]
            for j in svg_cluster[svg_cluster == i].index
        ]
    else:
        if isinstance(center_spots, list):
            selected_spots = center_spots
            if len(center_spots) <= 10:
                center_spots = ", ".join(center_spots)
            else:
                center_spots = ", ".join(center_spots[:3])
                center_spots = center_spots + "..."
        else:
            selected_spots = spot_type[spot_type["type_1"] ==
                                       center_spots].index
            if selected_genes is None:
                selected_genes = svg_cluster[svg_cluster == center_spots].index

    count_sub = hotspot_df.reindex(
        index=selected_spots,
        columns=selected_genes,
    )
    con_matrix = count_sub.T @ count_sub
    n_hotspots = np.diag(con_matrix.to_numpy())
    con_matrix = con_matrix - np.diag(n_hotspots)
    men_matrix = (1 - count_sub.T) @ count_sub

    con_ratio = con_matrix / len(selected_spots)
    men_ratio = men_matrix / len(selected_spots)

    gene_pairs_df = pd.DataFrame(columns=[
        "SVG_cluster", "gene_1", "gene_2", "colocalization_score",
        "exclusive_score"
    ])
    line = 0
    for gene_1 in con_ratio.columns:
        for gene_2 in con_ratio.index:
            if gene_1 == gene_2:
                continue
            write_dict = {
                "SVG_cluster": center_spots,
                "gene_1": gene_1,
                "gene_2": gene_2,
                "colocalization_score": con_ratio[gene_1][gene_2],
                "exclusive_score": men_ratio[gene_1][gene_2],
            }
            write_series = pd.Series(write_dict, name=line).to_frame().T
            line += 1
            gene_pairs_df = pd.concat([gene_pairs_df, write_series])

    gene_pairs_df.fillna(0, inplace=True)

    gmm_con = BayesianGaussianMixture(
        n_components=3,
        max_iter=500,
        n_init=5,
    )
    result_con = gmm_con.fit_predict(
        gene_pairs_df["colocalization_score"].to_numpy().reshape(-1, 1))
    rank_str = ["low", "middle", "high"]
    mean_argsort = gmm_con.means_.argsort(axis=0).flatten()
    rank_dict = {i: j for i, j in zip(mean_argsort, rank_str)}
    result_con = [rank_dict[i] for i in result_con]
    gene_pairs_df["colocalization_degree"] = result_con

    gmm_men = BayesianGaussianMixture(
        n_components=3,
        max_iter=500,
        n_init=5,
    )
    result_men = gmm_men.fit_predict(
        gene_pairs_df["exclusive_score"].to_numpy().reshape(-1, 1))
    mean_argsort = gmm_men.means_.argsort(axis=0).flatten()
    rank_dict = {i: j for i, j in zip(mean_argsort, rank_str)}
    result_men = [rank_dict[i] for i in result_men]
    gene_pairs_df["exclusive_degree"] = result_men

    return gene_pairs_df


def find_combinations(
    dataset: STDataset,
    center_spots: Union[int, list],
    selected_genes: Optional[list] = None,
    use_neighbor: bool = True,
) -> pd.DataFrame:
    """
    Find gene pairs in certain SVG cluster.

    Parameters
    ==========
    dataset : STDataset
        A STDataset with all steps finished.

    center_spots : int or list
        If a ``int`` is given, find gene pairs for this SVG cluster. If a
        ``list`` of spots is given, find gene pairs within those spots.
        ``selected_genes`` should be given if ``center_cluster`` is a ``list``.

    selected_genes : list or None, default None
        If a ``list`` of genes is given, find gene pairs within given genes.

    use_neighbor : bool, default True
        Whether to find gene pairs with neighbor SVG clusters.

    Returns
    =======
    gene_pairs_df : pd.DataFrame
        A pd.DataFrame for gene pairs in center_cluster with weights.
    """
    return _find_combinations(
        dataset.hotspot_df,
        dataset.coordinate_df,
        dataset.spot_type,
        dataset.svg_cluster,
        center_spots,
        selected_genes,
        use_neighbor,
    )
