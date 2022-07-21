from __future__ import annotations
from math import ceil
from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
from scipy.cluster import hierarchy as sch


def svg_heatmap(
    hotspot_df: pd.DataFrame,
    coordinate_df: pd.DataFrame,
    cluster_result: pd.Series,
    save_path: Union[str, Path],
    he_image: Image.Image = None,
) -> None:
    """
    Draw SVG distribution heatmap.

    Patameters
    ==========
    hotspot_df: pd.DataFrame
        A pd.DataFrame for hotspots.

    coordinate_df: pd.DataFrame
        A pd.DataFrame for coordinate files.

    cluster_result: pd.Series
        A pd.Series for cluster result.

    save_path: str or pathlib.Path
        Heatmap save path.

    he_image: PIL.Image.Image, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.
    """

    left, bottom = 0.1, 0.1
    width, height = 0.66, 0.66
    spacing, cluster_width = 0.03, 0.22
    rect_heatmap = [left + spacing, bottom + spacing * 2, width, height]

    fig = plt.figure(figsize=(10, 10))
    ax_heatmap = fig.add_axes(rect_heatmap)
    axes = [ax_heatmap]

    for i, j in enumerate(set(cluster_result.values)):
        rect_cluster = [
            left + width + spacing * (2 + i // 3) + cluster_width * (i // 3),
            bottom + spacing * (3 - i % 3) + cluster_width * (2 - i % 3),
            cluster_width,
            cluster_width,
        ]
        ax_cluster = fig.add_axes(rect_cluster)
        ax_cluster.set_title(f"Cluster {j} hotspot mean")
        ax_cluster.set_xticks([])
        ax_cluster.set_yticks([])
        ax_cluster.imshow(he_image) if he_image is not None else None
        None if he_image is None else ax_cluster.imshow(he_image)
        ax_cluster.axis("off")
        sc = ax_cluster.scatter(
            coordinate_df.iloc[:, 0],
            coordinate_df.iloc[:, 1],
            c=hotspot_df[cluster_result[cluster_result == j].index].T.mean(),
            cmap="autumn_r",
            vmin=0,
            vmax=1,
            s=4,
            alpha=0.7,
        )
        axes.append(ax_cluster)

    ii = ceil(len(set(cluster_result.values)) / 3)
    rect_cb = [
        left + width + spacing * (ii + 2) + cluster_width * ii,
        bottom + spacing * 2,
        0.02,
        height,
    ]
    ax_cb = fig.add_axes(rect_cb)
    fig.colorbar(sc, cax=ax_cb)
    axes.append(ax_cb)

    spot_distmat = sch.distance.pdist(hotspot_df, metric="jaccard")
    Z_spot = sch.linkage(spot_distmat, method="ward")
    spot_result = pd.Series(
        sch.fcluster(Z_spot, t=8, criterion="maxclust"),
        index=hotspot_df.index,
    ).sort_values()

    ax_heatmap.imshow(
        hotspot_df.reindex(
            index=spot_result.index,
            columns=cluster_result.index,
        ),
        cmap="Reds",
        aspect="auto",
    )

    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])
    ax_heatmap.set_xlabel("Genes")
    ax_heatmap.set_ylabel("Spots")

    fig.savefig(save_path, bbox_inches="tight")
