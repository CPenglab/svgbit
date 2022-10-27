from __future__ import annotations

from math import ceil
from pathlib import Path
from typing import Optional, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
from seaborn import color_palette

from .STDataset import STDataset


def _svg_heatmap(
    hotspot_df: pd.DataFrame,
    coordinate_df: pd.DataFrame,
    cluster_result: pd.Series,
    spot_type: pd.DataFrame,
    save_path: Union[str, Path],
    he_image: Optional[Union[str, Path]] = None,
    s: float = 4,
    dpi: float = 300,
) -> None:
    """
    Draw SVG distribution heatmap.

    Parameters
    ==========
    hotspot_df : pd.DataFrame
        A pd.DataFrame for hotspots.

    coordinate_df : pd.DataFrame
        A pd.DataFrame for coordinate files.

    cluster_result : pd.Series
        A pd.Series for cluster result.

    spot_type : pd.DataFrame
        A pd.DataFrame for assigned type_df.

    save_path : str or pathlib.Path
        Heatmap save path.

    he_image : str or pathlib.Path, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.
    """

    spot_type = spot_type.sort_values(by="type_1")

    left, bottom = 0.1, 0.1
    width, height = 0.66, 0.66
    spacing, cluster_width = 0.03, 0.22
    rect_heatmap = [left + spacing, bottom + spacing * 2, width, height]

    fig = plt.figure(figsize=(10, 10), dpi=dpi)
    ax_heatmap = fig.add_axes(rect_heatmap)

    if he_image is not None:
        he_image = Image.open(he_image)

    for i, j in enumerate(set(cluster_result.values)):
        rect_cluster = [
            left + width + spacing * (2 + i // 3) + cluster_width * (i // 3),
            bottom + spacing * (3 - i % 3) + cluster_width * (2 - i % 3),
            cluster_width,
            cluster_width,
        ]
        ax_cluster = fig.add_axes(rect_cluster)
        ax_cluster.set_title(f"Cluster {j} hotspot distribution")
        ax_cluster.set_xticks([])
        ax_cluster.set_yticks([])
        ax_cluster.imshow(he_image) if he_image is not None else None
        ax_cluster.axis("off")
        sc = ax_cluster.scatter(
            coordinate_df.iloc[:, 0],
            coordinate_df.iloc[:, 1],
            c=hotspot_df[cluster_result[cluster_result == j].index].T.mean(),
            cmap="autumn_r",
            vmin=0,
            vmax=1,
            s=s,
            alpha=0.7,
        )

    ii = ceil(len(set(cluster_result.values)) / 3)
    rect_cb = [
        left + width + spacing * (ii + 2) + cluster_width * ii,
        bottom + spacing * 2,
        0.02,
        height,
    ]
    ax_cb = fig.add_axes(rect_cb)
    fig.colorbar(sc, cax=ax_cb)

    ax_heatmap.imshow(
        hotspot_df.reindex(
            index=spot_type.index,
            columns=cluster_result.index,
        ),
        cmap="Reds",
        aspect="auto",
    )

    flag = 0
    for i in range(1, len(set(cluster_result.values))):
        flag += len(cluster_result[cluster_result == i])
        ax_heatmap.plot([flag, flag], [0, hotspot_df.shape[0] - 5])

    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])
    ax_heatmap.set_xlabel("Genes")
    ax_heatmap.set_ylabel("Spots")

    fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)
    he_image.close() if he_image is not None else None


def _hotspot_distribution_map(
    hotspot_df: pd.DataFrame,
    coordinate_df: pd.DataFrame,
    cluster_result: pd.Series,
    cluster: Union[str, int],
    save_path: Union[str, Path],
    he_image: Optional[Union[str, Path]] = None,
    s: float = 4,
    dpi: float = 300,
) -> None:
    """
    Draw hotspot distribution map for one SVG cluster.

    Parameters
    ==========
    hotspot_df : pd.DataFrame
        A pd.DataFrame for hotspots.

    coordinate_df : pd.DataFrame
        A pd.DataFrame for coordinate files.

    cluster_result : pd.Series
        A pd.Series for cluster result.

    cluster : str or int
        Specify drawing cluster.

    save_path : str or pathlib.Path
        Heatmap save path.

    he_image : str or pathlib.Path, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.
    """

    fig, ax = plt.subplots(figsize=(10, 10), dpi=dpi)
    cluster = int(cluster)

    ax.set_title(f"Cluster {cluster} hotspot distribution")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")
    sc = ax.scatter(
        coordinate_df.iloc[:, 0],
        coordinate_df.iloc[:, 1],
        c=hotspot_df[cluster_result[cluster_result == cluster].index].T.mean(),
        cmap="autumn_r",
        vmin=0,
        vmax=1,
        s=s,
        alpha=0.7,
    )
    fig.colorbar(sc)

    if he_image is not None:
        with Image.open(he_image) as he_image:
            ax.imshow(he_image)

    fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)


def _spot_type_map(
    hotspot_df: pd.DataFrame,
    coordinate_df: pd.DataFrame,
    type_df: pd.DataFrame,
    save_path: Union[str, Path],
    he_image: Optional[Union[str, Path]] = None,
    draw_uncertain: bool = True,
    s: float = 16,
    dpi: float = 300,
) -> None:
    """
    Draw SVG type map.

    Parameters
    ==========
    hotspot_df : pd.DataFrame
        A pd.DataFrame for hotspots.

    coordinate_df : pd.DataFrame
        A pd.DataFrame for coordinate files.

    type_df : pd.DataFrame
        A pd.DataFrame for spot type.

    save_path : str or pathlib.Path
        Type map save path.

    he_image : str or pathlib.Path, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    draw_uncertain : bool, default True
        Whether to draw uncertain spots.

    s : float, default 16
        Spot size.

    dpi : float, default 300
        DPI for saved figure.
    """

    if draw_uncertain:
        certain_spots = type_df.index
    else:
        certain_spots = type_df[type_df["spot_type"] != "uncertain"].index
    coordinate_df = coordinate_df.reindex(index=certain_spots)
    type_df = type_df.reindex(index=certain_spots)

    ncs = len(set(type_df["type_1"]))
    if ncs <= 10:
        cmap = "tab10"
        ncol = 1
    elif ncs <= 20:
        cmap = "tab20"
        ncol = 2
    else:
        cmap = color_palette(palette="hls", n_colors=ncs, as_cmap=True)
        ncol = 4

    if mpl.rcParams["legend.title_fontsize"] is None:
        legend_fontsize = 16
    else:
        legend_fontsize = mpl.rcParams["legend.title_fontsize"]

    fig, ax = plt.subplots(figsize=(10, 10), dpi=dpi)
    sc = ax.scatter(
        coordinate_df["X"],
        coordinate_df["Y"],
        s=s,
        c=type_df["type_1"],
        cmap=cmap,
    )
    leg = ax.legend(
        *sc.legend_elements(),
        ncol=ncol,
        title="Cluster",
        title_fontsize=legend_fontsize,
        fontsize=legend_fontsize,
    )
    ax.add_artist(leg)

    if he_image is not None:
        with Image.open(he_image) as he_image:
            ax.imshow(he_image)

    ax.axis("off")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Spot type map", fontsize=legend_fontsize)

    fig.savefig(save_path)
    plt.close(fig)


def svg_heatmap(
    dataset: STDataset,
    save_path: Union[str, Path],
    he_image: Optional[Union[str, Path]] = None,
    s: float = 4,
    dpi: float = 300,
) -> None:
    """
    Draw SVG distribution heatmap.

    Parameters
    ==========
    dataset : STDataset
        A STDataset with hotspot and SVG cluster estimation finished.

    save_path : str or pathlib.Path
        Heatmap save path.

    he_image : str or pathlib.Path, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.
    """
    coor_df = dataset.coordinate_df
    if he_image is None:
        if dataset._array_coordinate is not None:
            coor_df = dataset._array_coordinate
    _svg_heatmap(
        hotspot_df=dataset.hotspot_df,
        coordinate_df=coor_df,
        cluster_result=dataset.svg_cluster,
        spot_type=dataset.spot_type,
        save_path=save_path,
        he_image=he_image,
        s=s,
        dpi=dpi,
    )


def hotspot_distribution_map(
    dataset: STDataset,
    cluster: Union[str, int],
    save_path: Union[str, Path],
    he_image: Optional[Union[str, Path]] = None,
    s: float = 4,
    dpi: float = 300,
) -> None:
    """
    Draw hotspot distribution map for one SVG cluster.

    Parameters
    ==========
    dataset : STDataset
        A STDataset with hotspot and SVG cluster estimation finished.

    cluster : str or int
        Specify drawing cluster.

    save_path : str or pathlib.Path
        Heatmap save path.

    he_image : str or pathlib.Path, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.
    """
    coor_df = dataset.coordinate_df
    if he_image is None:
        if dataset._array_coordinate is not None:
            coor_df = dataset._array_coordinate
    _hotspot_distribution_map(
        hotspot_df=dataset.hotspot_df,
        coordinate_df=coor_df,
        cluster_result=dataset.svg_cluster,
        cluster=cluster,
        save_path=save_path,
        he_image=he_image,
        s=s,
        dpi=dpi,
    )


def spot_type_map(
    dataset: STDataset,
    save_path: Union[str, Path],
    he_image: Optional[Union[str, Path]] = None,
    draw_uncertain: bool = True,
    s: float = 16,
    dpi: float = 300,
) -> None:
    """
    Draw SVG type map.

    Parameters
    ==========
    dataset : STDataset
        A STDataset with hotspot and SVG cluster estimation finished.

    save_path : str or pathlib.Path
        Type map save path.

    he_image : str or pathlib.Path, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    draw_uncertain : bool, default True
        Whether to draw uncertain spots.

    s : float, default 16
        Spot size.

    dpi : float, default 300
        DPI for saved figure.
    """
    coor_df = dataset.coordinate_df
    if he_image is None:
        if dataset._array_coordinate is not None:
            coor_df = dataset._array_coordinate
    _spot_type_map(
        hotspot_df=dataset.hotspot_df,
        coordinate_df=coor_df,
        type_df=dataset.spot_type,
        save_path=save_path,
        he_image=he_image,
        draw_uncertain=draw_uncertain,
        s=s,
        dpi=dpi,
    )
