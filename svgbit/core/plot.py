from __future__ import annotations

from itertools import combinations
from math import ceil
from pathlib import Path
from typing import Optional, Union, List

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
from scipy.special import comb

from .STDataset import STDataset
from .utils import get_cmap


def _svg_heatmap(
    hotspot_df: pd.DataFrame,
    coordinate_df: pd.DataFrame,
    cluster_result: pd.Series,
    spot_type: pd.DataFrame,
    save_path: Optional[Union[str, Path]] = None,
    he_image: Optional[Union[str, Path, Image.Image]] = None,
    s: float = 4,
    dpi: float = 300,
) -> mpl.figure.Figure:
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

    save_path : str or pathlib.Path, default None
        Heatmap save path. If None, plot will not save.

    he_image : str, pathlib.Path or PIL.Image.Image, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.

    Returns
    =======
    fig : matplotlib.figure.Figure
        A matplotlib.figure.Figure plot.
    """
    spot_type = spot_type.sort_values(by="type_1")

    left, bottom = 0.1, 0.1
    width, height = 0.66, 0.66
    spacing, cluster_width = 0.03, 0.22
    rect_heatmap = [left + spacing, bottom + spacing * 2, width, height]

    he_close = False
    figsize = (10, 10)
    if he_image is not None:
        if not isinstance(he_image, Image.Image):
            he_close = True
            he_image = Image.open(he_image)
            figsize = (10, he_image.height / he_image.width * 10)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax_heatmap = fig.add_axes(rect_heatmap)

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

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight")

    he_image.close() if he_close else None

    return fig


def _hotspot_distribution_map(
    hotspot_df: pd.DataFrame,
    coordinate_df: pd.DataFrame,
    cluster_result: pd.Series,
    cluster: Union[str, int],
    save_path: Optional[Union[str, Path]] = None,
    he_image: Optional[Union[str, Path, Image.Image]] = None,
    s: float = 4,
    dpi: float = 300,
) -> mpl.figure.Figure:
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

    save_path : str or pathlib.Path, default None
        Heatmap save path. If None, plot will not save.

    he_image : str, pathlib.Path or PIL.Image.Image, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.

    Returns
    =======
    fig : matplotlib.figure.Figure
        A matplotlib.figure.Figure plot.
    """
    he_close = False
    figsize = (10, 10)
    if he_image is not None:
        if not isinstance(he_image, Image.Image):
            he_close = True
            he_image = Image.open(he_image)
            figsize = (10, he_image.height / he_image.width * 10)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

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

    ax.imshow(he_image) if he_image is not None else None

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight")

    he_image.close() if he_close else None

    return fig


def _hotspot_expression(
    hotspot_df: pd.DataFrame,
    coordinate_df: pd.DataFrame,
    gene: str,
    save_path: Optional[Union[str, Path]] = None,
    he_image: Optional[Union[str, Path, Image.Image]] = None,
    s: float = 4,
    dpi: float = 300,
) -> mpl.figure.Figure:
    """
    Draw hotspot expression for one gene.

    Parameters
    ==========
    hotspot_df : pd.DataFrame
        A pd.DataFrame for hotspots.

    coordinate_df : pd.DataFrame
        A pd.DataFrame for coordinate files.

    gene : str
        Which gene to draw.

    save_path : str or pathlib.Path, default None
        Heatmap save path. If None, plot will not save.

    he_image : str, pathlib.Path or PIL.Image.Image, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.

    Returns
    =======
    fig : matplotlib.figure.Figure
        A matplotlib.figure.Figure plot.
    """
    he_close = False
    figsize = (10, 10)
    if he_image is not None:
        if not isinstance(he_image, Image.Image):
            he_close = True
            he_image = Image.open(he_image)
            figsize = (10, he_image.height / he_image.width * 10)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    ax.set_title(gene)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")
    sc = ax.scatter(
        coordinate_df.iloc[:, 0],
        coordinate_df.iloc[:, 1],
        c=hotspot_df[gene],
        cmap="autumn_r",
        vmin=0,
        vmax=1,
        s=s,
        alpha=0.7,
    )
    fig.colorbar(sc)

    ax.imshow(he_image) if he_image is not None else None

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight")

    he_image.close() if he_close else None

    return fig


def _hotspot_colocalization_map(
    hotspot_df: pd.DataFrame,
    coordinate_df: pd.DataFrame,
    genes: list,
    save_path: Optional[Union[str, Path]] = None,
    he_image: Optional[Union[str, Path, Image.Image]] = None,
    s: float = 4,
    dpi: float = 300,
    colors: Optional[List] = None,
) -> mpl.figure.Figure:
    """
    Draw hotspot expression for multiple genes.

    Parameters
    ==========
    hotspot_df : pd.DataFrame
        A pd.DataFrame for hotspots.

    coordinate_df : pd.DataFrame
        A pd.DataFrame for coordinate files.

    genes : str
        Which genes to draw.

    save_path : str or pathlib.Path, default None
        Heatmap save path. If None, plot will not save.

    he_image : str, pathlib.Path or PIL.Image.Image, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.

    colors : list, default None
        A list for spot colors. If None, auto generate colors.

    Returns
    =======
    fig : matplotlib.figure.Figure
        A matplotlib.figure.Figure plot.
    """
    he_close = False
    figsize = (10, 10)
    if he_image is not None:
        if not isinstance(he_image, Image.Image):
            he_close = True
            he_image = Image.open(he_image)
            figsize = (10, he_image.height / he_image.width * 10)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    if colors is None:
        if len(genes) <= 3:
            colors = [
                "tab:red", "tab:green", "tab:blue", "tab:orange", "tab:pink",
                "tab:cyan", "k"
            ]
        else:
            i = len(genes)
            n_colors = 0
            while i > 0:
                n_colors += comb(len(genes), i)
                i -= 1
            colors = get_cmap(n_colors)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")

    for i, gene in enumerate(genes):
        draw_counts = hotspot_df[gene].loc[hotspot_df[gene] == 1]
        ax.scatter(
            coordinate_df.iloc[:, 0].reindex(index=draw_counts.index),
            coordinate_df.iloc[:, 1].reindex(index=draw_counts.index),
            s=s,
            color=colors[i],
            alpha=0.7,
            label=f"{gene}",
        )

    n_combinations = 2
    while n_combinations <= len(genes):
        for c in combinations(genes, n_combinations):
            i += 1
            draw_counts = pd.Series(index=coordinate_df.index, dtype=int)
            draw_counts = draw_counts.fillna(0)
            for gene in c:
                draw_counts += hotspot_df[gene]
            draw_counts = draw_counts[draw_counts == n_combinations]
            ax.scatter(
                coordinate_df.iloc[:, 0].reindex(index=draw_counts.index),
                coordinate_df.iloc[:, 1].reindex(index=draw_counts.index),
                c=colors[i],
                s=s,
                label=" & ".join(c),
            )
        n_combinations += 1

    ax.legend(markerscale=2)

    ax.imshow(he_image) if he_image is not None else None

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight")

    he_image.close() if he_close else None

    return fig


def _gene_expression(
    expression_df: pd.DataFrame,
    coordinate_df: pd.DataFrame,
    gene: str,
    save_path: Optional[Union[str, Path]] = None,
    he_image: Optional[Union[str, Path, Image.Image]] = None,
    s: float = 4,
    dpi: float = 300,
) -> mpl.figure.Figure:
    """
    Draw expression for one gene.

    Parameters
    ==========
    count_df : pd.DataFrame
        A pd.DataFrame for hotspots.

    coordinate_df : pd.DataFrame
        A pd.DataFrame for coordinate files.

    gene : str
        Which gene to draw.

    save_path : str or pathlib.Path, default None
        Heatmap save path. If None, plot will not save.

    he_image : str, pathlib.Path or PIL.Image.Image, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.

    Returns
    =======
    fig : matplotlib.figure.Figure
        A matplotlib.figure.Figure plot.
    """
    he_close = False
    figsize = (10, 10)
    if he_image is not None:
        if not isinstance(he_image, Image.Image):
            he_close = True
            he_image = Image.open(he_image)
            figsize = (10, he_image.height / he_image.width * 10)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    ax.set_title(gene)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")
    sc = ax.scatter(
        coordinate_df.iloc[:, 0],
        coordinate_df.iloc[:, 1],
        c=expression_df[gene],
        cmap="autumn_r",
        s=s,
        alpha=0.7,
    )
    fig.colorbar(sc)

    ax.imshow(he_image) if he_image is not None else None

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight")

    he_image.close() if he_close else None

    return fig


def _spot_type_map(
    hotspot_df: pd.DataFrame,
    coordinate_df: pd.DataFrame,
    type_df: pd.DataFrame,
    save_path: Optional[Union[str, Path]] = None,
    he_image: Optional[Union[str, Path, Image.Image]] = None,
    draw_uncertain: bool = True,
    s: float = 16,
    dpi: float = 300,
) -> mpl.figure.Figure:
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

    save_path : str or pathlib.Path, default None
        Heatmap save path. If None, plot will not save.

    he_image : str, pathlib.Path or PIL.Image.Image, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    draw_uncertain : bool, default True
        Whether to draw uncertain spots.

    s : float, default 16
        Spot size.

    dpi : float, default 300
        DPI for saved figure.

    Returns
    =======
    fig : matplotlib.figure.Figure
        A matplotlib.figure.Figure plot.
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
        cmap = get_cmap(ncs)
        ncol = 4

    if mpl.rcParams["legend.title_fontsize"] is None:
        legend_fontsize = 16
    else:
        legend_fontsize = mpl.rcParams["legend.title_fontsize"]

    he_close = False
    figsize = (10, 10)
    if he_image is not None:
        if not isinstance(he_image, Image.Image):
            he_close = True
            he_image = Image.open(he_image)
            figsize = (10, he_image.height / he_image.width * 10)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

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

    ax.axis("off")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Spot type map", fontsize=legend_fontsize)

    ax.imshow(he_image) if he_image is not None else None

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight")

    he_image.close() if he_close else None

    return fig


def svg_heatmap(
    dataset: STDataset,
    save_path: Optional[Union[str, Path]] = None,
    he_image: Optional[Union[str, Path, Image.Image]] = None,
    s: float = 4,
    dpi: float = 300,
) -> mpl.figure.Figure:
    """
    Draw SVG distribution heatmap.

    Parameters
    ==========
    dataset : STDataset
        A STDataset with hotspot and SVG cluster estimation finished.

    save_path : str or pathlib.Path, default None
        Heatmap save path. If None, plot will not save.

    he_image : str, pathlib.Path or PIL.Image.Image, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.

    Returns
    =======
    fig : matplotlib.figure.Figure
        A matplotlib.figure.Figure plot.
    """
    coor_df = dataset.coordinate_df
    if he_image is None:
        if dataset._array_coordinate is not None:
            coor_df = dataset._array_coordinate
    return _svg_heatmap(
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
    save_path: Optional[Union[str, Path]] = None,
    he_image: Optional[Union[str, Path, Image.Image]] = None,
    s: float = 4,
    dpi: float = 300,
) -> mpl.figure.Figure:
    """
    Draw hotspot distribution map for one SVG cluster.

    Parameters
    ==========
    dataset : STDataset
        A STDataset with hotspot and SVG cluster estimation finished.

    cluster : str or int
        Specify drawing cluster.

    save_path : str or pathlib.Path, default None
        Heatmap save path. If None, plot will not save.

    he_image : str, pathlib.Path or PIL.Image.Image, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.

    Returns
    =======
    fig : matplotlib.figure.Figure
        A matplotlib.figure.Figure plot.
    """
    coor_df = dataset.coordinate_df
    if he_image is None:
        if dataset._array_coordinate is not None:
            coor_df = dataset._array_coordinate
    return _hotspot_distribution_map(
        hotspot_df=dataset.hotspot_df,
        coordinate_df=coor_df,
        cluster_result=dataset.svg_cluster,
        cluster=cluster,
        save_path=save_path,
        he_image=he_image,
        s=s,
        dpi=dpi,
    )


def hotspot_expression(
    dataset: STDataset,
    gene: str,
    save_path: Optional[Union[str, Path]] = None,
    he_image: Optional[Union[str, Path, Image.Image]] = None,
    s: float = 4,
    dpi: float = 300,
) -> mpl.figure.Figure:
    """
    Draw hotspot expression for one gene.

    Parameters
    ==========
    dataset : STDataset
        A STDataset with hotspot estimation finished.

    gene : str
        Which gene to draw.

    save_path : str or pathlib.Path, default None
        Heatmap save path. If None, plot will not save.

    he_image : str, pathlib.Path or PIL.Image.Image, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.

    Returns
    =======
    fig : matplotlib.figure.Figure
        A matplotlib.figure.Figure plot.
    """
    coor_df = dataset.coordinate_df
    if he_image is None:
        if dataset._array_coordinate is not None:
            coor_df = dataset._array_coordinate
    return _hotspot_expression(
        hotspot_df=dataset.hotspot_df,
        coordinate_df=coor_df,
        gene=gene,
        save_path=save_path,
        he_image=he_image,
        s=s,
        dpi=dpi,
    )


def hotspot_colocalization_map(
    dataset: STDataset,
    genes: list,
    save_path: Optional[Union[str, Path]] = None,
    he_image: Optional[Union[str, Path, Image.Image]] = None,
    s: float = 4,
    dpi: float = 300,
    colors: Optional[list] = None,
) -> mpl.figure.Figure:
    """
    Draw hotspot expression for multiple genes.

    Parameters
    ==========
    dataset : STDataset
        A STDataset with hotspot estimation finished.

    genes : str
        Which genes to draw.

    save_path : str or pathlib.Path, default None
        Heatmap save path. If None, plot will not save.

    he_image : str, pathlib.Path or PIL.Image.Image, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.

    colors : list, default None
        A list for spot colors. If None, auto generate colors.

    Returns
    =======
    fig : matplotlib.figure.Figure
        A matplotlib.figure.Figure plot.
    """
    coor_df = dataset.coordinate_df
    if he_image is None:
        if dataset._array_coordinate is not None:
            coor_df = dataset._array_coordinate
    return _hotspot_colocalization_map(
        hotspot_df=dataset.hotspot_df,
        coordinate_df=coor_df,
        genes=genes,
        save_path=save_path,
        he_image=he_image,
        s=s,
        dpi=dpi,
        colors=colors,
    )


def gene_expression(
    dataset: STDataset,
    gene: str,
    save_path: Optional[Union[str, Path]] = None,
    he_image: Optional[Union[str, Path, Image.Image]] = None,
    s: float = 4,
    dpi: float = 300,
) -> mpl.figure.Figure:
    """
    Draw expression for one gene.

    Parameters
    ==========
    dataset : STDataset
        A STDataset with hotspot and SVG cluster estimation finished.

    gene : str
        Which gene to draw.

    save_path : str or pathlib.Path, default None
        Heatmap save path. If None, plot will not save.

    he_image : str, pathlib.Path or PIL.Image.Image, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    s : float, default 4
        Spot size.

    dpi : float, default 300
        DPI for saved figure.

    Returns
    =======
    fig : matplotlib.figure.Figure
        A matplotlib.figure.Figure plot.
    """
    coor_df = dataset.coordinate_df
    if he_image is None:
        if dataset._array_coordinate is not None:
            coor_df = dataset._array_coordinate
    return _gene_expression(
        expression_df=dataset.count_df,
        coordinate_df=coor_df,
        gene=gene,
        save_path=save_path,
        he_image=he_image,
        s=s,
        dpi=dpi,
    )


def spot_type_map(
    dataset: STDataset,
    save_path: Optional[Union[str, Path]] = None,
    he_image: Optional[Union[str, Path, Image.Image]] = None,
    draw_uncertain: bool = True,
    s: float = 16,
    dpi: float = 300,
) -> mpl.figure.Figure:
    """
    Draw SVG type map.

    Parameters
    ==========
    dataset : STDataset
        A STDataset with hotspot and SVG cluster estimation finished.

    save_path : str or pathlib.Path, default None
        Heatmap save path. If None, plot will not save.

    he_image : str, pathlib.Path or PIL.Image.Image, default None
        H&E image of tissue. If None is given (default), distribution map
        will not show tissue picture.

    draw_uncertain : bool, default True
        Whether to draw uncertain spots.

    s : float, default 16
        Spot size.

    dpi : float, default 300
        DPI for saved figure.

    Returns
    =======
    fig : matplotlib.figure.Figure
        A matplotlib.figure.Figure plot.
    """
    coor_df = dataset.coordinate_df
    if he_image is None:
        if dataset._array_coordinate is not None:
            coor_df = dataset._array_coordinate
    return _spot_type_map(
        hotspot_df=dataset.hotspot_df,
        coordinate_df=coor_df,
        type_df=dataset.spot_type,
        save_path=save_path,
        he_image=he_image,
        draw_uncertain=draw_uncertain,
        s=s,
        dpi=dpi,
    )
