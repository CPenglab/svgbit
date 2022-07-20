from __future__ import annotations
from typing import Tuple, Union, Optional
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from libpysal.weights import KNN
from libpysal.weights import W as libpysal_W

from . import cluster, density, moran, plot

DataFrames = Union[pd.DataFrame, np.ndarray, Path, str]


class STDataset(object):
    """
    STDataset: A meta class for discribing Spatial Transcriptomics data.

    Parameters
    ==========
    count_df : np.ndarray, pd.DataFrame, str or Path
       Expression matrix for Spatial Transcriptomics Data. If ``str`` or ``Path``
       is given, svgbit will try to read file with given path with pandas.

       Default shape: (spot * gene)

    coordinate_df : np.ndarray, pd.DataFrame, str or Path
       Coordinates for Spatial Transcriptomics Data. If ``str`` or ``Path``
       is given, svgbit will try to read file with given path with pandas.

       Default shape: (spot * 2)

    count_transpose : bool, default False
        Whether to transpose count matrix.

    coordinate_transpose : bool, default False
        Whether to transpose coordinate dataframe.

    count_df_kwargs : dict, default {}
        Keyword arguments pass to ``pandas.read_csv`` if ``str`` or ``Path`` is
        given to ``count_df``.

    coordinate_df_kwargs : dict, default {}
        Keyword arguments pass to ``pandas.read_csv`` if ``str`` or ``Path`` is
        given to ``coordinate_df``.

    """
    def __init__(
        self,
        count_df: DataFrames,
        coordinate_df: DataFrames,
        count_transpose: bool = False,
        coordinate_transpose: bool = False,
        count_df_kwargs: dict = {},
        coordinate_df_kwargs: dict = {},
    ) -> None:

        # attributes initial
        self._count_df: Optional[pd.DataFrame] = None
        self._coordinate_df: Optional[pd.DataFrame] = None
        self._weight: Optional[libpysal_W] = None
        self._weight_type: Tuple[Optional[str], Optional[str]] = (None, None)
        self._hotspot_df: Optional[pd.DataFrame] = None
        self._local_moran_i: Optional[pd.DataFrame] = None
        self._local_moran_p: Optional[pd.DataFrame] = None
        self._AI: Optional[pd.DataFrame] = None
        self._Di: Optional[pd.DataFrame] = None
        self._svg_cluster: Optional[pd.Series] = None

        # dataframes check
        if isinstance(count_df, pd.DataFrame):
            self._count_df = count_df
        elif isinstance(count_df, np.ndarray):
            self._count_df = pd.DataFrame(count_df)
        else:
            self._count_df = pd.read_csv(count_df, **count_df_kwargs)
        if count_transpose:
            self._count_df = self._count_df.T

        if isinstance(coordinate_df, pd.DataFrame):
            self._coordinate_df = coordinate_df
        elif isinstance(coordinate_df, np.ndarray):
            self._coordinate_df = pd.DataFrame(coordinate_df)
        else:
            self._coordinate_df = pd.read_csv(
                coordinate_df,
                **coordinate_df_kwargs,
            )
        if coordinate_transpose:
            self._coordinate_df = self._coordinate_df.T

        self._count_df.sort_index(inplace=True)
        self._coordinate_df.sort_index(inplace=True)
        self._coordinate_df.columns = ["X", "Y"]

        err = "Expression matrix and coordinate file have different number of spots."
        assert self._count_df.shape[0] == self._coordinate_df.shape[0], err
        err = "Spots' name mismatch!"
        assert all(self._count_df.index == self._coordinate_df.index), err

    def acquire_weight(self, k: int = 6, **kwargs) -> None:
        """
        Acquire weight for analysis.

        Parameters
        ==========
        k: int, default 6
            Number of nearest neighbors for KNN network.

        **kwargs
            Additional keyword arguments passed to the libpysal.weights.KNN call.

        """
        self._weight = KNN(self._coordinate_df, k=k, **kwargs)
        self._weight_type = ("KNN", str(k))

    def acquire_hotspot(self, **kwargs) -> None:
        """
        Acquire hotspot matrix.

        Parameters
        ==========
        **kwargs
            Additional keyword arguments passed to local_moran call.

        """
        if self._weight is None:
            self.acquire_weight()
        hotspot, i_value, p_value = moran.local_moran(
            gene_expression_df=self.count_df,
            weights=self.weight,
            **kwargs,
        )
        self._hotspot_df = hotspot.reindex(
            index=self.spots,
            columns=self.genes,
        )
        self._local_moran_i = i_value
        self._local_moran_p = p_value

    def acquire_density(self, cores: int = density.cpu_count()) -> None:
        """
        Acquire local Di and global AI value.

        Parameters
        ==========
        cores: int
            Number of threads to run svgbit. Use all available cpus by default.

        """
        if self._hotspot_df is None:
            self.acquire_hotspot()
        self._AI, self._Di = density.hotspot_AI(
            hotspot_df=self._hotspot_df,
            weight_df=self._local_moran_p,
            knn=self._weight,
            cores=cores,
        )

    def find_clusters(
        self,
        n_svgs: int = 1000,
        n_svg_clusters: int = 8,
    ) -> None:
        """
        Find SVG clusters.

        Parameters
        ==========
        n_svgs: int, default 1000
            Number of SVGs to find clusters.

        n_svg_clusters: int, default 8
            Number of SVG clusters to find.

        """
        self._svg_cluster = cluster.cluster(
            self._hotspot_df,
            self._AI,
            n_svgs=n_svgs,
            n_svg_clusters=n_svg_clusters,
        )

    def svg_heatmap(
        self,
        save_path: Union[str, Path],
        he_image=None,
    ) -> Tuple[mpl.figure.Figure, np.ndarray[mpl.axes.Axes]]:
        """
        Draw SVG distribution heatmap.

        Patameters
        ==========
        save_path: str or pathlib.Path
            Heatmap save path.

        he_image: PIL.Image.Image, default None
            H&E image of tissue. If None is given (default), distribution map
            will not show tissue picture.
        """
        plot.svg_heatmap(
            self._hotspot_df,
            self._coordinate_df,
            self._svg_cluster,
            save_path,
            he_image,
        )

    @property
    def count_df(self) -> pd.DataFrame:
        """Expression matrix."""
        return self._count_df

    @property
    def coordinate_df(self) -> pd.DataFrame:
        """Coordinate information."""
        return self._coordinate_df

    @property
    def n_spots(self) -> int:
        """Number of total spots."""
        return self._count_df.shape[0]

    @property
    def spots(self) -> pd.Index:
        """An Index for spots' names."""
        return self._count_df.index

    @property
    def n_genes(self) -> int:
        """Number of total genes."""
        return self._count_df.shape[1]

    @property
    def genes(self) -> pd.Index:
        """An Index for genes' names."""
        return self._count_df.columns

    @property
    def weight(self) -> libpysal_W:
        """Weight used by svgbit. Use KNN if not specified."""
        return self._weight

    @weight.setter
    def weight(self, value: libpysal_W):
        self._weight = value
        self._weight_type = ("User specified weight", None)

    @property
    def weight_type(self) -> Tuple[Optional[str], Optional[str]]:
        """
        What kind of weight is used. The second element indicates parameter k
        used by KNN in default.
        """
        return self._weight_type

    @property
    def hotspot_df(self) -> pd.DataFrame:
        """Hotspot matrix."""
        return self._hotspot_df

    @property
    def AI(self) -> pd.Series:
        """A Series for AI value."""
        return self._AI

    @property
    def Di(self) -> pd.DataFrame:
        """A DataFrame for local Di value."""
        return self._Di

    @property
    def svg_cluster(self) -> pd.Series:
        """SVG cluster result."""
        return self._svg_cluster


if __name__ == "__main__":
    pass
