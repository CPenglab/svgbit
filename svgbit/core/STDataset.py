from __future__ import annotations

import warnings
from collections import Counter
from copy import deepcopy
from typing import Optional, Tuple, Union
from pathlib import Path

import numpy as np
import pandas as pd
from libpysal.weights import KNN
from libpysal.weights import W as libpysal_W

from . import cluster, density, moran

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

    make_sparse : bool, default True
        Whether to use sparse DataFrame in order to save memory.

    check_duplicate_genes : bool, default True
        Whether to check duplicated gene names.

    sort_spots : bool, default True
        Whether to sort spots with spots' name.
    """
    def __init__(
        self,
        count_df: DataFrames,
        coordinate_df: DataFrames,
        count_transpose: bool = False,
        coordinate_transpose: bool = False,
        count_df_kwargs: dict = {},
        coordinate_df_kwargs: dict = {},
        make_sparse: bool = True,
        check_duplicate_genes: bool = True,
        sort_spots: bool = True,
    ) -> None:

        # attributes initial
        self._count_df: Optional[pd.DataFrame] = None
        self._coordinate_df: Optional[pd.DataFrame] = None
        self._normalizer: Optional[str] = None
        self._weight: Optional[libpysal_W] = None
        self._weight_type: Tuple[Optional[str], Optional[str]] = (None, None)
        self._hotspot_df: Optional[pd.DataFrame] = None
        self._local_moran_i: Optional[pd.DataFrame] = None
        self._local_moran_p: Optional[pd.DataFrame] = None
        self._AI: Optional[pd.Series] = None
        self._Di: Optional[pd.DataFrame] = None
        self._svg_cluster: Optional[pd.Series] = None
        self._spot_type: Optional[pd.DataFrame] = None
        self._array_coordinate: Optional[pd.DataFrame] = None

        # dataframes check
        if isinstance(count_df, pd.DataFrame):
            self._count_df = deepcopy(count_df)
        elif isinstance(count_df, np.ndarray):
            self._count_df = pd.DataFrame(count_df)
        else:
            self._count_df = pd.read_csv(count_df, **count_df_kwargs)
        if count_transpose:
            self._count_df = self._count_df.T

        if isinstance(coordinate_df, pd.DataFrame):
            self._coordinate_df = deepcopy(coordinate_df)
        elif isinstance(coordinate_df, np.ndarray):
            self._coordinate_df = pd.DataFrame(coordinate_df)
        else:
            self._coordinate_df = pd.read_csv(
                coordinate_df,
                **coordinate_df_kwargs,
            )
        if coordinate_transpose:
            self._coordinate_df = self._coordinate_df.T

        if sort_spots:
            self._count_df.sort_index(inplace=True)
            self._coordinate_df = self._coordinate_df.reindex(
                index=self.count_df.index)
        self._coordinate_df.columns = ["X", "Y"]

        err = "Expression matrix and coordinate file have different number of spots."
        assert self._count_df.shape[0] == self._coordinate_df.shape[0], err
        err = "Spots' name mismatch!"
        assert all(self._count_df.index == self._coordinate_df.index), err

        self._count_df.fillna(0, inplace=True)

        # Rename duplicated columns
        if check_duplicate_genes:
            genes = []
            flag = 0
            c = Counter(self._count_df.columns)
            gene_suffix = {}
            for gene_name in self._count_df.columns:
                if c[gene_name] > 1:
                    flag = 1
                    if gene_name in gene_suffix:
                        gene_suffix[gene_name] += 1
                        gene_name += f".{gene_suffix[gene_name]}"
                    else:
                        gene_suffix[gene_name] = 0
                genes.append(gene_name)
            if flag:
                self._count_df.columns = genes
                print("Duplicated column names found. Auto rename.")
                warnings.warn("Duplicated column names found. Auto rename.")

        if make_sparse:
            self.to_sparse()

    def __repr__(self) -> str:
        descr = f"STDataset with n_spots x n_genes = {self.n_spots} x {self.n_genes}"
        descr = f"{descr}\nApplied normalizers: {self._normalizer}"
        descr = f"{descr}\nAssigned attributes: "
        flag = 0
        for attr in ["weight", "hotspot_df", "AI", "svg_cluster"]:
            if getattr(self, attr) is not None:
                descr += f"{attr}, "
                flag = 1
        if flag:
            descr = descr[:-2]
        return descr

    def __str__(self) -> str:
        return self.__repr__()

    def __getitem__(self, pos) -> STDataset:
        """
        Return a sub STDataset instance with empty attributes.
        """
        try:
            count_sub = self.count_df.loc[pos]
        except TypeError:
            count_sub = self.count_df.iloc[pos]
        coor_sub = self.coordinate_df.loc[count_sub.index, ]
        return STDataset(
            count_sub,
            coor_sub,
            check_duplicate_genes=False,
            sort_spots=False,
            make_sparse=pd.api.types.is_sparse(self.count_df.iloc[:, 1]),
        )

    def __del__(self) -> None:
        del self._count_df
        del self._coordinate_df
        del self._normalizer
        del self._hotspot_df
        del self._weight
        del self._weight_type
        del self._local_moran_i
        del self._local_moran_p
        del self._AI
        del self._Di
        del self._svg_cluster
        del self._spot_type
        del self._array_coordinate

    def to_dense(self) -> None:
        """Convert count_df with sparse values to dense."""
        self._count_df = self._count_df.sparse.to_dense()

    def to_sparse(self) -> None:
        """Convert count_df with dense values to sparse."""
        is_int = [i == "int" for i in self.count_df.dtypes]
        if all(is_int):
            dt = "int64"
        else:
            dt = "float64"
        self._count_df = self._count_df.astype(pd.SparseDtype(dt, 0))

    def acquire_weight(self, k: int = 6, **kwargs) -> None:
        """
        Acquire weight for analysis.

        Parameters
        ==========
        k : int, default 6
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
        self._local_moran_i = i_value.astype(pd.SparseDtype("float", 0))
        self._local_moran_p = p_value.astype(pd.SparseDtype("float", 0))

    def acquire_density(self, cores: int = density.cpu_count()) -> None:
        """
        Acquire local Di and global AI value.

        Parameters
        ==========
        cores : int
            Number of threads to run svgbit. Use all available cpus by default.

        """
        if self._hotspot_df is None:
            self.acquire_hotspot()
        results = density.hotspot_AI(
            hotspot_df=self._hotspot_df,
            weight_df=-np.log(self._local_moran_p),
            knn=self._weight,
            cores=cores,
        )
        self._AI = results[0].reindex(index=self.genes)
        self._Di = results[1].astype(pd.SparseDtype("float", 0)).reindex(
            index=self.spots, columns=self.genes)

    def find_clusters(
        self,
        n_svgs: int = 1000,
        n_svg_clusters: int = 8,
        threshold: float = 0.3,
    ) -> None:
        """
        Find SVG clusters.

        Parameters
        ==========
        n_svgs : int, default 1000
            Number of SVGs to find clusters.

        n_svg_clusters : int, default 8
            Number of SVG clusters to find.

        threshold : float, dafault 0.3
            min value to identify multiple svg clusters to spot.

        """
        results = cluster.cluster(
            self._hotspot_df,
            self._AI,
            n_svgs=n_svgs,
            n_svg_clusters=n_svg_clusters,
            threshold=threshold,
        )
        self._svg_cluster = results[0]
        self._spot_type = results[1].reindex(index=self.spots)

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

    @property
    def spot_type(self) -> pd.DataFrame:
        """A pd.DataFrame for spot type."""
        return self._spot_type


if __name__ == "__main__":
    pass
