from typing import Tuple, Union, Optional
from pathlib import Path

import numpy as np
import pandas as pd
from libpysal.weights import KNN
from libpysal.weights import W as libpysal_W

from . import cluster, density, moran

DataFrames = Union[pd.DataFrame, np.ndarray, Path, str]


class STDataset(object):
    def __init__(
        self,
        count_df: DataFrames,
        coordinate_df: DataFrames,
        count_transpose: bool = False,
        coordinate_transpose: bool = False,
        count_df_kwargs: dict = None,
        coordinate_df_kwargs: dict = None,
    ) -> None:
        """
        STDataset
        """

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
        self._gene_cluster: Optional[pd.DataFrame] = None

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

        err = "Spots' name mismatch!"
        assert all(self._count_df.index == self._coordinate_df.index), err

    def acquire_weight(self, k: int = 6, **kwargs) -> None:
        self._weight = KNN(self._coordinate_df, k=k, **kwargs)
        self._weight_type = ("KNN", str(k))

    def acquire_hotspot(self, **kwargs) -> None:
        if self._weight is None:
            self.acquire_weight()
        result_dict = moran.local_moran(
            gene_expression_df=self.count_df,
            weights=self.weight,
            **kwargs,
        )
        self._hotspot_df = result_dict["hotspot"].reindex(
            index=self.spots,
            columns=self.genes,
        )
        self._local_moran_i = result_dict["i_value"]
        self._local_moran_p = result_dict["p_value"]

    def acquire_density(self, cores: int = density.cpu_count()) -> None:
        if self._hotspot_df is None:
            self.acquire_hotspot()
        result_dict = density.hotspot_AI(
            hotspot_df=self._hotspot_df,
            weight_df=self._local_moran_p,
            knn=self._weight,
        )
        self._AI = result_dict["AI"]
        self._Di = result_dict["Di"]

    def find_clusters(
        self,
        n_genes: int = 1000,
        n_gene_clusters: int = 8,
    ) -> None:
        self._gene_cluster = cluster.cluster(
            self._hotspot_df,
            self._AI,
            n_genes=n_genes,
            n_gene_clusters=n_gene_clusters,
        )

    @property
    def count_df(self) -> pd.DataFrame:
        return self._count_df

    @property
    def coordinate_df(self) -> pd.DataFrame:
        return self._coordinate_df

    @property
    def n_spots(self) -> int:
        return self._count_df.shape[0]

    @property
    def spots(self) -> pd.Series:
        return self._count_df.index

    @property
    def n_genes(self) -> int:
        return self._count_df.shape[1]

    @property
    def genes(self) -> pd.Series:
        return self._count_df.columns

    @property
    def weight(self) -> libpysal_W:
        return self._weight

    @weight.setter
    def weight(self, value: libpysal_W):
        self._weight = value
        self._weight_type = ("User specified weight", None)

    @property
    def weight_type(self) -> Tuple[Optional[str], Optional[str]]:
        return self._weight_type

    @property
    def hotspot_df(self) -> pd.DataFrame:
        return self._hotspot_df

    @property
    def AI(self) -> pd.Series:
        return self._AI

    @property
    def Di(self) -> pd.DataFrame:
        return self._Di

    @property
    def gene_cluster(self) -> pd.Series:
        return self._gene_cluster


if __name__ == "__main__":
    pass
