from functools import partial
from multiprocessing import Pool, cpu_count
from typing import Dict, Union

import numpy as np
import pandas as pd
from libpysal.weights import W as pysal_W

from .utils import pysal_to_pandas


def _hotspot_AI(
    gene: str,
    hotspot_df: pd.DataFrame,
    weight_df: pd.DataFrame,
    knn: Union[pd.DataFrame, pysal_W],
) -> Dict[str, pd.Series]:

    if isinstance(knn, pysal_W):
        knn = pysal_to_pandas(knn)
        knn.index = hotspot_df.index
        knn.columns = hotspot_df.index

    gene_hotspot_df = hotspot_df[hotspot_df[gene] != 0]
    gene_hotspot_df = pd.DataFrame(gene_hotspot_df[gene], columns=[gene])
    n_hotspots = len(gene_hotspot_df)
    if n_hotspots == 0:
        ai_series = pd.Series([0], index=[gene], name="AI")
        di_series = pd.Series(
            [0] * len(hotspot_df),
            index=hotspot_df.index,
            name=gene,
        )
        return {"AI": ai_series, "Di": di_series}
    hotspots = gene_hotspot_df.index.tolist()
    hs = {}
    for i in hotspots:
        hi_wnn_df = pd.DataFrame(knn[i], columns=[i])
        hi_wnn_df = pd.DataFrame(hi_wnn_df[hi_wnn_df[i] == 1])
        knn_coors = hi_wnn_df.index.tolist()
        inter_spots = set(hotspots).intersection(set(knn_coors))
        hs[i] = ((weight_df[gene][inter_spots] /
                  weight_df[gene][knn_coors].sum()).sum())

    di_array = np.array(list(hs.values()))
    di_sum = di_array.sum()
    di_mean = di_sum / n_hotspots
    ai_series = pd.Series([di_mean], index=[gene], name="AI")
    di_series = pd.Series(
        hs,
        name=gene,
    ).reindex(index=hotspot_df.index).fillna(0)

    return {"AI": ai_series, "Di": di_series}


def hotspot_AI(
        hotspot_df: pd.DataFrame,
        weight_df: pd.DataFrame,
        knn: Union[pd.DataFrame, pysal_W],
        cores: int = cpu_count(),
) -> Dict[str, pd.DataFrame]:
    if isinstance(knn, pysal_W):
        knn = pysal_to_pandas(knn)
        knn.index = hotspot_df.index
        knn.columns = hotspot_df.index
    partial_func = partial(
        _hotspot_AI,
        hotspot_df=hotspot_df,
        weight_df=weight_df,
        knn=knn,
    )
    pool = Pool(processes=cores)
    result_lists = pool.map(partial_func, hotspot_df.columns)
    pool.close()
    pool.join()
    ai_df = pd.concat([i["AI"] for i in result_lists], axis=0)
    di_df = pd.concat([i["Di"] for i in result_lists], axis=1)

    return {"AI": ai_df, "Di": di_df}


if __name__ == "__main__":
    pass
