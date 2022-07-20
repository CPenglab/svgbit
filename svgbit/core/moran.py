from functools import partial
from multiprocessing import Pool, cpu_count
from typing import Tuple

import numba
import pandas as pd
from libpysal.weights import W as libpysal_W
from pysal.explore import esda


def local_moran(
        gene_expression_df: pd.DataFrame,
        weights: libpysal_W,
        transformation: str = 'r',
        permutation: int = 999,
        cores: int = cpu_count(),
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Calculate local moran for hotspot identification.

    Parameters
    ==========
    gene_expression_df : pd.DataFrame
       Expression matrix for Spatial Transcriptomics Data.

    weights : libpysal.weights.W
        Spatial weight for calculating Moran's I.

    transformation : str, default 'r'
        Weights transformation passed to esda.moran.Moran_Local

    permutation : int, default 999
        Number of random permutations for calculation of pseudo-p_values.

    cores : int
        Number of threads to run svgbit. Use all available cpus by default.

    Returns
    =======
    hotspot : pd.DataFrame
        Identified hotspots.

    i_value : pd.DataFrame
        Local Moran's I values.

    p_value : pd.DataFrame
        Local Moran's I p values.

    """
    partial_func = partial(
        _local_moran,
        gene_expression_df=gene_expression_df,
        weights=weights,
        transformation=transformation,
        permutation=permutation,
    )
    pool = Pool(processes=cores)
    result_lists = pool.map(partial_func, gene_expression_df.columns)
    pool.close()
    pool.join()
    hotspot = pd.concat([i[0] for i in result_lists], axis=1)
    i_value = pd.concat([i[1] for i in result_lists], axis=1)
    p_value = pd.concat([i[2] for i in result_lists], axis=1)

    return hotspot, i_value, p_value


def _local_moran(
    gene: str,
    gene_expression_df: pd.DataFrame,
    weights: libpysal_W,
    transformation: str = "r",
    permutation: int = 999,
) -> Tuple[pd.Series, pd.Series, pd.Series]:
    """
    Calculate local moran for hotspot identification for a single gene.

    Parameters
    ==========
    gene : str
        Calculate which gene.

    gene_expression_df : pd.DataFrame
       Expression matrix for Spatial Transcriptomics Data.

    weights : libpysal.weights.W
        Spatial weight for calculating Moran's I.

    transformation : str, default 'r'
        Weights transformation passed to esda.moran.Moran_Local

    permutation : int, default 999
        Number of random permutations for calculation of pseudo-p_values.

    Returns
    =======
    hotspot : pd.DataFrame
        Identified hotspots.

    i_value : pd.DataFrame
        Local Moran's I values.

    p_value : pd.DataFrame
        Local Moran's I p values.

    """
    numba.config.THREADING_LAYER = "workqueue"
    gene_lisa = esda.moran.Moran_Local(
        gene_expression_df[gene],
        weights,
        transformation=transformation,
        permutations=permutation,
        n_jobs=1,
    )
    # select the significant hotspots (p < 0.05)
    local_moran_hotspot = 1 * ((gene_lisa.p_sim < 0.05) & (gene_lisa.q == 1))

    hotspot = pd.Series(
        local_moran_hotspot,
        index=gene_expression_df.index,
        name=gene,
    )
    if (hotspot == 1).all():
        hotspot = 0

    i_value = pd.Series(
        gene_lisa.EI_sim,
        index=gene_expression_df.index,
        name=gene,
    )
    p_value = pd.Series(
        gene_lisa.p_sim,
        index=gene_expression_df.index,
        name=gene,
    )

    return hotspot, i_value, p_value


def global_moran(
        gene_expression_df: pd.DataFrame,
        weights: libpysal_W,
        transformation: str = 'r',
        permutation: int = 999,
        cores: int = cpu_count(),
) -> pd.DataFrame:
    """
    Calculate global Moran's I value.

    Parameters
    ==========
    gene_expression_df : pd.DataFrame
       Expression matrix for Spatial Transcriptomics Data.

    weights : libpysal.weights.W
        Spatial weight for calculating Moran's I.

    transformation : str, default 'r'
        Weights transformation passed to esda.moran.Moran

    permutation : int, default 999
        Number of random permutations for calculation of pseudo-p_values.

    cores : int
        Number of threads to run svgbit. Use all available cpus by default.

    Returns
    =======
    global_moran_df : pd.DataFrame
        Global Moran's I result.

    """
    partial_func = partial(
        _global_moran,
        gene_expression_df=gene_expression_df,
        weights=weights,
        transformation=transformation,
        permutation=permutation,
    )
    pool = Pool(processes=cores)
    result_lists = pool.map(partial_func, gene_expression_df.columns)
    pool.close()
    pool.join()
    result_df = pd.concat(result_lists, axis=0)

    global_moran_df = result_df
    global_moran_df.sort_values(
        by="p_value_sim",
        inplace=False,
        ascending=True,
    )

    return global_moran_df


def _global_moran(
    gene: str,
    gene_expression_df: pd.DataFrame,
    weights: libpysal_W,
    transformation: str = 'r',
    permutation: int = 999,
) -> pd.DataFrame:
    """
    Calculate global Moran's I value for a single gene.

    Parameters
    ==========
    gene : str
        Calculate which gene.

    gene_expression_df : pd.DataFrame
       Expression matrix for Spatial Transcriptomics Data.

    weights : libpysal.weights.W
        Spatial weight for calculating Moran's I.

    transformation : str, default 'r'
        Weights transformation passed to esda.moran.Moran

    permutation : int, default 999
        Number of random permutations for calculation of pseudo-p_values.

    Returns
    =======
    global_moran_df : pd.DataFrame
        Global Moran's I result.

    """
    global_moran_df = pd.DataFrame(
        columns=[
            'I_value',
            'p_value_sim',
            'p_value_rand',
            'p_value_z_sim',
            'p_value_norm',
        ],
        index=[gene],
    )
    moran_value = esda.moran.Moran(
        gene_expression_df[gene],
        weights,
        transformation=transformation,
        permutations=permutation,
    )
    p_value_sim = moran_value.p_sim
    moranI = moran_value.I
    p_value_norm = moran_value.p_norm
    p_value_rand = moran_value.p_rand
    p_value_z_sim = moran_value.p_z_sim
    global_moran_df['p_value_sim'][gene] = p_value_sim
    global_moran_df['I_value'][gene] = moranI
    global_moran_df['p_value_rand'][gene] = p_value_rand
    global_moran_df['p_value_z_sim'] = p_value_z_sim
    global_moran_df['p_value_norm'] = p_value_norm

    return global_moran_df


if __name__ == "__main__":
    pass
