from functools import partial
from multiprocessing import Pool, cpu_count
from typing import Dict

import pandas as pd
from libpysal.weights import W as libpysal_W
from pysal.explore import esda


def local_moran(
    gene_expression_df: pd.DataFrame,
    weights: libpysal_W,
    transform_moran: str = 'r',
    permutation: int = 999,
    cores: int = cpu_count(),
) -> Dict[str, pd.DataFrame]:
    partial_func = partial(
        _local_moran,
        gene_expression_df=gene_expression_df,
        weights=weights,
        transform_moran=transform_moran,
        permutation=permutation,
    )
    pool = Pool(processes=cores)
    result_lists = pool.map(partial_func, gene_expression_df.columns)
    pool.close()
    pool.join()
    hotspot_df = pd.concat([i["hotspot"] for i in result_lists], axis=1)
    i_value_df = pd.concat([i["i_value"] for i in result_lists], axis=1)
    p_value_df = pd.concat([i["p_value"] for i in result_lists], axis=1)

    return {
        "hotspot": hotspot_df,
        "i_value": i_value_df,
        "p_value": p_value_df,
    }


def _local_moran(
    gene: str,
    gene_expression_df: pd.DataFrame,
    weights: libpysal_W,
    transform_moran: str = "r",
    permutation: int = 999,
) -> pd.DataFrame:
    indexs = gene_expression_df.index.tolist()
    gene_lisa = esda.moran.Moran_Local(
        gene_expression_df[gene],
        weights,
        transformation=transform_moran,
        permutations=permutation,
    )
    gene_moran_df = pd.DataFrame(0, index=indexs, columns=[gene])
    # select the significant hotspots (p < 0.05)
    local_moran_hotspot = 1 * ((gene_lisa.p_sim < 0.05) & (gene_lisa.q == 1))
    gene_moran_df[gene] = local_moran_hotspot
    if (gene_moran_df[gene] == 1).all():
        gene_moran_df[gene] = 0

    gene_i_df = pd.DataFrame(0, index=indexs, columns=[gene])
    gene_i_df[gene] = gene_lisa.EI_sim

    gene_p_df = pd.DataFrame(0, index=indexs, columns=[gene])
    gene_p_df[gene] = gene_lisa.p_sim

    gene_z_sim_df = pd.DataFrame(0, index=indexs, columns=[gene])
    gene_z_sim_df[gene] = gene_lisa.z_sim

    return {
        "hotspot": gene_moran_df,
        "i_value": gene_i_df,
        "p_value": gene_p_df,
        "z_sim": gene_z_sim_df,
    }


def global_moran(
    gene_expression_df: pd.DataFrame,
    weights: libpysal_W,
    transform_moran: str = 'r',
    permutation: int = 999,
    cores: int = cpu_count(),
) -> pd.DataFrame:
    partial_func = partial(
        _global_moran,
        gene_expression_df=gene_expression_df,
        weights=weights,
        transform_moran=transform_moran,
        permutation=permutation,
    )
    pool = Pool(processes=cores)
    result_lists = pool.map(partial_func, gene_expression_df.columns)
    pool.close()
    pool.join()
    result_df = pd.concat(result_lists, axis=0)

    global_genes_p_df = result_df
    global_genes_p_df.sort_values(
        by="p_value_sim",
        inplace=False,
        ascending=True,
    )

    return global_genes_p_df


def _global_moran(
    gene: str,
    gene_expression_df: pd.DataFrame,
    weights: libpysal_W,
    transform_moran: str = 'r',
    permutation: int = 999,
) -> pd.DataFrame:
    global_gene_p_df = pd.DataFrame(
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
        transformation=transform_moran,
        permutations=permutation,
    )
    p_value_sim = moran_value.p_sim
    moranI = moran_value.I
    p_value_norm = moran_value.p_norm
    p_value_rand = moran_value.p_rand
    p_value_z_sim = moran_value.p_z_sim
    global_gene_p_df['p_value_sim'][gene] = p_value_sim
    global_gene_p_df['I_value'][gene] = moranI
    global_gene_p_df['p_value_rand'][gene] = p_value_rand
    global_gene_p_df['p_value_z_sim'] = p_value_z_sim
    global_gene_p_df['p_value_norm'] = p_value_norm

    return global_gene_p_df


if __name__ == "__main__":
    pass
