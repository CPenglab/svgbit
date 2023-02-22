from __future__ import annotations

import pandas as pd

from .STDataset import STDataset


def low_variance_filter(dataset: STDataset, var: float = 0) -> STDataset:
    """
    Filter genes with low variance. This may also filter genes with 0 expressions.

    Parameters
    ==========
    dataset : STDataset
        STDataset to be filtered.

    var : int, default 0
        Only remain genes with variance greater than (but not equal to) var.

    Returns
    =======
    dataset : STDataset
        A new STDataset instance with filtered genes.
    """
    var_series = pd.Series()
    for gene in dataset.count_df.columns:
        try:
            var_series[gene] = dataset.count_df[gene].sparse.to_dense().var()
        except AttributeError:
            var_series[gene] = dataset.count_df[gene].var()
    var_series = var_series[var_series > var]
    return_dset = STDataset(
        dataset.count_df.reindex(columns=var_series.index),
        dataset.coordinate_df,
        check_duplicate_genes=False,
        sort_spots=False,
    )
    return_dset._array_coordinate = dataset._array_coordinate
    return return_dset


def high_expression_filter(
    dataset: STDataset,
    max_ratio: float = 0.99,
) -> STDataset:
    """
    Filter genes with high expression ratio.

    Parameters
    ==========
    dataset : STDataset
        STDataset to be filtered.

    max_ratio : int, default 0.99
        Only remain genes with expression ratio less than (but not equal to)
        max_ratio.

    Returns
    =======
    dataset : STDataset
        A new STDataset instance with filtered genes.
    """
    temp_df = dataset.count_df.where(dataset.count_df < 1, 1)
    temp_series = temp_df.sum() / temp_df.shape[0]
    drop_genes = []
    for i in temp_series.index:
        if temp_series[i] > max_ratio:
            drop_genes.append(i)
    count_df = dataset.count_df.drop(columns=drop_genes)
    return_dset = STDataset(
        count_df,
        dataset.coordinate_df,
        check_duplicate_genes=False,
        sort_spots=False,
    )
    return_dset._array_coordinate = dataset._array_coordinate
    return return_dset


def quantile_filter(dataset: STDataset, quantile: float = 0.99) -> STDataset:
    """
    Filter genes with quantile.

    Parameters
    ==========
    dataset : STDataset
        STDataset to be filtered.

    quantile : int, default 0.99
        Only remain genes with mean less than (but not equal to) quantile number.

    Returns
    =======
    dataset : STDataset
        A new STDataset instance with filtered genes.
    """
    mean_series = dataset.count_df.mean()
    q = mean_series.quantile(quantile)
    remain_genes = mean_series[mean_series < q].index
    count_df = dataset.count_df.reindex(columns=remain_genes)
    return_dset = STDataset(
        count_df,
        dataset.coordinate_df,
        check_duplicate_genes=False,
        sort_spots=False,
    )
    return_dset._array_coordinate = dataset._array_coordinate
    return return_dset
