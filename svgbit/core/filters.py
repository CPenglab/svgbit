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
        A STDataset instance with filtered genes.
    """
    var_series = dataset.count_df.var()
    var_series = var_series[var_series > var]
    return STDataset(
        dataset.count_df.reindex(columns=var_series.index),
        dataset.coordinate_df,
    )


def high_expression_filter(
    dataset: STDataset,
    max_ratio: float = 0.95,
) -> STDataset:
    temp_df = dataset.count_df.where(dataset.count_df < 1, 1)
    temp_series = temp_df.sum() / temp_df.shape[0]
    drop_genes = []
    for i in temp_series.index:
        if temp_series[i] > max_ratio:
            drop_genes.append(i)
    count_df = dataset.count_df.drop(columns=drop_genes)
    return STDataset(count_df, dataset.coordinate_df)