from __future__ import annotations

from numpy import log as nlog
from pandas import SparseDtype

from .STDataset import STDataset


def cpm_normalizer(dataset: STDataset) -> STDataset:
    """
    Perform CPM on dataset.

    Parameters
    ==========
    dataset : STDataset
        STDataset to be normalized.

    Returns
    =======
    dataset : STDataset
        A STDataset instance with normalized expression matrix.

    """
    count_df = dataset.count_df.sparse.to_dense()
    scale_df = (count_df.T * 10000 / count_df.T.sum())
    d = STDataset(
        scale_df.T,
        dataset.coordinate_df,
    )
    d._normalizer = "cpm"
    return d


def logcpm_normalizer(dataset: STDataset) -> STDataset:
    """
    Perform logcpm on dataset.

    Parameters
    ==========
    dataset : STDataset
        STDataset to be normalized.

    Returns
    =======
    dataset : STDataset
        A STDataset instance with normalized expression matrix.

    """
    count_df = dataset.count_df.sparse.to_dense()
    scale_df = nlog(count_df.T * 10000 / count_df.T.sum() + 1)
    d = STDataset(
        scale_df.T,
        dataset.coordinate_df,
    )
    d._normalizer = "logcpm"
    return d
