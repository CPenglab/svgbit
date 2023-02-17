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
        A new STDataset instance with normalized expression matrix.

    """
    try:
        count_df = dataset.count_df.sparse.to_dense()
    except AttributeError:
        pass
    try:
        scale_df = (count_df.T * 10000 / count_df.T.sum())
    except AttributeError:
        pass
    return_dset = STDataset(
        scale_df.T,
        dataset.coordinate_df,
        check_duplicate_genes=False,
        sort_spots=False,
    )
    return_dset._normalizer = "cpm"
    return_dset._array_coordinate = dataset._array_coordinate
    return return_dset


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
        A new STDataset instance with normalized expression matrix.

    """
    try:
        count_df = dataset.count_df.sparse.to_dense()
    except AttributeError:
        pass
    scale_df = nlog(count_df.T * 10000 / count_df.T.sum() + 1)
    return_dset = STDataset(
        scale_df.T,
        dataset.coordinate_df,
        check_duplicate_genes=False,
        sort_spots=False,
    )
    return_dset._normalizer = "logcpm"
    return_dset._array_coordinate = dataset._array_coordinate
    return return_dset
