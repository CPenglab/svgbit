from __future__ import annotations

from numpy import log as nlog

from .STDataset import STDataset


def cpm_normalizer(dataset: STDataset) -> STDataset:
    scale_df = (dataset.count_df.T * 10000 / dataset.count_df.T.sum())
    return STDataset(
        scale_df.T,
        dataset.coordinate_df,
    )


def logcpm_normalizer(dataset: STDataset) -> STDataset:
    scale_df = nlog(dataset.count_df.T * 10000 / dataset.count_df.T.sum() + 1)
    return STDataset(
        scale_df.T,
        dataset.coordinate_df,
    )
