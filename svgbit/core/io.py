from __future__ import annotations

import gzip
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scipy.io
from .STDataset import STDataset


def load_10X(read_path) -> STDataset:
    """
    Load 10X Genomics Space Ranger outputs and generate STDataset.

    Parameters
    ==========
    read_path : str or pathlib.Path
        A location points to 10X outs dir. Assume directories
        ``filtered_feature_bc_matrix`` and ``spatial`` are in this
        path.

    Returns
    =======
    dataset : STDataset
        A STDataset instance generated from read_dir.
    """
    read_path = Path(read_path)
    mat_dir = Path.joinpath(read_path, "filtered_feature_bc_matrix")
    mtx_path = Path.joinpath(mat_dir, "matrix.mtx.gz")
    features_path = Path.joinpath(mat_dir, "features.tsv.gz")
    barcodes_path = Path.joinpath(mat_dir, "barcodes.tsv.gz")
    position_path = Path.joinpath(read_path, "spatial",
                                  "tissue_positions_list.csv")

    gene_name = []
    with gzip.open(features_path, "rt") as f:
        for line in f:
            line = line.strip()
            gene_name.append(line.split("\t")[1])

    spot_name = []
    with gzip.open(barcodes_path, "rt") as f:
        for line in f:
            line = line.strip()
            spot_name.append(line)

    count_df = pd.DataFrame.sparse.from_spmatrix(
        scipy.io.mmread(mtx_path),
    ).T.sparse.to_dense()
    count_df.index = spot_name
    count_df.columns = gene_name
    coor_df = pd.read_csv(position_path, index_col=0, header=None)
    coor_df = coor_df[[5, 4]]
    coor_df.index.name = "barcode"
    coor_df.columns = ["X", "Y"]
    coor_df = coor_df.reindex(index=count_df.index)

    dataset = STDataset(count_df, coor_df)
    return dataset


def load_anndata_h5(read_path, **kwargs) -> STDataset:
    """
    Load anndata saved h5ad file and generate STDataset.

    .. note::
        ``load_anndata_h5`` will try to use anndata.X as expression matrix.

    Parameters
    ==========
    read_path : str or pathlib.Path
        File name to read from.

    **kwargs
        Additional keyword arguments passed to anndata.read_h5ad

    Returns
    =======
    dataset : STDataset
        A STDataset instance generated from read_dir.
    """
    adata = anndata.read_h5ad(read_path, **kwargs)
    count_df = pd.DataFrame.sparse.from_spmatrix(
        adata.X,
        index=adata.obs.index,
        columns=adata.var.index,
    ).sparse.to_dense()
    coor_df = pd.DataFrame(adata.obsm["spatial"])
    coor_df.index = count_df.index
    coor_df.columns = ["X", "Y"]

    if isinstance(count_df.iloc[0, 0], np.float32):
        count_df = count_df.astype(np.float64)

    dataset = STDataset(count_df, coor_df)
    return dataset
