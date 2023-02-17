from __future__ import annotations

import gzip
from pathlib import Path

import anndata
import h5py
import numpy as np
import pandas as pd
import scipy.io
from .STDataset import STDataset


def load_10X(read_path, make_sparse=True) -> STDataset:
    """
    Load 10X Genomics Space Ranger outputs and generate STDataset.

    Parameters
    ==========
    read_path : str or pathlib.Path
        A location points to 10X outs dir. Assume directories
        ``filtered_feature_bc_matrix`` and ``spatial`` are in this
        path.

    make_sparse : bool, default True
        Whether to use sparse DataFrame in order to save memory.

    Returns
    =======
    dataset : STDataset
        A STDataset instance generated from read_path.
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

    count_df = pd.DataFrame.sparse.from_spmatrix(scipy.io.mmread(mtx_path), ).T
    count_df.index = spot_name
    count_df.columns = gene_name

    try:
        coor_df = pd.read_csv(position_path, index_col=0, header=None)
        spaceranger_version = "v1"
    except FileNotFoundError:
        position_path = Path.joinpath(read_path, "spatial",
                                      "tissue_positions.csv")
        coor_df = pd.read_csv(position_path, index_col=0, header=0)
        spaceranger_version = "v2"
    if spaceranger_version == "v1":
        array_coor = coor_df[[3, 2]]
        coor_df = coor_df[[5, 4]]
    elif spaceranger_version == "v2":
        array_coor = coor_df[["array_row", "array_col"]]
        coor_df = coor_df[["pxl_row_in_fullres", "pxl_col_in_fullres"]]
    coor_df.index.name = "barcode"
    coor_df.columns = ["X", "Y"]
    array_coor = array_coor.reindex(index=count_df.index)
    array_coor.columns = ["X", "Y"]
    coor_df = coor_df.reindex(index=count_df.index)

    if not make_sparse:
        count_df = count_df.sparse.to_dense()
    dataset = STDataset(count_df, coor_df, make_sparse=make_sparse)
    dataset._array_coordinate = array_coor
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
        A STDataset instance generated from read_path.
    """
    adata = anndata.read_h5ad(read_path, **kwargs)
    count_df = pd.DataFrame.sparse.from_spmatrix(
        adata.X,
        index=adata.obs.index,
        columns=adata.var.index,
    )
    coor_df = pd.DataFrame(adata.obsm["spatial"])
    coor_df.index = count_df.index
    coor_df.columns = ["X", "Y"]

    if isinstance(count_df.iloc[0, 0], np.float32):
        count_df = count_df.astype(np.float64)

    dataset = STDataset(count_df, coor_df)
    return dataset


def load_table(read_path, **kwargs) -> STDataset:
    """
    Load text file and generate STDataset.

    Support tables in following format:

    ======== === === ======== ============
     geneID   X   Y   counts   spot_name
                               (optional)
    ======== === === ======== ============
     gene_1   1   1      1      spot_1
     gene_2   1   1      2      spot_1
     gene_1   1   2      1      spot_2
     gene_3   1   2      3      spot_2
      ...    ... ...    ...       ...
    ======== === === ======== ============

    Parameters
    ==========
    read_path : str or pathlib.Path
        File name to read from.

    **kwargs
        Additional keyword arguments passed to pd.read_csv

    Returns
    =======
    dataset : STDataset
        A STDataset instance generated from read_dir.
    """
    read_df = pd.read_csv(read_path, **kwargs)
    count_dict = {}
    coor_dict = {}
    for it in read_df.iterrows():
        try:
            spot_name = it[1].iloc[3]
        except IndexError:
            spot_name = f"{it[1].iloc[0]}x{it[1].iloc[1]}"
        if spot_name not in coor_dict.keys():
            coor_dict[spot_name] = {"X": it[1].iloc[0], "Y": it[1].iloc[1]}
        if spot_name not in count_dict.keys():
            count_dict[spot_name] = {}
        if it[0] not in count_dict[spot_name].keys():
            count_dict[spot_name][it[0]] = it[1][2]
        else:
            count_dict[spot_name][it[0]] += it[1][2]
    count_df = pd.DataFrame.from_dict(count_dict).fillna(0).T
    coor_df = pd.DataFrame.from_dict(coor_dict).T
    coor_df = coor_df.reindex(index=count_df.index)

    return STDataset(count_df, coor_df)


def load_gef(read_path, slot="bin20") -> STDataset:
    """
    Load gef file and generate STDataset.

    Parameters
    ==========
    read_path : str or pathlib.Path
        File name to read from.

    Returns
    =======
    dataset : STDataset
        A STDataset instance generated from read_dir.
    """
    count_dict = {}
    coor_dict = {}
    f = h5py.File(read_path)
    for record in f["geneExp"][slot]["gene"]:
        gene_symbol = record[0].decode("utf8")
        begin_line = record[1]
        end_line = record[1] + record[2] - 1
        exp = f["geneExp"][slot]["expression"][begin_line:end_line]
        for line in exp:
            line = list(map(int, line))
            spot_name = f"{line[0]}x{line[1]}"
            if spot_name not in coor_dict.keys():
                coor_dict[spot_name] = {"X": line[0], "Y": line[1]}
            if spot_name not in count_dict.keys():
                count_dict[spot_name] = {}
            if gene_symbol not in count_dict[spot_name].keys():
                count_dict[spot_name][gene_symbol] = line[2]
            else:
                count_dict[spot_name][gene_symbol] += line[2]
    count_df = pd.DataFrame.from_dict(count_dict).fillna(0).T
    coor_df = pd.DataFrame.from_dict(coor_dict).T
    coor_df = coor_df.reindex(index=count_df.index)
    f.close()

    return STDataset(count_df, coor_df)
