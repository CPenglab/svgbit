import gzip
from pathlib import Path

import pandas as pd
import scipy
from .STDataset import STDataset


def load_10X(read_dir) -> STDataset:
    """
    Load 10X Genomics Space Ranger outputs and generate STDataset.

    Patameters
    ==========
    read_dir: str or pathlib.Path
        A location points to 10X outs dir. Assume the directories
        ``filtered_feature_bc_matrix`` and ``spatial`` are in this
        path.

    Returns
    =======
    dataset: STDataset
        A STDataset instance generated from read_dir.
    """
    read_dir = Path(read_dir)
    mat_dir = Path.joinpath(read_dir, "filtered_feature_bc_matrix")
    mtx_path = Path.joinpath(mat_dir, "matrix.mtx.gz")
    features_path = Path.joinpath(mat_dir, "features.tsv.gz")
    barcodes_path = Path.joinpath(mat_dir, "barcodes.tsv.gz")
    position_path = Path.joinpath(read_dir, "spatial",
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
        index=gene_name,
        columns=spot_name,
    ).T
    coor_df = pd.read_csv(position_path, index_col=0, header=None)
    coor_df = coor_df[[5, 4]]
    coor_df.index.name = "barcode"
    coor_df.columns = ["X", "Y"]
    coor_df = coor_df.reindex(index=count_df.index)

    dataset = STDataset(count_df, coor_df)
    return dataset
