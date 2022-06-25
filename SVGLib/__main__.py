from argparse import ArgumentParser
from multiprocessing import cpu_count
from pathlib import Path

from SVGLib import STDataset
from SVGLib import run


def main():
    parser = ArgumentParser(
        prog="SVGLib",
        description="Find spatial variable genes for Spatial Trasncriptomics data.",
    )
    parser.add_argument(
        "count",
        help="path to count matrix file (shape: (spot * gene))",
    )
    parser.add_argument(
        "coordinate",
        help="path to coordinate file (shape: (spot * 2))",
    )
    parser.add_argument(
        "--count_transpose",
        action="store_true",
        help="transpose count matrix if specified",
    )
    parser.add_argument(
        "--coordinate_transpose",
        action="store_true",
        help="transpose coordinate file if specified",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=6,
        help="number of neighbors to build knn network (default: %(default)s)",
    )
    parser.add_argument(
        "--n_genes",
        type=int,
        default=1000,
        help="number of genes to find clusters (default: %(default)s)",
    )
    parser.add_argument(
        "--n_gene_clusters",
        type=int,
        default=8,
        help="number of gene clusters to find (default: %(default)s)",
    )
    parser.add_argument(
        "--savedir",
        default=".",
        help="path to save results (default: %(default)s)",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=cpu_count(),
        help="number of threads to run SVGLib (default: %(default)s)",
    )
    args = parser.parse_args()

    d = STDataset(
        args.count,
        args.coordinate,
        count_transpose=args.count_transpose,
        coordinate_transpose=args.coordinate_transpose,
    )
    d = run(
        d,
        k=args.k,
        n_genes=args.n_genes,
        n_gene_clusters=args.n_gene_clusters,
        cores=args.cores,
    )

    d.hotspot_df.to_csv(Path.joinpath(args.savedir, "hotspot_df.csv"))
    d.AI.to_csv(Path.joinpath(args.savedir, "AI.csv"))
    d.Di.to_csv(Path.joinpath(args.savedir, "Di.csv"))
    d.gene_cluster.to_csv(Path.joinpath(args.savedir, "gene_cluster.csv"))


if __name__ == "__main__":
    main()
