from __future__ import annotations
from argparse import ArgumentParser
from multiprocessing import cpu_count
from pathlib import Path

from PIL import Image

from svgbit import STDataset
from svgbit import run


def main():
    parser = ArgumentParser(
        prog="svgbit",
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
        help="number of nearest neighbors for KNN network (default: %(default)s)",
    )
    parser.add_argument(
        "--n_svgs",
        type=int,
        default=1000,
        help="number of SVGs to find clusters (default: %(default)s)",
    )
    parser.add_argument(
        "--n_svg_clusters",
        type=int,
        default=8,
        help="number of SVG clusters to find (default: %(default)s)",
    )
    parser.add_argument(
        "--he_image",
        default=None,
        help="path to H&E image. Only used for visualization (default: %(default)s)",
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
        help="number of threads to run svgbit (default: %(default)s)",
    )
    args = parser.parse_args()

    d = STDataset(
        args.count,
        args.coordinate,
        count_transpose=args.count_transpose,
        coordinate_transpose=args.coordinate_transpose,
        count_df_kwargs={"index_col": 0, "header": 0},
        coordinate_df_kwargs={"index_col": 0, "header": 0},
    )
    d = run(
        d,
        k=args.k,
        n_svgs=args.n_svgs,
        n_svg_clusters=args.n_svg_clusters,
        cores=args.cores,
    )

    savedir = Path(args.savedir)
    d.hotspot_df.to_csv(Path.joinpath(savedir, "hotspot_df.csv"))
    d.AI.to_csv(Path.joinpath(savedir, "AI.csv"))
    d.Di.to_csv(Path.joinpath(savedir, "Di.csv"))
    d.svg_cluster.to_csv(Path.joinpath(savedir, "svg_cluster.csv"))

    he_image = None if args.he_image is None else Image.open(args.he_image)
    d.svg_heatmap(Path.joinpath(savedir, "heatmap.jpg"), he_image)
    None if he_image is None else he_image.close()


if __name__ == "__main__":
    main()
