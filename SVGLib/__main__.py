from argparse import ArgumentParser
from multiprocessing import cpu_count
from pathlib import Path

from SVGLib import STDataset
from SVGLib import run


def main():
    parser = ArgumentParser(prog="SVGLib", description="")
    parser.add_argument(
        "--count",
        required=True,
        help="path to count matrix file (shape: (spot * gene))",
    )
    parser.add_argument(
        "--coordinate",
        required=True,
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
        "--neighbors",
        type=int,
        default=6,
        help="number of neighbors to build knn network (default: %(default)s)",
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
    run(d)
    d.acquire_weight(k=args.neighbors)
    d.acquire_hotspot(cores=args.cores)
    d.acquire_density(cores=args.cores)
    d.hotspot_df.to_csv(Path.joinpath(args.savedir, "hotspot_df.csv"))
    d.AI.to_csv(Path.joinpath(args.savedir, "AI.csv"))
    d.Di.to_csv(Path.joinpath(args.savedir, "Di.csv"))


if __name__ == "__main__":
    main()
