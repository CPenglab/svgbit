from __future__ import annotations
from argparse import ArgumentParser
from multiprocessing import cpu_count
from pathlib import Path

from svgbit import load_10X, load_anndata_h5, run, plot
from svgbit.filters import low_variance_filter, quantile_filter
from svgbit.normalizers import logcpm_normalizer


def main() -> None:
    parser = ArgumentParser(
        prog="svgbit",
        description='''Find spatial variable genes for Spatial Trasncriptomics
                       data.''',
    )
    parser.add_argument(
        "read_path",
        help='''Read Spatial Transcriptomics data. Support format in 10X Space
                Ranger(dir named ``outs``) and anndata hdf5''',
    )
    parser.add_argument(
        "--k",
        type=int,
        default=6,
        help='''number of nearest neighbors for KNN network
                (default: %(default)s)''',
    )
    parser.add_argument(
        "--n_svgs",
        type=int,
        default=1000,
        help='''number of SVGs to find clusters (default: %(default)s)''',
    )
    parser.add_argument(
        "--n_svg_clusters",
        type=int,
        default=8,
        help='''number of SVG clusters to find (default: %(default)s)''',
    )
    parser.add_argument(
        "--he_image",
        default=None,
        help='''path to H&E image. Only used for visualization
                (default: %(default)s)''',
    )
    parser.add_argument(
        "--savedir",
        default=".",
        help='''path to save results (default: %(default)s)''',
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=cpu_count(),
        help='''number of threads to run svgbit (default: %(default)s)''',
    )
    args = parser.parse_args()

    read_path = Path(args.read_path)
    if read_path.is_dir():
        load_func = load_10X
    else:
        load_func = load_anndata_h5

    d = load_func(read_path)
    d = low_variance_filter(d)
    d = quantile_filter(d, 0.99)
    d = logcpm_normalizer(d)
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

    plot.svg_heatmap(d, Path.joinpath(savedir, "heatmap.jpg"), args.he_image)
    plot.spot_type_map(d, Path.joinpath(savedir, "type_map.jpg"),
                       args.he_image)


if __name__ == "__main__":
    main()
