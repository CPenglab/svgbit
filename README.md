# svgbit

## Installation

Install with `pip`:

```sh
pip install svgbit
```

## Usage

### Command-line interface
svgbit has a command line version. Just tape

```sh
svgbit --help
```

after installation, and you may get a short help massage

```
usage: svgbit [-h] [--count_transpose] [--coordinate_transpose] [--k K] [--n_svgs N_SVGS]
              [--n_svg_clusters N_SVG_CLUSTERS] [--he_image HE_IMAGE] [--savedir SAVEDIR]
              [--cores CORES]
              count coordinate

Find spatial variable genes for Spatial Trasncriptomics data.

positional arguments:
  count                 path to count matrix file (shape: (spot * gene))
  coordinate            path to coordinate file (shape: (spot * 2))

optional arguments:
  -h, --help            show this help message and exit
  --count_transpose     transpose count matrix if specified
  --coordinate_transpose
                        transpose coordinate file if specified
  --k K                 number of nearest neighbors for KNN network (default: 6)
  --n_svgs N_SVGS       number of SVGs to find clusters (default: 1000)
  --n_svg_clusters N_SVG_CLUSTERS
                        number of SVG clusters to find (default: 8)
  --he_image HE_IMAGE   path to H&E image. Only used for visualization (default: None)
  --savedir SAVEDIR     path to save results (default: .)
  --cores CORES         number of threads to run svgbit (default: 8)
```

Follow the introduction and results will save to --savedir. 

### Python API
svgbit has a set of python API. You may run svgbit through command line or python.
We recommend the usage of python API for more feature and full control of your input
data.

#### Load data

```python
import svgbit
dataset = svgbit.STDataset(
    count_df="Data/count_df.csv",
    coordinate_df="Data/coor_df.csv",
    count_df_kwargs={"index_col": 0, "header": 0},
    coordinate_df_kwargs={"index_col": 0, "header": 0},
)
```

`svgbit.STDataset` also accept `pd.DataFrame`and `np.array` as `count_df` and 
`coordinate_df`. If `str` or `pathlib.Path` are given, svgbit would try to load
data with `pandas`.

**Notice**: svgbit assume that the shape of count matrix and coordinate file is 
(spot * gene) and (spot * 2). Specify `count_transpose` and `coordinate_transpose`
as `True` when necessary. 

#### Run svgbit with one function

```python
svgbit.run(dataset)
```

Visit our API reference for further detail.

## API References

## Citation
