# svgene

## Installation

Install with `pip`:

```sh
pip install svgene
```

## Usage

### Command-line interface
svgene has a command line version. Just tape

```sh
svgene --help
```

after installation, and you may get a short help massage

```
usage: svgene [-h] [--count_transpose] [--coordinate_transpose] [--k K]
              [--n_genes N_GENES] [--n_gene_clusters N_GENE_CLUSTERS] [--savedir SAVEDIR]
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
  --n_genes N_GENES     number of genes to find clusters (default: 1000)
  --n_gene_clusters N_GENE_CLUSTERS
                        number of gene clusters to find (default: 8)
  --savedir SAVEDIR     path to save results (default: .)
  --cores CORES         number of threads to run svgene (default: 8)
```

Follow the introduction and results will save to --savedir. 

### Python API

#### Load data

```python
import svgene
dataset = svgene.STDataset(
    count_df="Data/count_df.csv",
    coordinate_df="Data/coor_df.csv",
    count_df_kwargs={"index_col": 0, "header": 0},
    coordinate_df_kwargs={"index_col": 0, "header": 0},
)
```

`svgene.STDataset` also accept `pd.DataFrame`and `np.array` as `count_df` and 
`coordinate_df`. If `str` or `pathlib.Path` are given, svgene would try to load
data with `pandas`.

**Notice**: svgene assume that the shape of count matrix and coordinate file is 
(spot * gene) and (spot * 2). Specify `count_transpose` and `coordinate_transpose`
as `True` when necessary. 


#### Run svgene with one function

```python
svgene.run(dataset)
```

Visit our API reference for further detail.

#### Visualization

#TODO: still under processing...

## API References

## Citation
