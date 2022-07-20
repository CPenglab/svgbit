#################
svgene quickstart
#################



Installation
============
Install svgene with ``pip``::

    pip install svgene



Command-line Interface
======================
svgene has a command line version. Just tape::

    svgene --help

after installation, and you may get a short help massage::

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
      --he_image HE_IMAGE   path to H&E image. Only used for visualization (default: None)
      --savedir SAVEDIR     path to save results (default: .)
      --cores CORES         number of threads to run svgene (default: 8)

Follow the introduction and results will save to --savedir.



Python Interface
================
svgene has a set of python API. You may run svgene through command line or
python. We recommend the usage of python API for more feature and convient
control of your input data.


Run svgene with one function
----------------------------
Load data::
    
    import svgene
    dataset = svgene.STDataset(
        count_df="Data/count_df.csv",
        coordinate_df="Data/coor_df.csv",
        count_df_kwargs={"index_col": 0, "header": 0},
        coordinate_df_kwargs={"index_col": 0, "header": 0},
    )

``svgene.STDataset`` also accept ``pd.DataFrame`` and ``np.array`` as 
``count_df`` and ``coordinate_df``. If ``str`` or ``pathlib.Path`` are 
given, svgene would try to load data with ``pandas``.

**Notice**: svgene assume that the shape of count matrix and coordinate 
file is  (spot * gene) and (spot * 2). Specify ``count_transpose`` or
``coordinate_transpose`` as ``True`` when necessary. 

After data loading, run::

    svgene.run(dataset)

to perform full pipeline of svgene. Results will save as attributes of ``dataset``.

Visit our :any:`API references <../references/index>` for further detail.


Details about svgene.run()
--------------------------
When you perform ``svgene.run()``, sevaral steps will be done as below.
For further detail of calculation, please refer to our publication. 

Acquire weight
::::::::::::::

To calculate hotspot matrix, svgene needs a weight network which discribes
association across spots. svgene uses k-nearest neighbors with 6 neighbors
as a default. You may pass key word argument ``k`` to ``svgene.run()`` to
change this behavior.

In this step, ``svgene.run()`` will execute ``STDataset.acquire_weight()``
method with given parameters. You may also perform this step by::
    
    dataset.acquire_weight()

Weight will save as attribute ``weight`` of ``STDataset`` and detailed
discription of weight is saved to ``weight_type`` attribute. Users may 
provide a ``libpysal.weights.W`` instance as user-specified weight::

    dataset.weight = user_specified_weight

Acquire hotspot
:::::::::::::::

Hotspot matrix is estimated by::
    
    dataset.acquire_hotspot()

and save to ``hotspot_df`` attribute.

Density
:::::::

AI and Di value discribed in our paper will be calculate by::

    dataset.acquire_density()

and save to ``AI`` and ``Di`` attribute as ``pd.Series``. 

Find SVG clusters
:::::::::::::::::



Citation
========
