#################
svgbit quickstart
#################
.. image:: https://badge.fury.io/py/svgbit.svg
    :target: https://badge.fury.io/py/svgbit

For further detail, please visit our `API reference`_ on readthedocs.org

Folder ``slides`` contains 10X Visium slide layout GPR files used in our
publication. Please refer to our paper's supplementary table 1 for sample
index.



Installation
============
Install svgbit with ``pip``::

    pip install svgbit



Command-line Interface
======================
svgbit has a command line version. Just tape::

    svgbit --help

after installation, and you may get a short help massage::

    usage: svgbit [-h] [--k K] [--n_svgs N_SVGS] [--n_svg_clusters N_SVG_CLUSTERS]
                  [--he_image HE_IMAGE] [--save_dir SAVE_DIR] [--cores CORES] read_path

    Find spatial variable genes for Spatial Trasncriptomics data.

    positional arguments:
    read_path             Read Spatial Transcriptomics data. Support format in 10X
                          Space Ranger(dir named ``outs``) and anndata hdf5


    optional arguments:
    -h, --help            show this help message and exit
    --normalization NORMALIZATION
                          apply which normalization on read data. If None (default),
                          neither normalization will apply. Supported values: None, cpm,
                          logcpm (default: None)

    --k K                 number of nearest neighbors for KNN network (default: 6)
    --n_svgs N_SVGS       number of SVGs to find clusters (default: 1000)
    --n_svg_clusters N_SVG_CLUSTERS
    number of SVG clusters to find (default: 8)
    --he_image HE_IMAGE   path to H&E image. Only used for visualization (default: None)
    --save_dir SAVE_DIR   path to save results (default: .)
    --cores CORES         number of threads to run svgbit (default: 8)

Follow the introduction and results will save to --savedir.

.. note::
    svgbit will use all available CPUs as default. While python's child process
    will inherit all resources from parent, this may consume much memory. Specify
    ``cores`` keyword argument to avoid memory error.

.. note::
   svgbit may consume ~35 Gib memory when running with a 2980 spots, 32285 genes
   matrix with 8 cores.



Python Interface
================
svgbit has a set of python API. You may run svgbit through command line or
python. We recommend the usage of python API for more feature and convient
control of your input data.


Load Data
---------
svgbit could load data from Space Ranger output directory::

    import svgbit as sb
    dataset = sb.load_10X("spaceranger_output/outs")

Or load data from csv files::

    import svgbit as sb
    dataset = sb.STDataset(
        count_df="Data/count_df.csv",
        coordinate_df="Data/coor_df.csv",
        count_df_kwargs={"index_col": 0, "header": 0},
        coordinate_df_kwargs={"index_col": 0, "header": 0},
    )

.. note::
    ``sb.STDataset`` also accept ``pd.DataFrame`` and ``np.array`` as
    ``count_df`` and ``coordinate_df``. If ``str`` or ``pathlib.Path`` are
    given, svgbit would try to load data with ``pandas``.

.. note::
   When init STDataset instance in this way, svgbit assume that the shape
   of count matrix and coordinate file is  (spot * gene) and (spot * 2).
   Specify ``count_transpose`` or ``coordinate_transpose`` as ``True``
   when necessary.


Data Preprocessing
------------------
Genes with low variance may filter with::

    dataset = sb.filters.low_variance_filter(dataset, var=0)

.. note::
   With var=0, genes with no expression may be filtered. We recommend
   apply a var=0 filter first for better computational performing.

Genes with extremely high expressions usually show no pattern and may
distrub performing. Filter with::

    dataset = sb.filters.quantile_filter(dataset, 0.99)

svgbit alse has count normalization functions::

    dataset = sb.normalizers.logcpm_normailzer(dataset)

Feel free for choosing gene filters and data normalizers. Other filters
and normalizers are also provided. Visit our `API reference`_ for further
detail.


Run svgbit
----------
To perform full pipeline of svgbit, running::

    sb.run(dataset)


Visualization
-------------
Draw SVG heatmap and spot type distribution map with::

    sb.plot.svg_heatmap(dataset, save_path="heatmap.jpg", he_image="he_image.jpg")
    sb.plot.spot_type_map(dataset, save_path="spot_type.jpg", he_image="he_image.jpg")

Parameter ``he_image`` is optional. If not specified, hotspot discription
map will show without morphological information.


Gene combinations
-----------------
Users may find gene combinations with::

    sb.find_combinations(dataset, center_spots=1)

to find gene combinations for SVG cluster 1. A ``pd.DataFrame`` with gene colocalization
and exclusive result will be returned.


Details about sb.run()
----------------------
When you perform ``sb.run()``, sevaral steps will be done as below.
For further detail of calculation, please refer to our publication.

Acquire weight
::::::::::::::

To calculate hotspot matrix, svgbit needs a weight network which discribes
association across spots. svgbit uses k-nearest neighbors with 6 neighbors
as a default. You may pass key word argument ``k`` to ``sb.run()`` to
change this behavior.

In this step, ``sb.run()`` will execute ``STDataset.acquire_weight()``
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

SVG clusters is estimated by::

    dataset.find_clusters()

and save to ``svg_cluster`` attribute.

For further discription of hotspot, AI, Di and SVG cluster, please refer to
our manuscript.



Citation
========

If you use ``SVGbit`` in your academic publication, please cite:

> Hong, Y., Song, K., Zhang, Z. et al. The spatiotemporal dynamics of spatially variable genes in developing mouse brain revealed by a novel computational scheme. Cell Death Discov. 9, 264 (2023). https://doi.org/10.1038/s41420-023-01569-w

For test codes used in our paper, please visit
https://github.com/fuyu-sama/svgbit-test



.. _API reference: https://svgbit.readthedocs.io/en/latest/
