#################
svgbit quickstart
#################



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
                  [--he_image HE_IMAGE] [--savedir SAVEDIR] [--cores CORES] read_dir

    Find spatial variable genes for Spatial Trasncriptomics data.

    positional arguments:
    read_dir              a location points to 10X outs dir. Assume directories
                          ``filtered_feature_bc_matrix`` and ``spatial`` are in this path.

    optional arguments:
    -h, --help            show this help message and exit
    --k K                 number of nearest neighbors for KNN network (default: 6)
    --n_svgs N_SVGS       number of SVGs to find clusters (default: 1000)
    --n_svg_clusters N_SVG_CLUSTERS
    number of SVG clusters to find (default: 8)
    --he_image HE_IMAGE   path to H&E image. Only used for visualization (default: None)
    --savedir SAVEDIR     path to save results (default: .)
    --cores CORES         number of threads to run svgbit (default: 8)

Follow the introduction and results will save to --savedir.

.. note::
    svgbit will use all available CPUs as default. While python's child process
    will inherit all resources from parent, this may consume much memory. Specify
    ``cores`` keyword argument to avoid ``Cannot allocate memory`` error.

.. note::
   svgbit may consume ~35 Gib memory when running with a 2980 spots, 32285 genes
   matrix with 8 cores.



Python Interface
================
svgbit has a set of python API. You may run svgbit through command line or
python. We recommend the usage of python API for more feature and convient
control of your input data.


Run svgbit with one function
----------------------------
svgbit could load data from Space Ranger output directory::

    import svgbit
    dataset = svgbit.load_10X("spaceranger_output/outs")

Or load data from csv files::

    import svgbit
    dataset = svgbit.STDataset(
        count_df="Data/count_df.csv",
        coordinate_df="Data/coor_df.csv",
        count_df_kwargs={"index_col": 0, "header": 0},
        coordinate_df_kwargs={"index_col": 0, "header": 0},
    )

.. note::
    ``svgbit.STDataset`` also accept ``pd.DataFrame`` and ``np.array`` as
    ``count_df`` and ``coordinate_df``. If ``str`` or ``pathlib.Path`` are
    given, svgbit would try to load data with ``pandas``.

.. note::
   When init STDataset instance in this way, svgbit assume that the shape
   of count matrix and coordinate file is  (spot * gene) and (spot * 2).
   Specify ``count_transpose`` or ``coordinate_transpose`` as ``True``
   when necessary.

After data loading, run::

    svgbit.run(dataset)

to perform full pipeline of svgbit. Results will save as attributes of
``dataset``.

Visit our API references for further detail.


Visualization
-------------
Draw SVG heatmap with::

    svgbit.svg_heatmap(dataset, save_path="heatmap.jpg", he_image="he_image.jpg")

Parameter ``he_image`` is optional. If not specified, hotspot discription
map will show without morphological information.


Details about svgbit.run()
--------------------------
When you perform ``svgbit.run()``, sevaral steps will be done as below.
For further detail of calculation, please refer to our publication.

Acquire weight
::::::::::::::

To calculate hotspot matrix, svgbit needs a weight network which discribes
association across spots. svgbit uses k-nearest neighbors with 6 neighbors
as a default. You may pass key word argument ``k`` to ``svgbit.run()`` to
change this behavior.

In this step, ``svgbit.run()`` will execute ``STDataset.acquire_weight()``
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
