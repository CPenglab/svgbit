#############
API reference
#############



STDataset
=========
.. currentmodule:: svgbit
.. autosummary::
    :toctree: _autosummary
    :maxdepth: 5
    :recursive:

    STDataset
    load_10X



Run svgbit within one function
==============================
.. currentmodule:: svgbit
.. autosummary::
    :toctree: _autosummary
    :maxdepth: 5
    :recursive:

    run



Filters
=======
.. currentmodule:: svgbit
.. autosummary::
    :toctree: _autosummary
    :maxdepth: 5
    :recursive:

    filters.low_variance_filter
    filters.high_expression_filter
    filters.quantile_filter



Normalizers
===========
.. currentmodule:: svgbit
.. autosummary::
    :toctree: _autosummary
    :maxdepth: 5
    :recursive:

    normalizers.logcpm_normalizer
    normalizers.cpm_normalizer



Visualization
=============
.. currentmodule:: svgbit
.. autosummary::
    :toctree: _autosummary
    :maxdepth: 5
    :recursive:

    plot.svg_heatmap
    plot.spot_type_map



STDataset Methods
=================
.. currentmodule:: svgbit
.. autosummary::
    :toctree: _autosummary
    :maxdepth: 5
    :recursive:

    STDataset.acquire_weight
    STDataset.acquire_hotspot
    STDataset.acquire_density
    STDataset.find_clusters




STDataset Attributes
====================
.. currentmodule:: svgbit
.. autosummary::
    :toctree: _autosummary
    :maxdepth: 5
    :recursive:

    STDataset.AI
    STDataset.Di
    STDataset.coordinate_df
    STDataset.count_df
    STDataset.svg_cluster
    STDataset.spot_type
    STDataset.genes
    STDataset.hotspot_df
    STDataset.n_genes
    STDataset.n_spots
    STDataset.spots
    STDataset.weight
    STDataset.weight_type




Inner functions, not recommended use directly
=============================================
.. currentmodule:: svgbit
.. autosummary::
    :toctree: _autosummary
    :maxdepth: 5
    :recursive:

    core.cluster.cluster
    core.density.hotspot_AI
    core.moran.local_moran
    core.moran.global_moran



