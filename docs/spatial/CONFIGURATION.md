# Configuration: Spatial

The spatial module has no configuration namespace: no `spatial.*` config keys, no
`SPATIAL_*` environment variables and no YAML section. Nothing under
`src/metainformant/spatial` reads the environment or `metainformant.core.utils.config`
(which has no `set()` or `validate_schema()`), so there are no layers or precedence to
resolve. Every setting is a keyword argument of the function that uses it; the tables
below list the real parameters and their defaults.

## Core graph settings

| Parameter | Function | Type | Default | Notes |
|-----------|----------|------|---------|-------|
| `method` | `build_spatial_graph` | `str` | `'knn'` | `'knn'`, `'delaunay'` or `'radius'` |
| `n_neighbors` | `build_spatial_graph` | `int` | `6` | neighbours per spot for `'knn'` (6 matches the Visium hexagonal grid) |
| `radius` | `build_spatial_graph` | `float` | `None` | required for `method='radius'` (`ValueError` otherwise); in coordinate units |
| `method` | `spatial_weights_matrix` | `str` | `'knn'` | `'knn'` (binary), `'distance'` (inverse distance within `bandwidth`) or `'binary'` (within `bandwidth`) |
| `k` | `spatial_weights_matrix` | `int` | `6` | neighbours for `'knn'` |
| `bandwidth` | `spatial_weights_matrix` | `float` | `None` | cut-off for `'distance'` and `'binary'`; `None` = 1.5 × the median nearest-neighbour distance |
| `row_standardize` | `spatial_weights_matrix` | `bool` | `True` | rows sum to 1 (a spot with no neighbours keeps an all-zero row) |

Edge weights are fixed: a `'knn'` graph edge weighs `1/(1+d)`, where `d` is the distance
between its two endpoints, and the graph is symmetrised by the entrywise maximum. An edge
therefore has the same weight whether the two spots are mutual neighbours or only one lists
the other (a mutual pair at distance 10 weighs 1/11), and no cap is needed. Delaunay and
radius edges are binary, and graphs have no self-loops. Coordinates keep the units of the
loaded dataset (Visium: pixels), so a `radius` or `bandwidth` is in those units.

## Clustering

| Parameter | Function | Type | Default | Notes |
|-----------|----------|------|---------|-------|
| `method` | `spatial_cluster` | `str` | `'leiden'` | `'leiden'`, `'louvain'` or `'kmeans'` |
| `n_clusters` | `spatial_cluster` | `int` | `None` | required for `'kmeans'`; ignored by the graph methods |
| `resolution` | `spatial_cluster` | `float` | `1.0` | graph methods; higher = more clusters |
| `spatial_weight` | `spatial_cluster` | `float` | `0.5` | 0 = expression only, 1 = spatial only |
| `n_neighbors` | `spatial_cluster` | `int` | `6` | neighbours in the expression graph and, for `graph_method='knn'`, the spatial graph; unused by `'kmeans'` |
| `graph_method` | `spatial_cluster` | `str` | `'knn'` | `'knn'`, `'delaunay'` or `'radius'` (`'kmeans'` builds no graph and ignores it) |
| `radius` | `spatial_cluster` | `float` | `None` | keyword-only; the neighbourhood radius in coordinate units, required for `graph_method='radius'` with the graph methods (`ValueError` otherwise) and ignored for the other graph methods |
| `n_pcs` | `spatial_cluster` | `int` | `15` | PCA components of the expression |
| `seed` | `spatial_cluster` | `int` | `42` | random seed |
| `n_domains`, `n_pcs`, `n_neighbors`, `max_iterations`, `seed` | `spatial_domains` | | `None`, `15`, `6`, `100`, `42` | BayesSpace-inspired domain finder. `n_domains=None` chooses 2 to 15 domains with the gap statistic (10 reference datasets drawn with `seed`); `result.metadata["n_domains_selection"]` is `"gap_statistic"` or `"user"` |
| `resolution`, `n_iterations`, `seed` | `leiden_clustering` | | `1.0`, `-1`, `42` | takes an adjacency matrix, returns `(labels, modularity)`. `n_iterations=-1` (any negative value) runs until an iteration brings no improvement, `0` raises `ValueError`; the greedy fallback used without leidenalg ignores it (at most 50 sweeps) |
| `resolution`, `seed` | `louvain_clustering` | | `1.0`, `42` | same input and return value; falls back to the greedy routine without `python-louvain` |

Both `leiden_clustering` and `louvain_clustering` return the same `(labels, modularity)` tuple
whether the published package or the built-in greedy fallback ran; the choice is logged at INFO
level only, and `spatial_cluster` keeps the requested method name in `result.method`.

## Autocorrelation

| Parameter | Function | Type | Default | Notes |
|-----------|----------|------|---------|-------|
| `significance` | `local_morans_i`, `getis_ord_g` | `float` | `0.05` | level for the cluster / hot-spot labels |
| `n_bins`, `max_distance` | `spatial_variogram` | | `20`, `None` | |
| `method`, `k`, `min_mean_expression`, `min_variance` | `detect_spatially_variable_genes` | | `'knn'`, `6`, `None`, `0.0` | per-gene Moran's I; the expression array may be dense or any SciPy sparse matrix |
| `n_permutations` | `neighborhood_enrichment` | `int` | `1000` | the only configurable permutation count; `compute_ligand_receptor_interactions` always runs 100 permutations |
| `n_simulations` | `ripley_k` | `int` | `99` | complete-spatial-randomness simulations |
| `edge_correction` | `ripley_k` | `str` | `'isotropic'` | keyword-only; `'isotropic'` (Ripley's exact correction for the bounding-box window, applied at every number of points) or `'none'` (every pair weighs 1); any other value raises `ValueError`. `area` is the caller's and sets the intensity only |
| `seed` | `neighborhood_enrichment`, `ripley_k`, `niche_detection` | `int` | `42` | |
| `gene_names` | `compute_ligand_receptor_interactions` | `Sequence[str] \| None` | `None` | keyword-only; one name per expression column, required for a plain matrix and taken from a dataset or DataFrame otherwise (explicit names win). Ligand and receptor symbols are matched by exact, case-sensitive name; pairs with an absent gene are skipped and counted in `summary` |
| `seed` | `compute_ligand_receptor_interactions` | `int \| None` | `None` | keyword-only; `None` draws the permutations from a fresh generator seeded from operating-system entropy, so the p-values differ between calls; pass an integer for reproducible p-values |

`morans_i`, `gearys_c`, `local_morans_i` and `getis_ord_g` use analytic (normal
approximation) statistics: they take only the values, the weights and, for the local
ones, `significance`, so there is no permutation count, seed or edge-correction option.

### Seeds

The random draws of the spatial functions come from local generators: numpy's global random
state is never read or advanced. Most functions default their seed to `42`
(`spatial_cluster`, `spatial_domains`, `leiden_clustering`, `louvain_clustering`,
`neighborhood_enrichment`, `niche_detection`, `ripley_k`, `anchor_based_transfer`,
`impute_spatial_genes`), `identify_niches` names it `random_state` (default `42`), and two functions
default to `None`, which is not reproducible:

| Function | Seed parameter | Default | With the default |
|----------|----------------|---------|------------------|
| `communication.cell_communication.compute_ligand_receptor_interactions` | `seed` | `None` | permutation p-values differ between calls (the scores do not) |
| `deconvolution.spatial_deconvolution.niche_identification` | `seed` | `None` | the k-means initialisation, and therefore the niches, can differ between calls |

Pass an integer to either function for reproducible output.

## Deconvolution

| Parameter | Function | Type | Default | Notes |
|-----------|----------|------|---------|-------|
| `method` | `deconvolve_spots` | `str` | `'nnls'` | `'nnls'` or `'nmf'`; both solve one non-negative least-squares problem per spot |
| `alpha` | `deconvolve_spots` | `float` | `0.0` | L2 strength; used only by `'nmf'` |
| `cell_type_names` | `deconvolve_spots` | `list[str]` | `None` | defaults to `type_0`, `type_1`, … |
| `spatial_gene_names`, `reference_gene_names` | `deconvolve_spots` | `list[str]` | `None` | when both are given, spots and reference are restricted to the shared genes |
| `method` | `create_reference_profiles` | `str` | `'mean'` | `'mean'` or `'median'` over each cell type's cells |

`deconvolve_spots` in `spatial.deconvolution.spatial_deconvolution` is a separate
function that takes a reference dict and `method` `'nnls'` (default) or `'regression'`.
There are no Stereoscope or Tangram backends, no device selection and no training
parameters (epochs, learning rate, batch size).

| Parameter | Function | Type | Default | Notes |
|-----------|----------|------|---------|-------|
| `spatial_markers` | `validate_deconvolution` | `dict` | `None` | cell type → marker gene names; checking markers needs `expression` and gene names |
| `expression`, `gene_names` | `validate_deconvolution` | | `None`, `None` | keyword-only; the spots × genes matrix the proportions came from (or a dataset or DataFrame carrying its names) and one name per column; markers are matched by exact name |
| `correlation_method` | `validate_deconvolution` | `str` | `'pearson'` | keyword-only; `'pearson'` or `'spearman'` (needs scipy); any other value raises `ValueError` |

## Visualisation

| Parameter | Function | Type | Default | Notes |
|-----------|----------|------|---------|-------|
| `cmap` | `plot_spatial_scatter` | `str` | `'viridis'` | other plot functions have their own default (`'Reds'` in `plot_gene_expression_map`, `'hot'` in `plot_tissue_overlay`, `'RdBu_r'` in `plot_interaction_heatmap`) |
| `point_size` | `plot_spatial_scatter` | `float` | `10.0` | `plot_cell_type_map` and `plot_spatial_autocorrelation` default to `15.0` |
| `alpha` | `plot_spatial_scatter` | `float` | `0.8` | |
| `figsize` | `plot_spatial_scatter` | `tuple` | `(8, 8)` | |
| `coordinate_order` | every plot that draws coordinates | `str` | `'row_col'` | `'row_col'` (Visium: column 0 is the vertical axis) or `'xy'` (Xenium, MERFISH); anything else raises `ValueError`. `plot_gene_expression_map` and `plot_cell_type_map` default to `None`: the dataset's own `coordinate_order`, else `'row_col'` |

Every plot function takes `output_path`, saves the figure there at a fixed 150 dpi
(`dpi` is not a parameter) and returns the matplotlib Figure.

### Not implemented

None of the following exists in the module:

- Any `SPATIAL_*` environment variable (`SPATIAL_GRAPH_METHOD`, `SPATIAL_N_NEIGHBORS`,
  `SPATIAL_RADIUS`, `SPATIAL_EDGE_WEIGHT`, `SPATIAL_SIGMA`, `SPATIAL_INCLUDE_SELF`,
  `SPATIAL_USE_GEOMETRY`, `SPATIAL_DECONV_BACKEND`, `SPATIAL_DECONV_DEVICE`) or
  `spatial.*` config key, and `config.set()` / `config.validate_schema()`.
- Graph options `edge_weight`, `sigma`, `include_self` and `use_geometry`.
- Clustering options `min_cluster_size` and `random_state` (the clustering functions name
  the seed `seed`; only `niche.identification.identify_niches` takes `random_state`). Every
  seed defaults to 42 except the two `seed=None` defaults listed under [Seeds](#seeds).
- Autocorrelation options `n_permutations`, `seed` and `edge_correction` on Moran's I and
  its relatives.
- Deconvolution options `backend`, `device`, `max_epochs`, `lr` and `batch_size`.

`metainformant.core.utils.config.apply_env_overrides(cfg, prefix='SPATIAL')` would map
`SPATIAL_THREADS`, `SPATIAL_WORK_DIR` and `SPATIAL_LOG_DIR` onto `threads`, `work_dir`
and `log_dir` in a dict you pass it, but no spatial function reads those keys.

To keep settings in a file, load it yourself and pass the values as keyword arguments:

```yaml
# config/spatial.yaml — your own file
n_neighbors: 20
resolution: 1.2
spatial_weight: 0.5
```

```python-snippet
from metainformant.core.utils.config import load_mapping_from_file
from metainformant.spatial.analysis.clustering import spatial_cluster

settings = load_mapping_from_file('config/spatial.yaml')
result = spatial_cluster(expression, coordinates, **settings)
```
