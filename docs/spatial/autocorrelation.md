# Spatial Autocorrelation

Quantify how gene expression clusters across tissue coordinates.

The statistics in `metainformant.spatial.analysis.autocorrelation` take plain arrays: the values of one variable (length `n`, for example the expression of one gene across spots) and an `n × n` spatial weights matrix (dense or any SciPy sparse format). They do not read AnnData objects. The significance tests are analytic (a normal approximation), so there is no permutation count, seed or edge-correction option.

## Spatial weights

```python
from metainformant.spatial.analysis.autocorrelation import spatial_weights_matrix

# Binary k-nearest-neighbour weights, symmetrised, rows standardised to sum to 1
weights = spatial_weights_matrix(coords, method='knn', k=6)

# Binary weights within a distance band, left unstandardised (for Getis-Ord G*)
raw_weights = spatial_weights_matrix(coords, method='binary', bandwidth=150.0, row_standardize=False)
```

`method` is `'knn'`, `'distance'` (inverse distance within `bandwidth`) or `'binary'` (within `bandwidth`); `bandwidth=None` uses 1.5 times the median nearest-neighbour distance, in the units of `coords`. The result is a SciPy CSR matrix. Row-standardised weights suit Moran's I and local Moran's I; Getis-Ord G* expects unstandardised ones.

## Moran's I (global)

Range roughly [-1, +1]: negative (dispersion), about zero (random), positive (clustering).

```python
from metainformant.spatial.analysis.autocorrelation import morans_i

values = expression[:, gene_names.index('EPCAM')]    # one value per spot
res = morans_i(values, weights)
print(res.I, res.expected_I, res.z_score, res.p_value, res.n)
```

`MoransIResult` holds `I`, `expected_I` (-1 / (n - 1)), `variance_I`, `z_score`, `p_value` (two-sided, from the normal approximation) and `n`. Constant values, or weights that sum to zero, give `I = 0` and `p_value = 1`. The normal approximation is less reliable for small `n` (a few hundred spots or fewer) and for extremely skewed values.

## Local Moran (LISA)

Spot-level contribution to global Moran: reveals hot and cold spots.

```python
from metainformant.spatial.analysis.autocorrelation import local_morans_i, spatial_weights_matrix

# coords: (n_spots, 2) array; mki67: 1-D expression of one gene per spot
weights = spatial_weights_matrix(coords, method='knn', k=6)
lisa = local_morans_i(mki67, weights, significance=0.05)
adata.obs['lisa_quad'] = lisa.cluster_labels  # 'HH', 'LL', 'HL', 'LH' or 'NS'
adata.obs['lisa_p'] = lisa.p_values           # per-spot p-values
```

`LocalMoransResult` holds `local_I`, `expected_I`, `z_scores`, `p_values`, `cluster_labels` (a list of strings) and `significance_level`. The p-values are per spot and are not adjusted for multiple testing.

Quadrants:
- HH = high expression surrounded by high  (hotspot core)
- LL = low expression surrounded by low  (coldspot core)
- HL = high surrounded by low  (potential outlier)
- LH = low surrounded by high  (edge of hotspot)
- NS = not significant at the chosen `significance`

Visualise:

```python
from metainformant.spatial.visualization.plots import plot_spatial_autocorrelation

# Visium coordinates are [row, col]; pass coordinate_order='xy' for Xenium and MERFISH coordinates
plot_spatial_autocorrelation(coords, lisa.local_I, 'lisa_map.png', cluster_labels=lisa.cluster_labels)
```

The figure is saved to the given path and the matplotlib `Figure` is returned. With `cluster_labels` it draws the categorical LISA map; without them it draws `local_I` on a colour scale.

## Geary's C

Focuses on dissimilarity between neighbours; range roughly 0 to 2:

```python
from metainformant.spatial.analysis.autocorrelation import gearys_c

gc = gearys_c(values, weights)
print(f"C = {gc.C:.3f}  (z = {gc.z_score:.2f}, p = {gc.p_value:.3g})")
```

`GearyCResult` holds `C`, `expected_C` (1), `variance_C`, `z_score`, `p_value` and `n`. `C < 1` indicates positive spatial autocorrelation, `C = 1` none and `C > 1` negative. The p-value is analytic.

## Getis-Ord Gi*

Hot-spot detection; the observation itself is included in its neighbourhood. Use unstandardised weights:

```python
from metainformant.spatial.analysis.autocorrelation import getis_ord_g, spatial_weights_matrix

raw_weights = spatial_weights_matrix(coords, method='binary', bandwidth=150.0, row_standardize=False)
gi = getis_ord_g(cxcl12, raw_weights, significance=0.05)
hotspots = gi.hot_spots   # boolean mask of significant hot spots (gi.cold_spots for cold spots)
print(f"Hot spots: {hotspots.sum()} / {len(hotspots)}")
```

`GetisOrdResult` holds `g_star` (on a z-score scale), `p_values`, `hot_spots`, `cold_spots` and `significance_level`.

## Variogram

```python
from metainformant.spatial.analysis.autocorrelation import spatial_variogram

vg = spatial_variogram(values, coords, n_bins=20, max_distance=None)
print(vg.bin_centers, vg.semivariance, vg.n_pairs)
print(vg.nugget, vg.sill, vg.range_param)
```

`max_distance=None` considers distances up to about half the extent of the point cloud. `VariogramResult` holds `bin_centers`, `semivariance`, `n_pairs` and the estimated `nugget`, `sill` and `range_param`.

## Gene panels: spatially variable genes

`detect_spatially_variable_genes` computes Moran's I for every gene on one shared weights matrix and ranks the genes:

```python
from metainformant.spatial.analysis.spatially_variable_genes import detect_spatially_variable_genes

svg = detect_spatially_variable_genes(
    expression,              # spots x genes (or genes x spots), dense or any SciPy sparse matrix
    coords,
    genes=gene_names,
    method='knn',
    k=6,
    min_mean_expression=0.1,
    min_variance=0.0,
)
print(svg.top_genes(10))                         # genes ranked by Moran's I
significant = [g for g, q in zip(svg.genes, svg.q_values) if q < 0.05]
```

`SVGResult` holds `genes` (sorted by Moran's I, descending), `morans_i`, `z_scores`, `p_values`, `q_values` (Benjamini-Hochberg adjusted p-values, `None` when unavailable) and `n_spots`. Genes below `min_mean_expression` or with a variance not above `min_variance` are skipped. A sparse matrix is densified one gene at a time.

## Performance notes

Every statistic is analytic, so the cost is one pass over the weights per gene. `detect_spatially_variable_genes` is a single-threaded loop with one Moran's I per gene: pass a gene subset for large panels (see [Performance](PERFORMANCE.md) for measured timings).
