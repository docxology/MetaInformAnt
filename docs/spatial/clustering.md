# Spatial Clustering

Group spots into spatially coherent domains using expression + geometry.

## Available algorithms

| Algorithm | `method` | Key idea |
|-----------|----------|----------|
| Spatially‑constrained Leiden | `'leiden'` | Leiden on a graph that blends expression similarity with spatial proximity |
| Spatially‑constrained Louvain | `'louvain'` | Louvain on the same blended graph |
| Spatially‑constrained k‑means | `'kmeans'` | K‑means on expression PCs combined with scaled coordinates |

All three are selected through `spatial_cluster()`.

Default: **leiden** (`resolution=1.0`); robust and widely adopted.

## Usage

```python-snippet
from metainformant.spatial.analysis.clustering import spatial_cluster

# expression: (n_spots, n_genes) matrix; coords: (n_spots, 2) spatial coordinates
result = spatial_cluster(
    expression,
    coords,
    method='leiden',
    resolution=0.8,
    seed=42,
)
adata.obs['domain'] = result.labels
```

`spatial_cluster()` builds its own graph: a k‑nearest‑neighbour graph in PCA space
(expression) is blended with a spatial neighbour graph, so no precomputed
neighbour graph is needed. `spatial_weight` (0–1, default 0.5) sets the weight of the
spatial graph. The returned `SpatialClusterResult` also holds `n_clusters`; for the
graph‑based methods it additionally carries `modularity` and the combined graph
(`spatial_graph`), and `metadata` records the settings that were used (`spatial_weight`,
`n_pcs`, `resolution`, `graph_method` and `radius` for the graph methods; `spatial_weight`,
`n_pcs` and `inertia` for k‑means).

## Spatial graph options

`graph_method` chooses how the spatial neighbour graph is built:

| `graph_method` | Graph | Edge weight |
|----------------|-------|-------------|
| `'knn'` (default) | the `n_neighbors` nearest spots of every spot | `1 / (1 + distance)` |
| `'delaunay'` | Delaunay triangulation of the coordinates | 1 |
| `'radius'` | all spots within `radius` of each other | 1 |

```python-snippet
from metainformant.spatial.analysis.clustering import spatial_cluster

result = spatial_cluster(expression, coords, graph_method='delaunay')
result = spatial_cluster(expression, coords, graph_method='radius', radius=150.0)   # same units as coords
```

- `radius` is required with `graph_method='radius'` (there is no data-derived default): without it the graph methods raise `ValueError`. It is ignored by the other graph methods and by `method='kmeans'`, which builds no graph.
- `n_neighbors` (default 6, which matches the Visium hexagonal grid) sets the size of the expression graph and, for `'knn'`, of the spatial graph; `method='kmeans'` does not use it.
- The `'knn'` graph is symmetrised with the entrywise maximum: an edge's weight depends only on the distance between its two spots, whether both list each other as neighbours or only one does. Two mutual neighbours at distance 10 therefore have weight 1/11, the same as a pair at distance 10 that only one of them lists.

## Which implementation runs

`method='leiden'` uses `leidenalg` (with `igraph`) when it can be imported, and `method='louvain'` uses `python-louvain` (the `community` module, with `networkx`). When the package is missing, both fall back to a built-in greedy modularity routine, a simplified version of the first phase of Louvain. The fallback is only logged at INFO level: the returned `SpatialClusterResult` does not record which implementation ran, and its `method` keeps the requested name (`'leiden'` or `'louvain'`). The same holds for the lower-level `leiden_clustering` and `louvain_clustering`, whose `(labels, modularity)` return value does not show it either. Install `leidenalg` and `python-louvain` (the `singlecell` extra installs both) when you need the published algorithms.

## Lower-level functions

```python-snippet
from metainformant.spatial.analysis.clustering import build_spatial_graph, leiden_clustering, louvain_clustering

graph = build_spatial_graph(coords, method='knn', n_neighbors=6)      # SciPy CSR adjacency
labels, modularity = leiden_clustering(graph, resolution=1.0, n_iterations=-1, seed=42)
labels, modularity = louvain_clustering(graph, resolution=1.0, seed=42)
```

For `leiden_clustering`, `n_iterations` is the number of leidenalg iterations: `-1` (any negative value) iterates until an iteration brings no improvement, and `0` raises `ValueError`. The greedy fallback ignores it and sweeps the nodes until none moves, at most 50 times.

`spatial_domains(expression, coords, n_domains=None, *, n_pcs=15, n_neighbors=6, max_iterations=100, seed=42)` is a simplified, BayesSpace-inspired alternative: it smooths the PCA coordinates over the spatial neighbours and re-clusters with k-means until the labels stop changing (no MCMC). It returns a `SpatialClusterResult` (`method='spatial_domains'`) whose `metadata` holds `iterations`, `n_neighbors`, `smoothing_factor` and `n_domains_selection`.

### Choosing the number of domains

With `n_domains=None` the number of domains is chosen with the gap statistic (Tibshirani, Walther and Hastie, 2001), on the same PCA-plus-scaled-coordinates features that seed the clustering:

- For k = 1, 2, ... the log within-cluster sum of squares of a k-means fit is compared with its expectation on 10 reference datasets drawn uniformly over a box aligned with the principal axes of the features. The reference datasets are drawn with `seed`, so the choice is reproducible.
- The answer is the smallest k whose gap is within one standard error of the next one's (the "firstSEmax" rule of R's `clusGap`), between 2 and 15 domains (and at most the number of spots minus one). The rule can find k = 1, which means no domain structure; the function still returns at least 2.
- `result.metadata['n_domains_selection']` is `'gap_statistic'` when the number was chosen this way and `'user'` when you passed `n_domains`.
- With only tens of spots per domain the statistic can prefer a few extra domains. `tests/spatial/test_spatial_domains_gap_statistic.py` recovers 2 to 6 planted, well-separated domains at 40 spots per domain, and notes that at 30 spots per domain it can keep preferring a few extra ones. Pass `n_domains` when the number is known.

## Resolution tuning

The Leiden `resolution` controls cluster granularity:

```python
for res in [0.2, 0.5, 1.0, 1.5]:
    adata.obs[f'domain_r{res}'] = spatial_cluster(expression, coords, resolution=res).labels
```

- Lower (0.2–0.5) → few large domains (tissue‑level compartments)
- Higher (1.0–2.0) → many small domains (fine‑grained niches)

Pick the smallest resolution where each domain has ≥ 50 spots and corresponds to
a biologically plausible area on the tissue image.

## Marker gene discovery

```python-snippet
from metainformant.singlecell.analysis.clustering import find_markers

# Rank genes per domain by in-domain vs out-of-domain mean expression
markers = find_markers(expression, result.labels, n_markers=10, gene_names=gene_names)
# markers is a dict {domain: [{gene, score, mean_expr, mean_expr_other, log_fold_change, ...}, ...]}
top_per_domain = {d: [m['gene'] for m in rows] for d, rows in markers.items()}
```

`find_markers()` is shared with the single‑cell module: it works on any
(observations × genes) matrix plus labels. It ranks genes by the difference in mean
expression and does not run a statistical test.

## Visualisation

```python-snippet
from metainformant.spatial.visualization.plots import plot_spatial_scatter
plot_spatial_scatter(coords, result.labels.astype(str), 'domains.png', title='Spatial domains')
```

String labels get one colour per domain (`tab20`); integer labels would be drawn on a
continuous colour scale, hence the `astype(str)`. The figure is saved to the given
path and the matplotlib `Figure` is returned. `plot_spatial_scatter` reads `coords` as
`[row, col]` by default (Visium); for Xenium or MERFISH coordinates pass
`coordinate_order='xy'`. To draw the domains on the tissue
image, use `plot_tissue_overlay()` (see [Visualization](visualization.md)).

## Spatial k‑means (alternative)

Fast but less accurate for irregular tissues:

```python-snippet
from metainformant.spatial.analysis.clustering import spatial_cluster
adata.obs['kmeans'] = spatial_cluster(
    expression, coords, n_clusters=8, method='kmeans', spatial_weight=0.3
).labels
```

`spatial_weight` balances expression vs coordinates (0 = expression only, 1 = pure
spatial distance). `n_clusters` is required for `method='kmeans'` (`ValueError` otherwise) and
not used by the graph methods.

## Pitfalls

- **Resolution too high** → oversegmentation; combine similar domains manually.
- **Platform scale differences** — the default `n_neighbors=6` matches the Visium hexagonal grid. For single-cell platforms such as Xenium, tune `n_neighbors` or use `graph_method='radius'` with a `radius` in the units of the coordinates.
