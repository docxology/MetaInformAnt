# Spatial Clustering

Group spots into spatially coherent domains using expression + geometry.

## Available algorithms

| Algorithm | `method` | Key idea |
|-----------|----------|----------|
| Spatially‑constrained Leiden | `'leiden'` | Leiden on a graph that blends expression similarity with spatial proximity |
| Spatially‑constrained Louvain | `'louvain'` | Louvain variant |
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
(`spatial_graph`).

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
path and the matplotlib `Figure` is returned. To draw the domains on the tissue
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
spatial distance).

## Pitfalls

- **Resolution too high** → oversegmentation; combine similar domains manually.
- **Platform scale differences** — Xenium subcellular spots require `n_neighbors`
  increased to 30–50 to connect neighbours across cell boundaries.
