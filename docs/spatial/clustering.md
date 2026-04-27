# Spatial Clustering

Group spots into spatially coherent domains using expression + geometry.

## Available algorithms

| Algorithm | Function | Key idea |
|-----------|----------|----------|
| Spatially‑constrained Leiden | `spatial_leiden()` | Leiden with spatial edge weights |
| Spatially‑constrained Louvain | `spatial_louvain()` | Louvain variant |
| Spatially‑constrained k‑means | `spatial_kmeans()` | K‑means with coordinate penalty |
| Hierarchical (Ward) | `spatial_ward()` | Agglomerative on spatial graph |

Default: **spatial_leiden** with resolution 0.8; robust and widely adopted.

## Usage

```python
from metainformant.spatial.analysis.clustering import spatial_leiden
adata.obs['domain'] = spatial_leiden(
    adata,
    resolution=0.8,
    min_cluster_size=50,
    random_state=42,
)
```

`spatial_leiden()` works on an existing neighbour graph (`adata.obsp['connectivities']`),
so call `spatial.neighbors(adata)` first.

## Resolution tuning

The Leiden `resolution_parameter` controls cluster granularity:

```python
for res in [0.2, 0.5, 1.0, 1.5]:
    adata.obs[f'domain_r{res}'] = spatial_leiden(adata, resolution=res)
```

- Lower (0.2–0.5) → few large domains (tissue‑level compartments)
- Higher (1.0–2.0) → many small domains (fine‑grained niches)

Pick the smallest resolution where each domain has ≥ 50 spots and corresponds to
a biologically plausible area on the tissue image.

## Marker gene discovery

```python
from metainformant.spatial.analysis.clustering import find_markers
markers = find_markers(adata, groupby='domain', method='wilcoxon')
# markers is a dict {domain: DataFrame with genes, pvals, logfc}
top_per_domain = {d: m.head(10) for d,m in markers.items()}
```

Wilcoxon rank‑sum (Scanpy implementation) is fast; for >20 k spots consider
`method='t-test'`.

## Visualisation

```python
from metainformant.spatial.visualization import spatial_scatter
spatial_scatter(adata, color='domain', palette='tab20', show=True)

# Overlay boundaries
from metainformant.spatial.visualization import domain_outlines
domain_outlines(adata, linewidth=0.8, show=True)
```

## Merging similar domains

Post‑hoc merge domains with similar expression profiles:

```python
from metainformant.spatial.analysis.clustering import merge_domains
adata.obs['domain_merged'] = merge_domains(
    adata,
    threshold=0.8,          # cosine similarity ≥ 0.8
    min_size=100,           # keep small domains separate
)
```

## Spatial k‑means (alternative)

Fast but less accurate for irregular tissues:

```python
from metainformant.spatial.analysis.clustering import spatial_kmeans
adata.obs['kmeans'] = spatial_kmeans(adata, n_clusters=8, coord_weight=0.3)
```

`coord_weight` balances expression vs coordinates (0 = expression only, 1 = pure
spatial distance).

## Pitfalls

- **Too small `min_cluster_size`** → noisy singleton clusters. Increase to ≥ 30.
- **Resolution too high** → oversegmentation; combine similar domains manually.
- **Graph not built** → `spatial_leiden()` raises `ValueError('no graph')`; run
  `spatial.neighbors(adata)` first.
- **Platform scale differences** — Xenium subcellular spots require `n_neighbors`
  increased to 30–50 to connect neighbours across cell boundaries.
