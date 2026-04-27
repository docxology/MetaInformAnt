# Configuration: Spatial

All keys are under the `spatial` namespace. Precedence: function kwargs →
`config.set()` → `config.yaml` → environment variables → defaults.

## Core graph settings

| Key | Env | Type | Default | Notes |
|-----|-----|------|---------|-------|
| `spatial.graph.method` | `SPATIAL_GRAPH_METHOD` | `str` | `knn` | `knn` or `radius` |
| `spatial.graph.n_neighbors` | `SPATIAL_N_NEIGHBORS` | `int` | `15` | KNN neighbours; range [3, 50] |
| `spatial.graph.radius` | `SPATIAL_RADIUS` | `float` | `100.0` | µm if method=radius |
| `spatial.graph.edge_weight` | `SPATIAL_EDGE_WEIGHT` | `str` | `gaussian` | `gaussian` or `binary` |
| `spatial.graph.sigma` | `SPATIAL_SIGMA` | `float` | `auto` | Gaussian σ; auto = 0.5 × mean NN dist |
| `spatial.graph.include_self` | `SPATIAL_INCLUDE_SELF` | `bool` | `False` | self‑loops in adjacency |
| `spatial.graph.use_geometry` | `SPATIAL_USE_GEOMETRY` | `bool` | `True` | include (x,y) in KNN metric |

## Clustering

| Key | Type | Default |
|-----|------|---------|
| `spatial.clustering.method` | `str` | `leiden` |
| `spatial.clustering.resolution` | `float` | `0.8` |
| `spatial.clustering.min_cluster_size` | `int` | `50` |
| `spatial.clustering.random_state` | `int` | `0` |

## Autocorrelation

| Key | Type | Default |
|-----|------|---------|
| `spatial.autocorrelation.n_permutations` | `int` | `999` |
| `spatial.autocorrelation.seed` | `int` | `42` |
| `spatial.autocorrelation.edge_correction` | `str` | `weights` |

## Deconvolution

| Key | Type | Default |
|-----|------|---------|
| `spatial.deconvolution.backend` | `str` | `stereoscope` |
| `spatial.deconvolution.device` | `str` | `cpu` |
| `spatial.deconvolution.max_epochs` | `int` | `400` |
| `spatial.deconvolution.lr` | `float` | `1e-3` |
| `spatial.deconvolution.batch_size` | `int` | `256` |

## Visualisation

| Key | Type | Default |
|-----|------|---------|
| `spatial.visualization.point_size` | `float` | `1.2` |
| `spatial.visualization.cmap` | `str` | `viridis` |
| `spatial.visualization.dpi` | `int` | `300` |

All keys validated at import; use `config.validate_schema('spatial')` to check.

### Environment variable table

| Variable | Expands to |
|----------|------------|
| `SPATIAL_GRAPH_METHOD` | `spatial.graph.method` |
| `SPATIAL_N_NEIGHBORS` | `spatial.graph.n_neighbors` |
| `SPATIAL_RADIUS` | `spatial.graph.radius` |
| `SPATIAL_DECONV_BACKEND` | `spatial.deconvolution.backend` |
| `SPATIAL_DECONV_DEVICE` | `spatial.deconvolution.device` |

Example YAML:

```yaml
spatial:
  graph:
    method: knn
    n_neighbors: 20
  clustering:
    resolution: 1.2
    min_cluster_size: 100
  deconvolution:
    backend: tangram
    max_epochs: 600
```
