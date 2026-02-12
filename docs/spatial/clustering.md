# Spatial Clustering

The spatial clustering module implements spatially-aware clustering for spatial transcriptomics data. It combines gene expression similarity with physical spatial proximity through graph construction and community detection algorithms.

## Key Concepts

### Spatial Graph Construction

Before clustering, a spatial neighborhood graph is built from spot/cell coordinates. Three methods are supported:

- **KNN (K-Nearest Neighbors)**: Connects each spot to its k nearest spatial neighbors. Default k=6 for Visium hexagonal grids.
- **Delaunay triangulation**: Connects spots based on Delaunay triangulation of coordinates. Produces a natural spatial tessellation.
- **Radius**: Connects all spots within a fixed distance threshold.

### Community Detection

The module supports two community detection algorithms on the spatial-expression graph:

- **Leiden**: Modern algorithm with guaranteed well-connected communities. Generally preferred for its better modularity optimization.
- **Louvain**: Classic algorithm that maximizes modularity through iterative community merging.

### Spatial Weighting

Clustering balances expression similarity and spatial proximity via a `spatial_weight` parameter (0.0 = expression only, 1.0 = spatial only). Typical values are 0.3-0.7 for spatial transcriptomics data.

### Spatial Domains

A BayesSpace-inspired approach identifies spatially coherent transcriptomic domains using PCA-reduced expression data combined with coordinate information.

## Data Structures

### SpatialClusterResult

```python
@dataclass
class SpatialClusterResult:
    labels: Any             # Cluster label per spot (integer array)
    n_clusters: int         # Number of clusters found
    method: str             # Clustering method used
    modularity: float       # Modularity score (graph-based methods)
    spatial_graph: Any      # Adjacency matrix used
    metadata: dict          # Additional result metadata
```

## Function Reference

### build_spatial_graph

```python
def build_spatial_graph(
    coordinates: Any,
    method: Literal["knn", "delaunay", "radius"] = "knn",
    n_neighbors: int = 6,
    radius: float | None = None,
) -> Any  # scipy sparse matrix
```

Build a spatial neighborhood graph from coordinates. Returns a scipy sparse adjacency matrix.

### spatial_cluster

```python
def spatial_cluster(
    expression: Any,
    coordinates: Any,
    method: str = "leiden",
    spatial_weight: float = 0.5,
    n_neighbors: int = 6,
    resolution: float = 1.0,
    n_components: int = 20,
) -> SpatialClusterResult
```

Full spatial clustering pipeline: builds spatial graph, constructs expression similarity graph, combines them with `spatial_weight`, and applies community detection.

### leiden_clustering

```python
def leiden_clustering(
    adjacency: Any,
    resolution: float = 1.0,
) -> SpatialClusterResult
```

Apply Leiden community detection to an adjacency matrix.

### louvain_clustering

```python
def louvain_clustering(
    adjacency: Any,
    resolution: float = 1.0,
) -> SpatialClusterResult
```

Apply Louvain community detection to an adjacency matrix.

### spatial_domains

```python
def spatial_domains(
    expression: Any,
    coordinates: Any,
    n_domains: int = 7,
    n_components: int = 15,
    spatial_weight: float = 0.5,
) -> SpatialClusterResult
```

BayesSpace-inspired spatial domain identification using PCA + KMeans with spatial coordinates.

## Usage Examples

```python
from metainformant.spatial import io, analysis

# Load Visium data
dataset = io.load_visium("path/to/spaceranger_output/")

# Full spatial clustering pipeline
result = analysis.spatial_cluster(
    dataset.expression,
    dataset.coordinates,
    method="leiden",
    spatial_weight=0.5,
    resolution=1.0,
)
print(f"Found {result.n_clusters} clusters (modularity: {result.modularity:.3f})")

# Build spatial graph independently
graph = analysis.build_spatial_graph(
    dataset.coordinates, method="knn", n_neighbors=6
)

# Apply Leiden directly on a custom graph
leiden_result = analysis.leiden_clustering(graph, resolution=0.8)

# Spatial domain identification
domains = analysis.spatial_domains(
    dataset.expression, dataset.coordinates, n_domains=7
)
```

## Configuration

- **Environment prefix**: `SPATIAL_`
- **Optional dependencies**: numpy, scipy (KDTree, Delaunay), scikit-learn (KMeans, PCA, NearestNeighbors)
- Resolution parameter controls cluster granularity (higher = more clusters)
- n_components controls PCA dimensionality before clustering

## Related Modules

- `spatial.io` -- Data loading that produces inputs for clustering
- `spatial.analysis.autocorrelation` -- Spatial statistics on cluster assignments
- `spatial.visualization` -- `plot_spatial_scatter` for visualizing cluster labels
- `spatial.deconvolution` -- Cell type composition within spatial clusters
