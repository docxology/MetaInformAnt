# Single-Cell Analysis

Clustering, dimensionality reduction, and trajectory inference for single-cell transcriptomics data.

## Contents

| File | Purpose |
|------|---------|
| `clustering.py` | Leiden, Louvain, k-means, hierarchical clustering and marker genes |
| `dimensionality.py` | PCA, t-SNE, UMAP, diffusion maps, ICA, factor analysis |
| `trajectory.py` | Pseudotime inference: diffusion pseudotime, PAGA, Slingshot |

## Key Functions

| Function | Description |
|----------|-------------|
| `leiden_clustering()` | Graph-based Leiden clustering on kNN graph |
| `louvain_clustering()` | Graph-based Louvain community detection |
| `kmeans_clustering()` | K-means clustering on reduced dimensions |
| `find_marker_genes()` | Differential expression to identify cluster markers |
| `pca_reduction()` | PCA dimensionality reduction |
| `umap_reduction()` | UMAP nonlinear embedding |
| `tsne_reduction()` | t-SNE nonlinear embedding |
| `compute_diffusion_map()` | Diffusion map for trajectory analysis |
| `dpt_trajectory()` | Diffusion pseudotime trajectory ordering |
| `paga_trajectory()` | PAGA graph abstraction for trajectory |
| `slingshot_trajectory()` | Slingshot lineage-based trajectory inference |

## Usage

```python
from metainformant.singlecell.analysis.clustering import leiden_clustering, find_marker_genes
from metainformant.singlecell.analysis.trajectory import dpt_trajectory

data = leiden_clustering(data, resolution=1.0)
markers = find_marker_genes(data, groupby="cluster")
data = dpt_trajectory(data, root_cell=0)
```
