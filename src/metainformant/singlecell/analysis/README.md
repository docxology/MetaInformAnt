# Single-Cell Analysis

Clustering, dimensionality reduction, and trajectory inference for single-cell transcriptomics data.

## Contents

| File | Purpose |
|------|---------|
| `clustering.py` | Leiden, Louvain, k-means, hierarchical clustering, marker genes, cluster evaluation |
| `pca_methods.py` | PCA, ICA, factor analysis, highly variable gene selection |
| `nonlinear_methods.py` | t-SNE, UMAP, diffusion maps, MDS, neighbor graphs, embedding quality metrics |
| `dimensionality.py` | Backward-compatibility shim re-exporting `pca_methods` and `nonlinear_methods` |
| `trajectory.py` | Pseudotime inference: diffusion pseudotime, PAGA, Slingshot, branch detection |

## Key Functions

| Function | Description |
|----------|-------------|
| `leiden_clustering()` | Graph-based Leiden clustering on a kNN graph |
| `louvain_clustering()` | Graph-based Louvain community detection |
| `kmeans_clustering()` | K-means clustering on the expression matrix |
| `hierarchical_clustering()` | Agglomerative clustering via scipy linkage |
| `find_marker_genes()` | Differential expression to identify cluster markers |
| `compute_cluster_silhouette()` | Silhouette and cluster-quality metrics |
| `evaluate_clustering_performance()` | Cluster evaluation against optional ground truth |
| `select_hvgs()` | Highly variable gene selection (`seurat`, `cell_ranger`, `variance`) |
| `pca_reduction()` / `compute_pca()` | PCA dimensionality reduction |
| `ica_reduction()` | Independent component analysis |
| `factor_analysis_reduction()` | Factor analysis |
| `tsne_reduction()` / `umap_reduction()` | t-SNE and UMAP nonlinear embeddings |
| `mds_reduction()` | Multidimensional scaling on precomputed cell-cell distances |
| `diffusion_map_reduction()` / `compute_diffusion_map()` | Diffusion maps |
| `compute_neighbors()` | Neighbor graph stored in `uns["neighbors"]` |
| `compute_dimensionality_metrics()` | Embedding quality metrics (distance correlation, trustworthiness) |
| `dpt_trajectory()` | Diffusion pseudotime trajectory ordering |
| `paga_trajectory()` | PAGA graph abstraction for trajectory |
| `slingshot_trajectory()` | Slingshot lineage-based trajectory inference |
| `compute_pseudotime_from_dimensionality_reduction()` | Pseudotime along an existing embedding |
| `find_trajectory_branches()` | Branch detection in pseudotime-sorted cells |
| `compute_trajectory_entropy()` | Expression entropy along a trajectory |

## Usage

```python
from metainformant.singlecell.analysis.clustering import leiden_clustering, find_marker_genes
from metainformant.singlecell.analysis.trajectory import dpt_trajectory

data = leiden_clustering(data, resolution=1.0)
markers = find_marker_genes(data, groupby="leiden_cluster")
data = dpt_trajectory(data, root_cell=0)
```
