# Single-Cell Analysis

Clustering, dimensionality reduction, and trajectory inference for single-cell transcriptomics data.

## Contents

| File | Purpose |
|------|---------|
| `clustering.py` | Leiden, Louvain, k-means, hierarchical clustering, marker genes, cluster evaluation |
| `pca_methods.py` | PCA, ICA, factor analysis, highly variable gene selection |
| `nonlinear_methods.py` | t-SNE, UMAP, diffusion maps, MDS, neighbor graphs, embedding quality metrics |
| `dimensionality.py` | Backward-compatibility shim re-exporting `pca_methods` and `nonlinear_methods` |
| `trajectory.py` | Pseudotime inference: diffusion pseudotime, embedding pseudotime, simplified PAGA and Slingshot approximations, pseudotime grouping |

## Conventions

- Functions that add results to a `SingleCellData` return a modified copy and leave their input untouched. `select_hvgs`, `compute_neighbors` and `compute_diffusion_map` also take a keyword-only `inplace=False`; with `inplace=True` they write into the argument and return that same object.
- Where a function draws random numbers itself (the hierarchical-clustering subset, the silhouette and embedding-metric subsamples), it takes a `random_state` and uses a local generator; numpy's global random state is never read or advanced. The default `None` is not reproducible.
- `compute_diffusion_map` and diffusion pseudotime start ARPACK from a fixed vector, so repeated calls on the same data agree.

## Key Functions

| Function | Description |
|----------|-------------|
| `leiden_clustering()` | Leiden clustering of a kNN graph. Needs networkx, igraph and leidenalg and raises `ImportError` without them; a plain matrix can opt into a warned k-means fallback with `allow_kmeans_fallback=True`. `uns["leiden_clustering"]["method"]` records the method; `use_weights` is recorded but edge weights are not passed to leidenalg |
| `louvain_clustering()` | Graph-based Louvain community detection (needs networkx and python-louvain) |
| `kmeans_clustering()` | K-means clustering on the expression matrix |
| `hierarchical_clustering()` | Agglomerative clustering via scipy linkage (`linkage_method=`); above 1000 cells it uses a random subset seeded by `random_state` |
| `find_marker_genes()` | Per-group marker statistics (`t-test` or `wilcoxon`) for a `SingleCellData` |
| `find_markers()` | Mean-difference marker ranking for any cells × genes matrix plus labels (no statistical test) |
| `compute_cluster_silhouette()` | Silhouette and cluster-quality metrics (`random_state` seeds the optional subsample) |
| `evaluate_clustering_performance()` | Cluster evaluation against optional ground truth |
| `select_hvgs()` | Highly variable gene selection (`seurat`, `cell_ranger`, `seurat_v3`, `variance`; `seurat_v3` expects raw counts); flags genes in `var["highly_variable"]` and stores a flavor-dependent standardised variance in `var["variances_norm"]` |
| `pca_reduction()` / `compute_pca()` | PCA dimensionality reduction (`use_hvgs=True` restricts it to the flagged genes) |
| `ica_reduction()` | Independent component analysis |
| `factor_analysis_reduction()` | Factor analysis |
| `tsne_reduction()` / `umap_reduction()` | t-SNE and UMAP nonlinear embeddings |
| `mds_reduction()` | Multidimensional scaling on precomputed cell-cell distances |
| `diffusion_map_reduction()` / `compute_diffusion_map()` | Diffusion maps; `diffusion_map_reduction` raises `numpy.linalg.LinAlgError` when the eigendecomposition fails |
| `compute_neighbors()` | Neighbor graph stored in `uns["neighbors"]` |
| `compute_dimensionality_metrics()` | Embedding quality metrics (distance correlation, trustworthiness; `random_state` seeds the subsample above 1000 cells) |
| `compute_diffusion_pseudotime()` / `dpt_trajectory()` | Diffusion pseudotime (Haghverdi et al. 2016) on the neighbor graph: scaled to [0, 1], `inf` for cells not connected to the root, root recorded in `uns["dpt"]` |
| `paga_trajectory()` | Simplified PAGA-style approximation: group connectivity counted from the row order of the cells |
| `slingshot_trajectory()` | Simplified Slingshot-style approximation: straight lines between cluster centroids |
| `compute_pseudotime_from_dimensionality_reduction()` | Euclidean distance from a root cell in an existing embedding |
| `find_trajectory_branches()` | Gap-based grouping of cells ordered by pseudotime |
| `compute_trajectory_entropy()` | Expression entropy along a trajectory |

## Usage

```python
from metainformant.singlecell.analysis.clustering import leiden_clustering, find_marker_genes
from metainformant.singlecell.analysis.trajectory import dpt_trajectory

data = leiden_clustering(data, resolution=1.0)
markers = find_marker_genes(data, groupby="leiden_cluster")
data = dpt_trajectory(data, root_cell=0)
```
