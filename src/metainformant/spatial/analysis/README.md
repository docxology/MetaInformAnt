# Analysis

Spatial statistics and analysis algorithms for spatial transcriptomics data, providing autocorrelation testing, spatially-aware clustering, cell type deconvolution, and neighborhood enrichment analysis.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports autocorrelation, clustering, deconvolution, neighborhood, spatially_variable_genes |
| `autocorrelation.py` | Moran's I, Geary's C, LISA, Getis-Ord G*, variograms |
| `clustering.py` | Spatial graph construction, Leiden/Louvain clustering, BayesSpace-style domains |
| `deconvolution.py` | NNLS and NMF-based cell type deconvolution for multi-cell spots |
| `neighborhood.py` | Neighborhood enrichment, co-localization, Ripley's K, niche detection |
| `spatially_variable_genes.py` | Per-gene Moran's I ranking (`detect_spatially_variable_genes`) |

## Key Functions

| Function | Description |
|----------|-------------|
| `autocorrelation.spatial_weights_matrix()` | Build spatial weights matrix (KNN, distance-band) |
| `autocorrelation.morans_i()` | Global Moran's I spatial autocorrelation test (values, weights; analytic p-value) |
| `autocorrelation.local_morans_i()` | Local Moran's I (LISA) for hotspot detection |
| `autocorrelation.getis_ord_g()` | Getis-Ord G* statistic for clustering |
| `clustering.build_spatial_graph()` | Build KNN/Delaunay/radius spatial graph; the KNN graph weighs an edge `1 / (1 + distance)` and is symmetrised with the entrywise maximum |
| `clustering.spatial_cluster()` | Spatially-aware clustering combining expression and location; `graph_method='radius'` requires `radius=` |
| `clustering.leiden_clustering()` / `louvain_clustering()` | Community detection on an adjacency matrix, returning `(labels, modularity)`; they fall back to a greedy modularity routine without `leidenalg` / `python-louvain` and the return value does not say so. `leiden_clustering(n_iterations=-1)` runs until convergence and `0` raises `ValueError` |
| `clustering.spatial_domains()` | BayesSpace-inspired spatial domain identification; without `n_domains` the number (2 to 15) is chosen with the gap statistic and recorded in `metadata["n_domains_selection"]` (`"gap_statistic"` or `"user"`) |
| `deconvolution.deconvolve_spots()` | Deconvolve spot expression into cell type proportions |
| `neighborhood.neighborhood_enrichment()` | Compute cell type co-localization enrichment |
| `neighborhood.ripley_k()` | Ripley's K function for spatial point patterns, with a CSR envelope; `edge_correction="isotropic"` (default) applies Ripley's exact correction over the bounding box of the points at every number of points, `"none"` weighs every pair 1 |
| `neighborhood.niche_detection()` | Identify spatial niches from cell type composition |

## Usage

```python
from metainformant.spatial.analysis import autocorrelation, clustering, neighborhood

W = autocorrelation.spatial_weights_matrix(coordinates, method="knn", k=6)
result = autocorrelation.morans_i(expression_values, W)
clusters = clustering.spatial_cluster(expression, coordinates, n_clusters=8, method="kmeans")
enrichment = neighborhood.neighborhood_enrichment(cell_types, coordinates)
```
