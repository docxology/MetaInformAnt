# Analysis

Spatial statistics and analysis algorithms for spatial transcriptomics data, providing autocorrelation testing, spatially-aware clustering, cell type deconvolution, and neighborhood enrichment analysis.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports autocorrelation, clustering, deconvolution, neighborhood |
| `autocorrelation.py` | Moran's I, Geary's C, LISA, Getis-Ord G*, variograms |
| `clustering.py` | Spatial graph construction, Leiden/Louvain clustering, BayesSpace-style domains |
| `deconvolution.py` | NNLS and NMF-based cell type deconvolution for multi-cell spots |
| `neighborhood.py` | Neighborhood enrichment, co-localization, Ripley's K, niche detection |

## Key Functions

| Function | Description |
|----------|-------------|
| `autocorrelation.spatial_weights_matrix()` | Build spatial weights matrix (KNN, distance-band) |
| `autocorrelation.morans_i()` | Global Moran's I spatial autocorrelation test |
| `autocorrelation.local_morans_i()` | Local Moran's I (LISA) for hotspot detection |
| `autocorrelation.getis_ord_g()` | Getis-Ord G* statistic for clustering |
| `clustering.build_spatial_graph()` | Build KNN/Delaunay/radius spatial graph |
| `clustering.spatial_cluster()` | Spatially-aware clustering combining expression and location |
| `clustering.spatial_domains()` | BayesSpace-style spatial domain identification |
| `deconvolution.deconvolve_spots()` | Deconvolve spot expression into cell type proportions |
| `neighborhood.neighborhood_enrichment()` | Compute cell type co-localization enrichment |
| `neighborhood.ripley_k()` | Ripley's K function for spatial point patterns |
| `neighborhood.niche_detection()` | Identify spatial niches from cell type composition |

## Usage

```python
from metainformant.spatial.analysis import autocorrelation, clustering, neighborhood

W = autocorrelation.spatial_weights_matrix(coordinates, method="knn", k=6)
result = autocorrelation.morans_i(expression_values, W)
clusters = clustering.spatial_cluster(expression, coordinates, n_clusters=8)
enrichment = neighborhood.neighborhood_enrichment(coordinates, cell_types)
```
