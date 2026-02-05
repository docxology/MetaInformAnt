# Agent Directives: spatial

## Role
Documentation agent for the spatial transcriptomics module.

## Module Scope
- Multi-platform I/O (10x Visium, MERFISH, 10x Xenium)
- Spatial clustering and domain detection
- Cell type deconvolution (NNLS, NMF)
- Spatial autocorrelation statistics (Moran's I, Geary's C, Getis-Ord G)
- Neighborhood enrichment and niche detection
- Ligand-receptor signaling analysis
- scRNA-seq integration and anchor-based mapping
- Spatial visualization (tissue overlays, gene maps)

## Key Source Files
- `src/metainformant/spatial/io/` - Platform-specific data loaders
- `src/metainformant/spatial/analysis/` - Clustering, deconvolution, autocorrelation, neighborhood
- `src/metainformant/spatial/integration/` - scRNA-seq mapping
- `src/metainformant/spatial/visualization/` - Spatial plots
