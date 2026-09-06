# Specification: singlecell

## Scope

Single-cell analysis module for METAINFORMANT. Provides preprocessing,
dimensionality reduction (PCA, UMAP, t-SNE, diffusion maps), clustering,
cell type annotation, doublet detection, and RNA velocity analysis.

## Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## Data Structures

- **Sub-packages**: analysis, celltyping, data, differential, doublet, velocity, visualization
- **Core Class**: `SingleCellData` — expression matrix with `obs`, `var`, `obsm`, `varm`, `uns`
- **Key Concepts**: scRNA-seq analysis, clustering, celltyping, doublet detection, RNA velocity

## API Definition

### Exports — `analysis/`

- `compute_pca` — PCA with auto-clamping of n_components
- `compute_umap` — UMAP embedding
- `compute_tsne` — t-SNE with automatic perplexity adjustment
- `compute_neighbors` — k-NN graph (euclidean, cosine, manhattan)
- `compute_diffusion_map` — diffusion map for trajectory analysis
- `select_hvgs` — highly variable gene selection (seurat, variance, cell_ranger)


### Exports — `data/`

- `SingleCellData` — core data container (`obs`, `var`, `uns`, `obsm`, `varm`, `obsp`, `varp`, `layers`)
- `load_count_matrix` — load h5ad, csv/tsv, or Matrix Market count matrices
- `calculate_qc_metrics` — per-cell and per-gene QC metrics (counts, genes, mito/ribo %)
- `filter_cells` / `filter_genes` — QC filtering (preserves `obsm`/`varm`/`layers`)
- `normalize_counts` — total-count normalization (also `median`, `size_factors`)
- `log_transform` — log1p transformation (configurable base)
- `scale_data` — z-score standardization with optional clipping
- `identify_highly_variable_genes` — HVG selection (seurat, cell_ranger flavors)
- `remove_batch_effects` — per-batch mean-centering or ComBat-style standardization

### Exports — `data/integration.py`

- `bbknn_integration` / `harmony_integration` / `combat_integration` — single-object batch correction on a `batch_key`
- `scanorama_integration` / `mnn_integration` — multi-dataset integration of `SingleCellData` lists
- `integrate_multiple_batches` — dispatcher for the methods above
- `evaluate_integration_quality` — batch mixing and distance-ratio metrics

### Exports — `doublet/`

- `detect_doublets` — simulation-based doublet scoring (synthetic doublets + KNN, Scrublet-style)
- `DoubletResult` — scores, predicted doublets, threshold, doublet rate

### Exports — `velocity/`

- `compute_velocity` — steady-state velocity (`u - gamma * s`) with R² quality filtering
- `velocity_embedding` — project velocity onto a low-dimensional embedding
- `velocity_pseudotime` — transition-probability diffusion into pseudotime
- `velocity_confidence` — per-gene and per-cell confidence metrics
- `fit_dynamical_model` — EM-like fit of alpha/beta/gamma kinetics per gene

### Exports — `io.py`

- `download_geo_supplementary`, `download_sra_fastqs`, `fetch_atlas_datasets` — dataset acquisition (wget/fasterq-dump)
- `run_salmon_alevin` — single-cell quantification via Salmon alevin

### Exports — `visualization/`

- `plot_umap`, `plot_tsne`, `plot_pca` — embedding plots with categorical/numeric coloring
- `plot_trajectory` — pseudotime trajectory overlay on an embedding
- `plot_marker_expression` — dotplot, heatmap, or violin for marker genes
- `plot_qc_metrics` — QC distributions and counts-vs-genes scatter
- `plot_cluster_comparison` — side-by-side comparison of clusterings
