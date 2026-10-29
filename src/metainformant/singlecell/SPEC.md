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

### Conventions

- Analysis functions return a modified copy and leave their input untouched. `select_hvgs`, `compute_neighbors` and `compute_diffusion_map` take a keyword-only `inplace=False`; with `inplace=True` they write into the argument and return it.
- `SingleCellData.copy()`, `filter_cells` and `filter_genes` copy `uns` recursively.
- Argument checks raise `ValidationError` (`metainformant.core.utils.errors`, not a `ValueError` or `TypeError` subclass); a few functions raise the built-in `ValueError` instead (for example `select_hvgs`); missing optional packages raise `ImportError`.
- Random draws made by the module itself use a local generator seeded from the function's `random_state`; numpy's global random state is never read or advanced.

### Exports — `analysis/`

- `compute_pca` — PCA with auto-clamping of n_components; `use_hvgs=True` restricts it to the genes flagged by `select_hvgs`
- `compute_umap` — UMAP embedding
- `compute_tsne` — t-SNE with automatic perplexity adjustment
- `compute_neighbors` — k-NN graph in `uns["neighbors"]` (any scikit-learn metric; `use_rep` is `"X"` or any `obsm` key, an unknown one raises, and the default `"X_pca"` falls back to `X` without a PCA, recorded in `params["use_rep"]`; `inplace` supported)
- `compute_diffusion_map` — diffusion map for trajectory analysis (`obsm["X_diffmap"]`; kNN edges that both cells list weigh once, the graph is symmetrized by the entrywise maximum; the solver is chosen by size, with a fixed ARPACK start vector; a stale neighbor graph raises `ValidationError`; `inplace` supported)
- `diffusion_map_reduction` — dense Gaussian-kernel diffusion map: `obs["DC1"...]` are the right eigenvectors of the Markov matrix, from the symmetric `D^-1/2 W D^-1/2`; raises `ValidationError` for non-finite values or identical cells instead of substituting random vectors
- `select_hvgs` — highly variable gene selection (`seurat`, `cell_ranger`, `seurat_v3`, `variance`; `seurat_v3` expects raw counts and ranks genes by a standardised variance stored in `var["variances_norm"]`); `inplace` supported; raises the built-in `ValueError` for an unknown flavor
- `leiden_clustering` — Leiden on a kNN graph; needs networkx, igraph and leidenalg and raises `ImportError` without them. Matrix input returns labels and can opt into a warned k-means fallback with `allow_kmeans_fallback=True`; `uns["leiden_clustering"]["method"]` records the method. `use_weights` passes the inverse-distance edge weights to leidenalg (unit weights when False); `louvain_clustering` honours it too
- `louvain_clustering`, `kmeans_clustering`, `hierarchical_clustering` (seeded subset above 1000 cells), `find_marker_genes` (exact Welch t-test or rank-sum; an undefined test has `p_value = NaN` and a `test_status`), `find_markers`, `compute_cluster_silhouette`, `evaluate_clustering_performance` — clustering and evaluation
- `compute_diffusion_pseudotime` / `dpt_trajectory` — diffusion pseudotime (Haghverdi et al. 2016) on the neighbor graph: scaled to [0, 1], `inf` for cells not connected to the root, root recorded in `uns["dpt"]`
- `paga_trajectory` — PAGA connectivity (Wolf et al. 2019) from the kNN graph: observed over expected inter-group edges, capped at 1 (the group pseudotime it adds is a hop-count heuristic); `slingshot_trajectory` — a simplified approximation (straight lines between cluster centroids)
- `compute_pseudotime_from_dimensionality_reduction`, `find_trajectory_branches`, `compute_trajectory_entropy` — embedding-distance pseudotime and pseudotime-axis helpers


### Exports — `data/`

- `SingleCellData` — core data container (`obs`, `var`, `uns`, `obsm`, `varm`, `obsp`, `varp`, `layers`)
- `load_count_matrix` — load h5ad, csv/tsv, or Matrix Market count matrices
- `calculate_qc_metrics` — per-cell and per-gene QC metrics (counts, genes, mito/ribo %)
- `filter_cells` / `filter_genes` — QC filtering (preserves `obsm`/`varm`/`layers`); `filter_cells` recomputes the per-gene QC columns in `var` on the retained cells and drops the cell-indexed `uns["neighbors"]` when cells are removed, recording why in `uns["neighbors_dropped"]`
- `normalize_counts` — total-count normalization (also `median`, `size_factors`); `uns["normalization"]["size_factors"]` holds one divisor per cell for every method, so `X_norm == X / size_factors * target_sum`
- `log_transform` — log1p transformation (configurable base)
- `scale_data` — z-score standardization with optional clipping
- `identify_highly_variable_genes` — HVG selection (seurat, cell_ranger flavors)
- `remove_batch_effects` — per-batch mean-centering or per-batch z-scoring (`"combat"`, which is not ComBat)

### Exports — `data/integration.py`

- `bbknn_integration` (a real batch-balanced kNN graph) / `harmony_integration` (per-batch centering in PCA space) / `combat_integration` (location-scale adjustment; `covariates` refused) — single-object batch correction on a `batch_key`
- `scanorama_integration` / `mnn_integration` — simplified multi-dataset integration of `SingleCellData` lists (the datasets must share genes in the same order)
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

- `download_geo_supplementary`, `download_sra_fastqs`, `fetch_atlas_datasets` — dataset acquisition (wget/fasterq-dump); `fetch_atlas_datasets` raises `NotImplementedError` for an ArrayExpress dataset (`Li-2022`) before creating or downloading anything
- `run_salmon_alevin` — single-cell quantification via Salmon alevin
- `run_seurat_integration`, `run_differential_expression`, `run_scenic_grnboost2`, `run_scenic_cistarget`, `run_scenic_aucell`, `filter_regulons_post_scenic`, `run_pseudotime_trajectory`, `run_go_enrichment`, `run_shinycell_export` — not implemented: each raises `NotImplementedError` with guidance and creates no file or directory at its output path

### Exports — `visualization/`

- `plot_umap`, `plot_tsne`, `plot_pca` — embedding plots with categorical/numeric coloring
- `plot_trajectory` — pseudotime trajectory overlay on an embedding
- `plot_marker_expression` — dotplot, heatmap, or violin for marker genes
- `plot_qc_metrics` — QC distributions and counts-vs-genes scatter
- `plot_cluster_comparison` — side-by-side comparison of clusterings
