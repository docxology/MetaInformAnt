# Specification: singlecell

## 🎯 Scope

Single-cell analysis module for METAINFORMANT. Provides preprocessing,
dimensionality reduction (PCA, UMAP, t-SNE, diffusion maps), clustering,
cell type annotation, doublet detection, and RNA velocity analysis.

## 🧱 Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures

- **Sub-packages**: analysis, celltyping, data, differential, doublet, velocity, visualization
- **Core Class**: `SingleCellData` — expression matrix with `obs`, `var`, `obsm`, `varm`, `uns`
- **Key Concepts**: scRNA-seq analysis, clustering, celltyping, doublet detection, RNA velocity

## 🔌 API Definition

### Exports — `analysis/`

- `compute_pca` — PCA with auto-clamping of n_components
- `compute_umap` — UMAP embedding
- `compute_tsne` — t-SNE with automatic perplexity adjustment
- `compute_neighbors` — k-NN graph (euclidean, cosine, manhattan)
- `compute_diffusion_map` — diffusion map for trajectory analysis
- `select_hvgs` — highly variable gene selection (seurat, variance, cell_ranger)

### Exports — `data/`

- `SingleCellData` — core data container
- `normalize_counts` — total-count normalization
- `log_transform` — log1p transformation
- `scale_data` — z-score standardization
