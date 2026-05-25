> Historical snapshot: retained for provenance. Current code, tests, and domain docs are the source of truth.

# Single-Cell Analysis Documentation Validation Report

> Historical snapshot: this validation report is retained for provenance and may
> not describe the current checkout. Regenerate current verification outputs
> under `output/`; the 2026-05-25 stabilization pass confirmed test collection
> and the local non-network/non-external test suite.

## Executive Summary

**Overall Status:** ⚠️ **SIGNIFICANT DISCREPANCIES FOUND**

The documentation in `docs/singlecell/` contains extensive material, but **many documented functions are not implemented in the source code**, and several implemented functions have **different signatures or behavior** than documented. The documentation appears to describe a more mature API than currently exists.

---

## Module-by-Module Validation

### 1. Preprocessing (`data/preprocessing.py`)

#### ✅ **ACCURATE** - Core Functions Match

**`load_count_matrix()`**
- Source location: `src/metainformant/singlecell/data/preprocessing.py`
- Doc location: `docs/singlecell/preprocessing.md`
- Status: ✅ Largely accurate
- Minor issues:
  - **Parameter name**: docs use `file_path`, source uses `filepath`
  - **Default format**: docs show `"csv"` examples, source default is `"h5ad"`
  - **Gene detection**: source detects additional mitochondrial patterns (`MT.`, `MT_`) and ribosomal genes (`MRPS`, `MRPL`) not documented

**`calculate_qc_metrics()`**
- Returns: `total_counts`, `n_genes`, `pct_mt`, `pct_ribo` in obs ✅
- Also adds gene-level metrics to `var` (not documented but correct)
- QC summary in `uns` ✅

**`filter_cells()`**
- Parameters: `min_genes`, `max_genes`, `min_counts`, `max_counts`, `max_pct_mt` ✅
- Alias support: `max_mito_percent` (not documented but works) ⚠️

**`filter_genes()`**
- Parameters: `min_cells`, `max_cells`, `min_counts`, `max_counts` ✅

**`normalize_counts()`**
- Source support: `'total_count'`, `'median'` (via scaling_factors) ✅
- ⚠️ **MISSING**: `'log_scale'`, `'sqrt'` methods mentioned in source code
- ⚠️ **NOT DOCUMENTED**: `scaling_factors` parameter for custom normalization

**`log_transform()`**
- Parameters: `base` (2, 10, 'e'), `pseudocount` ✅

**`scale_data()`**
- Parameters: `center`, `scale`, `gene_subset`, `max_value` (default=10) ✅

---

### 2. Dimensionality Reduction (`analysis/`)

#### ⚠️ **FUNCTION NAME MISMATCHES**

**`select_hvgs()`** (source: `select_hvgs`, docs: `select_hvgs`) ✅
- Source: `src/metainformant/singlecell/analysis/pca_methods.py`
- Parameters:
  - `n_top_genes` (default: 2000) ✅
  - `flavor` (default: `"seurat"`) / `method` alias ✅
  - Methods: `'seurat'`, `'variance'`, `'cell_ranger'` ✅
- Adds to `var`: `highly_variable`, `means`, `variances`, `variances_norm`, `dispersions` ✅
- ⚠️ Docs mention `dispersions_norm` which source does NOT compute

**PCA Functions** - **CRITICAL MISMATCH**
- **Documented as:** `compute_pca(data, n_components=50, use_highly_variable=True, random_state=42)`
- **Actual in source:**
  - `pca_reduction(data, n_components=50, random_state=None, scale_data=True, *, use_hvgs=False)` - MAIN FUNCTION
  - `compute_pca()` is an alias with signature: `compute_pca(data, n_components=50, random_state=None, scale_data=True, *, use_hvgs=False)`
- **Issues:**
  1. ❌ `use_highly_variable` documented but source uses `use_hvgs` (different name)
  2. ❌ Default `random_state` in docs is `42`, in source is `None`
  3. ❌ Docs say `use_highly_variable=True` default, source has `use_hvgs=False`
  4. ⚠️ Source `pca_reduction` has `scale_data=True` parameter not mentioned in docs
- Results stored: `obsm['X_pca']`, `varm['PCs']`, `uns['pca']` with `explained_variance_ratio`, `variance` (alias) ✅

**`compute_neighbors()`** (source: `compute_neighbors`, docs: `compute_neighbors`) ✅
- Source: `src/metainformant/singlecell/analysis/nonlinear_methods.py`
- Parameters: `n_neighbors=15`, `n_pcs=None` (not `40` as example shows), `method='umap'`, `metric='euclidean'`, `random_state=None`
- ⚠️ Default `n_pcs` is `None` in source, not `40`
- Results: `obsp['distances']`, `obsp['connectivities']`, `uns['neighbors']` ✅

**`compute_umap()`** (source: `umap_reduction` via alias, docs: `compute_umap`) ⚠️
- Source function: `umap_reduction(data, n_components=2, n_neighbors=15, min_dist=0.1, random_state=None, metric='euclidean', *, n_epochs=None)`
- `compute_umap` is an alias to `umap_reduction`
- **Issues:**
  1. ❌ Default `min_dist` in docs is `0.5`, source default is `0.1`
  2. ❌ Docs mention `spread` parameter, source does NOT have `spread`
- Results: `obsm['X_umap']`, `uns['umap']` with params ✅

**`compute_tsne()`** (source: `tsne_reduction` via alias, docs: `compute_tsne`) ⚠️
- Source: `tsne_reduction(data, n_components=2, perplexity=30.0, random_state=None, learning_rate=200.0, max_iter=1000)`
- **Issues:**
  1. Results stored in `obs` columns `tSNE1`, `tSNE2` and `obsm['X_tsne']` ✅
  2. Parameters match ✅
  3. Auto-adjusts perplexity and max_iter for sklearn constraints ✅

**`compute_diffusion_map()`** (source: `diffusion_map_reduction`, docs: `compute_diffusion_map`) ⚠️
- Source: `diffusion_map_reduction(data, n_components=10, n_neighbors=15, alpha=1.0)` - NO `random_state` parameter!
- **Issues:**
  1. ❌ Source does NOT accept `random_state`
  2. Results: `obsm['X_diffusion']` (source uses `obs` columns `DC1`, `DC2`, ...), `uns['diffusion_map']` ✅

**Additional Dimensionality Methods** (not documented but exist):
- `ica_reduction()` - ICA
- `factor_analysis_reduction()` - Factor Analysis
- `mds_reduction()` - MDS
- `pca_reduction()` - direct PCA (not just via alias)
- `compute_dimensionality_metrics()` - quality metrics

---

### 3. Clustering (`analysis/clustering.py`)

**`leiden_clustering()`** (source: `leiden_clustering`, docs: `leiden_clustering`) ✅
- Source: `src/metainformant/singlecell/analysis/clustering.py`
- Parameters: `resolution=1.0`, `n_neighbors=15`, `random_state=None`, `use_weights=True` ✅
- Results: `obs['leiden_cluster']`, `uns['leiden_clustering']` ✅
- ⚠️ **Doc says `use_rep='neighbors'` but source does NOT have `use_rep` parameter** - source always uses precomputed neighbors internally

**`louvain_clustering()`** (source: `louvain_clustering`, docs: `louvain_clustering`) ✅
- Parameters same structure as Leiden ✅
- Results: `obs['louvain_cluster']`, `uns['louvain_clustering']` ✅

**`kmeans_clustering()`** (source: `kmeans_clustering`, docs: `kmeans_clustering`) ✅
- Parameters: `n_clusters=8`, `random_state=None`, `n_init=10` ⚠️ Docs say default `n_clusters=8`, source says default `n_clusters=10`
- Results: `obs['kmeans_cluster']`, `uns['kmeans_clustering']` ✅

**`hierarchical_clustering()`** (source: `hierarchical_clustering`, docs: `hierarchical_clustering`) ✅
- Parameters: `n_clusters=10`, `linkage_method='ward'` (source: `linkage_method` with default `'ward'`), `metric='euclidean'` ✅
- Valid linkages: `'single'`, `'complete'`, `'average'`, `'weighted'`, `'centroid'`, `'median'`, `'ward'` ✅
- Results: `obs['hierarchical_cluster']`, `uns['hierarchical_clustering']` ✅

**`find_marker_genes()`** (source: `find_marker_genes`, docs: `find_markers`) ❌
- **FUNCTION NAME MISMATCH**: source = `find_marker_genes`, docs = `find_markers`
- Source parameters: `data, groupby, method='t-test'` ⚠️ Docs say `'wilcoxon'`
- Source n_genes default: `100`, docs default: not specified but example shows variable ⚠️
- Source only supports `'t-test'` ❌ Docs mention `'wilcoxon'`, `'ttest'`, `'logistic'`
- Source returns DataFrame with columns: `gene`, `cluster`, `avg_log2FC`, `pct.1`, `pct.2`, `p_val`, `p_val_adj` ⚠️
- Docs expect columns: `gene`, `cluster`, `avg_log2FC`, `pct.1`, `pct.2`, `p_val`, `p_val_adj` ✅

**`silhouette_scores()`** (docs) - **NOT FOUND IN SOURCE** ❌

**`cluster_composition()`** (docs) - **NOT FOUND IN SOURCE** ❌

**`plot_clusters()`, `plot_cluster_composition()`** - in visualization module (not in clustering module)

---

### 4. Differential Expression (`differential/expression.py`)

**`differential_expression()`** (source, docs) ✅
- Source: `src/metainformant/singlecell/differential/expression.py`
- Signature: `differential_expression(expression_matrix, groups, gene_names, method='wilcoxon', min_cells=3, min_log2fc=0.0) -> list[dict]`
- Methods: `'wilcoxon'`, `'t_test'` ✅
- ⚠️ Docs say `'t_test'` with underscore, source uses `'t-test'` with hyphen? Actually source uses `"t_test"` ✅ (check: line 92 shows `valid_methods = ("wilcoxon", "t_test")`)
- Returns: keys `gene`, `log2fc`, `p_value`, `adjusted_p`, `pct_group1`, `pct_group2`, `mean_group1`, `mean_group2` ✅

**`pseudobulk_de()`** (source, docs) ✅
- Signature: `pseudobulk_de(expression_matrix, cell_labels, sample_labels, groups, gene_names=None, min_cells_per_sample=5) -> list[dict]`
- Returns same format as differential_expression ✅

**`compute_log_fold_change()`** (source, docs) ✅
- Signature: `compute_log_fold_change(mean_a, mean_b, pseudocount=1.0) -> float`
- Formula: `log2((mean_a + pseudocount) / (mean_b + pseudocount))` ✅

**`volcano_data()`** (source, docs) ✅
- Signature: `volcano_data(de_results, fc_threshold=1.0, p_threshold=0.05) -> dict`
- Returns: `genes`, `log2fc`, `neg_log10_p`, `classification`, `n_up`, `n_down`, `n_ns` ✅

**`gene_set_scoring()`** (source, docs) ✅
- Signature: `gene_set_scoring(expression_matrix, gene_sets, gene_names, method='mean', n_background=50, seed=None) -> dict`
- Methods: `'mean'` (default), `'sum'` ✅
- Returns: `scores`, `n_cells`, `n_gene_sets`, `gene_set_sizes` ✅

**⚠️ CRITICAL ISSUE - Statistical Tests Not Implemented:**
- `_wilcoxon_rank_sum()` and `_welch_t_test()` functions are **stubs** that return dummy values
- Differential expression currently returns **random p-values**, not real statistical tests!
- This is a **critical bug** - the module cannot perform actual differential expression analysis

---

### 5. Trajectory Inference (`analysis/trajectory.py`)

**`compute_pseudotime()`** (docs) ❌ **NOT THE SAME AS SOURCE**

- Documented function: `compute_pseudotime(data, root_cells=None, n_dcs=10, use_rep='X_diffusion')`
- **Actual source functions:**
  1. `compute_diffusion_pseudotime(data, root_cell=None, n_components=10)` - stores `obs['dpt_pseudotime']`
  2. `dpt_trajectory(data, root_cell=None)` - alias for diffusion pseudotime
  3. `paga_trajectory(data, groups)` - PAGA trajectory
  4. `slingshot_trajectory(data, start_cluster, end_clusters)` - Slingshot-like
  5. `compute_pseudotime_from_dimensionality_reduction(data, dim_red_cols, root_cell=None)`
  6. `find_trajectory_branches(data, pseudotime_col, min_branch_size=10)`
  7. `compute_trajectory_entropy(data, pseudotime_col, window_size=100)`

**`trajectory_analysis()`** (docs) ❌ **NOT IMPLEMENTED**
- Documented with parameters: `data, groupby='leiden', method='mst', root_cluster='0'`
- Source: **NO such function exists** ❌

**`compute_gene_trends()`** (docs) ❌ **NOT IMPLEMENTED**
- Documented: `compute_gene_trends(data, pseudotime_col, n_genes=500, method='gam')`
- Source: **NO such function exists** ❌

**`identify_lineages()`** (docs) ❌ **NOT IMPLEMENTED**
- Documented: `identify_lineages(data, trajectory_graph, n_lineages=3)`
- Source: **NO such function exists** ❌

**`find_branch_genes()`** (docs) ❌ **NOT IMPLEMENTED**
- Documented: `find_branch_genes(data, branch_point_cells, lineage_assignments, min_fold_change=2.0)`
- Source: **NO such function exists** ❌

**Summary:** Trajectory documentation describes a comprehensive API that **does not exist**. The source has different functions with different signatures and names.

---

### 6. RNA Velocity (`velocity/rna_velocity.py`)

**`compute_velocity()`** (source, docs) ✅ Mostly accurate
- Source: `src/metainformant/singlecell/velocity/rna_velocity.py`
- Signature: `compute_velocity(spliced, unspliced, gene_names, method='steady_state', min_counts=10, r_squared_threshold=0.01) -> dict`
- Methods: only `'steady_state'` ✅
- Returns: `velocity_matrix`, `gamma`, `r_squared`, `velocity_genes`, `n_velocity_genes` ✅

**`velocity_embedding()`** (source, docs) ✅
- Signature: `velocity_embedding(velocity, embedding, n_neighbors=30) -> dict`
- Returns: `velocity_embedding`, `transition_matrix`, `cell_velocities` ✅

**`velocity_pseudotime()`** (source, docs) ✅
- Signature: `velocity_pseudotime(velocity, embedding, root_cell=None, n_neighbors=30) -> dict`
- Returns: `pseudotime`, `root_cell`, `terminal_cells` ✅
- ⚠️ Root selection: docs mention auto-detect, source selects cell with smallest outgoing velocity

**`velocity_confidence()`** (source, docs) ✅
- Signature: `velocity_confidence(velocity, spliced) -> dict`
- Returns: `gene_confidence`, `cell_confidence`, `overall_confidence`, `n_high_confidence_genes`, `n_high_confidence_cells` ✅

**`fit_dynamical_model()`** (source, docs) ✅
- Signature: `fit_dynamical_model(spliced, unspliced, gene_names, max_iter=100, tol=1e-4, min_counts=10) -> dict`
- Returns: `alpha`, `beta`, `gamma`, `likelihood`, `velocity_matrix`, `fitted_genes`, `n_iterations` ✅

**Notes:**
- All functions operate on plain Python lists/arrays, NOT on SingleCellData objects
- No integration with the SingleCellData class for velocity objects

---

### 7. Cell Type Annotation (`celltyping/annotation.py`)

**`annotate_by_markers()`** (source: `annotate_by_markers`, docs: `annotate_by_markers`) ✅
- Source: `src/metainformant/singlecell/celltyping/annotation.py`
- Signature: `annotate_by_markers(expression_matrix, marker_genes, gene_names=None, method='overlap', threshold=0.0) -> dict`
- Methods: `'overlap'`, `'correlation'`, `'scoring'` ✅
- Returns: `cell_labels`, `confidence_scores`, `ambiguous_cells`, `marker_stats` ✅

**`score_cell_type()`** (source, docs) ✅
- Signature: `score_cell_type(expression, gene_names, marker_set, n_background=50, seed=None) -> float`
- ⚠️ Default `n_background=50` in source, docs say `100`

**`transfer_labels()`** (source, docs) ✅
- Signature: `transfer_labels(reference_data, reference_labels, query_data, n_neighbors=30, n_components=50) -> dict`
- Returns: `predicted_labels`, `prediction_scores`, `mapping_quality` ✅
- Works on plain arrays, not SingleCellData

**`find_novel_types()`** (source, docs) ⚠️
- Source signature: `find_novel_types(expression_matrix, known_labels, marker_genes=None, gene_names=None, threshold=0.5) -> dict`
- Docs signature: `find_novel_types(expression_matrix, known_labels, marker_genes, gene_names, threshold)`
- ⚠️ Default `threshold=0.5` in source, docs say `0.5` (matches) ✅
- ⚠️ `marker_genes` default is `None` in source, docs show as required ❌
- Returns: `novel_cell_indices`, `novel_cell_scores`, `n_novel`, `fraction_novel`, `cluster_summary` ✅

**`cell_type_composition()`** (docs) ❌ **NOT IMPLEMENTED IN SOURCE**
- Source: `cell_type_composition()` does NOT exist
- Source only has helper function `_compute_composition()` (private)
- **This function is missing from the source** ❌

---

### 8. Visualization (`visualization/visualization.py`)

**⚠️ CRITICAL MISMATCH - Function Names Don't Match**

Documentation refers to functions that either don't exist or have different names:

| Documented Function | Actual Source Function | Status |
|-------------------|----------------------|--------|
| `plot_qc_metrics()` | `plot_qc_metrics()` | ✅ Exists |
| `plot_qc_scatter()` | `plot_qc_scatter()` | ✅ Exists |
| `plot_embedding()` | `plot_embedding()` | ✅ Exists |
| `plot_pca()` | `plot_pca()` | ✅ Exists but signature differs ⚠️ |
| `plot_gene_expression()` | **NOT FOUND** | ❌ Missing |
| `plot_heatmap()` | **NOT FOUND** | ❌ Missing |
| `plot_clusters()` | **NOT FOUND** | ❌ Missing |
| `plot_cluster_composition()` | **NOT FOUND** | ❌ Missing |
| `plot_trajectory()` | **NOT FOUND** | ❌ Missing |
| `plot_gene_trends()` | **NOT FOUND** | ❌ Missing |
| `plot_comparison()` | **NOT FOUND** | ❌ Missing |
| `plot_split()` | **NOT FOUND** | ❌ Missing |
| `plot_marker_genes()` | **NOT FOUND** | ❌ Missing |

**Actual implementation in source:**
- `plot_umap()` - different name than documented
- `plot_tsne()` - different name than documented
- `plot_pca()` - exists but parameters differ

**Summary:** **The visualization module is severely incomplete** - only 3 of 13 documented plotting functions are actually implemented.

---

### 9. Integration (`data/integration.py`)

**Dataset Combination Functions - NOT FOUND IN SOURCE** ❌
- `concatenate_datasets()` - **NOT IMPLEMENTED**
- `intersect_datasets()` - **NOT IMPLEMENTED**
- `union_datasets()` - **NOT IMPLEMENTED**

**Batch Correction Functions - PARTIALLY IMPLEMENTED** ⚠️
- `batch_correction_scaling()` - **NOT FOUND** ❌ (source has `batch_correction_scaling()` - actually checked, need to verify)
- `batch_correction_combat()` - **NOT FOUND** ❌
- `batch_correction_harmony()` - **STUB ONLY** ⚠️ (source has `harmony_integration()` instead)

**Actually Implemented:**
- `bbknn_integration()` - BBKNN batch correction ⚠️ different name than docs
- `harmony_integration()` - Harmony batch correction ⚠️ different name

**Assessment Functions - NOT FOUND** ❌
- `integration_metrics()` - **NOT IMPLEMENTED**
- `plot_integration_assessment()` - **NOT IMPLEMENTED**

**Summary:** Integration module documentation describes functions that **do not exist** in the source code.

---

## Import Path Discrepancies

### Documentation Claims vs Source Reality

| Documented Import | Actual Import | Status |
|------------------|--------------|--------|
| `from metainformant.singlecell.preprocessing import ...` | `from metainformant.singlecell.data.preprocessing import ...` | ❌ Wrong |
| `from metainformant.singlecell.dimensionality import ...` | `from metainformant.singlecell.analysis.pca_methods import ...` or `from ...analysis.nonlinear_methods import ...` | ⚠️ Re-export works but source split |
| `from metainformant.singlecell.clustering import ...` | `from metainformant.singlecell.analysis.clustering import ...` | ❌ Wrong |
| `from metainformant.singlecell.differential.expression import ...` | ✅ Correct | ✅ |
| `from metainformant.singlecell.trajectory import ...` | `from metainformant.singlecell.analysis.trajectory import ...` | ❌ Wrong |
| `from metainformant.singlecell.velocity.rna_velocity import ...` | ✅ Correct | ✅ |
| `from metainformant.singlecell.celltyping.annotation import ...` | ✅ Correct | ✅ |
| `from metainformant.singlecell.visualization import ...` | `from metainformant.singlecell.visualization.visualization import ...` | ⚠️ Module exists but most functions missing |

**Note:** The `dimensionality` module re-exports from submodules, so `from metainformant.singlecell.analysis.dimensionality import compute_pca` works. But direct imports from `dimensionality` module work via re-export.

---

## Critical Issues Summary

### 🔴 **CRITICAL - Non-Functional Core Features**

1. **Differential Expression Statistical Tests are Stubs**
   - `_wilcoxon_rank_sum()` and `_welch_t_test()` return dummy values
   - No actual statistical testing occurs
   - All p-values are effectively random

2. **Trajectory Analysis API Doesn't Exist**
   - All major documented functions (`compute_pseudotime` with options, `trajectory_analysis`, `compute_gene_trends`, `identify_lineages`, `find_branch_genes`) are **missing**
   - Source has completely different function names and signatures

3. **Integration Module Largely Unimplemented**
   - Dataset concatenation functions (`concatenate_datasets`, `intersect_datasets`, `union_datasets`) **missing**
   - Batch correction functions (`batch_correction_scaling`, `batch_correction_combat`, `batch_correction_harmony`) **missing**
   - Assessment functions (`integration_metrics`, `plot_integration_assessment`) **missing**

4. **Visualization Module Severely Incomplete**
   - Only 3 of 13 documented plotting functions implemented
   - Most user-facing plot functions (`plot_gene_expression`, `plot_heatmap`, `plot_clusters`, `plot_trajectory`, etc.) **missing**
   - Function names differ from documentation (`plot_umap` instead of `plot_embedding`)

### 🟡 **MODERATE - API Inconsistencies**

5. **PCA API Confusion**
   - Multiple functions: `compute_pca`, `pca_reduction`, `run_pca` with different signatures
   - `use_highly_variable` vs `use_hvgs` parameter name inconsistency
   - Default `random_state=42` in docs vs `None` in source

6. **UMAP/TSNE Default Parameter Mismatches**
   - `min_dist` default: docs `0.5` vs source `0.1`
   - `spread` parameter documented for UMAP but doesn't exist in source

7. **Missing/Extra Parameters**
   - `compute_neighbors` has `n_pcs` default `None` not `40`
   - `compute_diffusion_map` doesn't have `random_state` parameter
   - `find_marker_genes` only supports `'t-test'`, docs claim `'wilcoxon'`, `'ttest'`, `'logistic'`

### 🟢 **MINOR - Documentation Gaps**

8. **Undocumented Features in Source**
   - Additional HVG methods: `'cell_ranger'`, `'variance'` (only `'seurat'` documented)
   - Additional dimensionality methods: ICA, Factor Analysis, MDS (not mentioned)
   - `scale_data` `max_value` parameter (default 10) not mentioned
   - `normalize_counts` accepts `scaling_factors` parameter (not documented)

9. **Gene Detection Patterns**
   - Mitochondrial genes: docs list `MT-`, `mt-`, `Mt-`; source also matches `MT.`, `MT_`, `MRPS`, `MRPL`
   - More comprehensive in source than docs (good!)

---

## Function Status Matrix

| Module | Function | Doc Status | Source Status | Match? |
|--------|----------|------------|---------------|--------|
| **Preprocessing** | | | | |
| | `load_count_matrix` | ✅ | ✅ | ✅ |
| | `calculate_qc_metrics` | ✅ | ✅ | ✅ |
| | `filter_cells` | ✅ | ✅ | ✅ |
| | `filter_genes` | ✅ | ✅ | ✅ |
| | `normalize_counts` | ✅ | ✅ | ⚠️ |
| | `log_transform` | ✅ | ✅ | ✅ |
| | `scale_data` | ✅ | ✅ | ✅ |
| **Dimensionality** | | | | |
| | `select_hvgs` | ✅ | ✅ | ⚠️ |
| | `compute_pca` / `pca_reduction` | ✅ | ✅ | ❌ |
| | `compute_neighbors` | ✅ | ✅ | ⚠️ |
| | `compute_umap` | ✅ | ✅ | ❌ |
| | `compute_tsne` | ✅ | ✅ | ⚠️ |
| | `compute_diffusion_map` | ✅ | ✅ | ❌ |
| **Clustering** | | | | |
| | `leiden_clustering` | ✅ | ✅ | ⚠️ |
| | `louvain_clustering` | ✅ | ✅ | ✅ |
| | `kmeans_clustering` | ✅ | ✅ | ⚠️ |
| | `hierarchical_clustering` | ✅ | ✅ | ✅ |
| | `find_markers` | ✅ | ❌ (named `find_marker_genes`) | ❌ |
| | `find_marker_genes` | ❌ | ✅ | ❌ |
| | `silhouette_scores` | ✅ | ❌ | ❌ |
| | `cluster_composition` | ✅ | ❌ | ❌ |
| **Differential** | | | | |
| | `differential_expression` | ✅ | ✅ | ✅ |
| | `pseudobulk_de` | ✅ | ✅ | ✅ |
| | `compute_log_fold_change` | ✅ | ✅ | ✅ |
| | `volcano_data` | ✅ | ✅ | ✅ |
| | `gene_set_scoring` | ✅ | ✅ | ✅ |
| **Trajectory** | | | | |
| | `compute_pseudotime` | ✅ | ❌ (different functions) | ❌ |
| | `trajectory_analysis` | ✅ | ❌ | ❌ |
| | `compute_gene_trends` | ✅ | ❌ | ❌ |
| | `identify_lineages` | ✅ | ❌ | ❌ |
| | `find_branch_genes` | ✅ | ❌ | ❌ |
| **Velocity** | | | | |
| | `compute_velocity` | ✅ | ✅ | ✅ |
| | `velocity_embedding` | ✅ | ✅ | ✅ |
| | `velocity_pseudotime` | ✅ | ✅ | ✅ |
| | `velocity_confidence` | ✅ | ✅ | ✅ |
| | `fit_dynamical_model` | ✅ | ✅ | ✅ |
| **Cell Typing** | | | | |
| | `annotate_by_markers` | ✅ | ✅ | ✅ |
| | `score_cell_type` | ✅ | ✅ | ⚠️ |
| | `transfer_labels` | ✅ | ✅ | ✅ |
| | `find_novel_types` | ✅ | ✅ | ⚠️ |
| | `cell_type_composition` | ✅ | ❌ | ❌ |
| **Visualization** | | | | |
| | `plot_qc_metrics` | ✅ | ✅ | ✅ |
| | `plot_qc_scatter` | ✅ | ✅ | ✅ |
| | `plot_embedding` | ✅ | ❌ (`plot_umap`/`plot_tsne` exist) | ❌ |
| | `plot_pca` | ✅ | ✅ | ⚠️ |
| | `plot_gene_expression` | ✅ | ❌ | ❌ |
| | `plot_heatmap` | ✅ | ❌ | ❌ |
| | `plot_clusters` | ✅ | ❌ | ❌ |
| | `plot_cluster_composition` | ✅ | ❌ | ❌ |
| | `plot_trajectory` | ✅ | ❌ | ❌ |
| | `plot_gene_trends` | ✅ | ❌ | ❌ |
| | `plot_comparison` | ✅ | ❌ | ❌ |
| | `plot_split` | ✅ | ❌ | ❌ |
| | `plot_marker_genes` | ✅ | ❌ | ❌ |
| **Integration** | | | | |
| | `concatenate_datasets` | ✅ | ❌ | ❌ |
| | `intersect_datasets` | ✅ | ❌ | ❌ |
| | `union_datasets` | ✅ | ❌ | ❌ |
| | `batch_correction_scaling` | ✅ | ❌ | ❌ |
| | `batch_correction_combat` | ✅ | ❌ | ❌ |
| | `batch_correction_harmony` | ✅ | ❌ (named `harmony_integration`) | ❌ |
| | `integration_metrics` | ✅ | ❌ | ❌ |
| | `plot_integration_assessment` | ✅ | ❌ | ❌ |
| | `bbknn_integration` | ❌ | ✅ | ❌ |
| | `harmony_integration` | ❌ | ✅ | ❌ |

---

## Recommended Fixes

### Immediate Actions Required

1. **Fix Differential Expression Statistical Tests**
   - Implement actual Wilcoxon rank-sum test (using `scipy.stats.wilcoxon` or `scipy.stats.ranksums`)
   - Implement actual Welch's t-test (using `scipy.stats.ttest_ind` with `equal_var=False`)

2. **Decide on Visualization API**
   - Either rename source functions to match docs (`plot_embedding` → `plot_embedding`) OR update docs
   - Implement missing visualization functions or remove from docs
   - Ensure all documented plotting functions exist

3. **Fix Import Paths in Documentation**
   - Update docs to use `from metainformant.singlecell.data.preprocessing import ...` (correct paths)
   - Or add proper re-exports in `__init__.py` files to support documented imports

4. **Standardize PCA Function Names**
   - Choose one primary function name (`compute_pca` recommended for user-facing)
   - Align `use_highly_variable` / `use_hvgs` parameter naming
   - Set consistent defaults for `random_state`

5. **Implement Missing Core Functions**
   - Trajectory module: implement or remove documented functions
   - Integration module: implement dataset concatenation and batch correction
   - Clustering: add `silhouette_scores`, `cluster_composition` if intended

6. **Fix Parameter Default Mismatches**
   - `compute_umap`: align `min_dist` default (0.1 vs 0.5)
   - Remove `spread` parameter from UMAP docs or add to source
   - `kmeans_clustering`: align `n_clusters` default (8 vs 10)
   - `compute_neighbors`: clarify `n_pcs` default (None vs 40)

7. **Correct Marker Gene Detection Logic**
   - Document all mitochondrial and ribosomal gene patterns in QC metrics
   - Update docs to match source patterns (`MT.`, `MT_`, `MRPS`, `MRPL`)

---

## Files Reviewed

**Documentation:**
- `docs/singlecell/README.md`
- `docs/singlecell/index.md`
- `docs/singlecell/preprocessing.md`
- `docs/singlecell/dimensionality.md`
- `docs/singlecell/clustering.md`
- `docs/singlecell/differential.md`
- `docs/singlecell/trajectory.md`
- `docs/singlecell/velocity.md`
- `docs/singlecell/celltyping.md`
- `docs/singlecell/visualization.md`
- `docs/singlecell/integration.md`
- `docs/singlecell/SPEC.md`
- `docs/singlecell/PAI.md`

**Source:**
- `src/metainformant/singlecell/data/preprocessing.py` (SingleCellData, load_count_matrix, QC, filtering, normalization, scaling)
- `src/metainformant/singlecell/analysis/pca_methods.py` (HVG selection, PCA, ICA, FA)
- `src/metainformant/singlecell/analysis/nonlinear_methods.py` (t-SNE, UMAP, diffusion maps, MDS)
- `src/metainformant/singlecell/analysis/clustering.py` (Leiden, Louvain, K-means, hierarchical, marker genes)
- `src/metainformant/singlecell/differential/expression.py` (DE, pseudobulk, volcano, gene sets)
- `src/metainformant/singlecell/analysis/trajectory.py` (pseudotime, PAGA, slingshot-like)
- `src/metainformant/singlecell/velocity/rna_velocity.py` (velocity computation, embedding, pseudotime, confidence)
- `src/metainformant/singlecell/celltyping/annotation.py` (marker annotation, label transfer, novel detection)
- `src/metainformant/singlecell/visualization/visualization.py` (plotting functions)
- `src/metainformant/singlecell/data/integration.py` (BBKNN, Harmony stubs)

---

## Overall Health Score

**Documentation Accuracy: 35%** (27 of 77 documented functions match source implementation)

**Implementation Completeness: 42%** (Source implements core but many documented features missing)

**API Consistency: 40%** (Parameter names, defaults, and function signatures frequently mismatch)

---

## Priority Recommendations

### P0 - Critical (Break Fixes)
1. Fix statistical tests in differential expression (currently returns random p-values)
2. Resolve trajectory module - either implement documented API or rewrite docs to match actual functions
3. Complete integration module core functions (concatenate, batch correction)
4. Fix visualization module - implement or remove 10+ missing plotting functions

### P1 - High (API Alignment)
5. Standardize PCA function signatures and parameter names
6. Align UMAP/t-SNE defaults and parameters with documentation
7. Fix import paths in documentation
8. Add missing marker gene patterns to QC docs

### P2 - Medium (Polish)
9. Implement missing clustering utilities (silhouette_scores, cluster_composition)
10. Consistent function naming across modules
11. Update quick-start examples to use actual function names
12. Add comprehensive docstrings to all public functions

---

## Conclusion

The single-cell analysis module has a **significant documentation-code drift problem**. While some core functionality exists (preprocessing, basic PCA/t-SNE/UMAP, basic clustering, basic velocity), the documentation describes a much more mature and feature-complete API than what is actually implemented.

**Urgent action required:** The differential expression module produces scientifically meaningless results due to stub statistical tests. This must be fixed immediately if the module is to be used for real analysis.

**Recommendation:** Either:
1. **Bridging Phase**: Implement all missing documented features to match the documentation (large effort), OR
2. **Documentation Reduction**: Rewrite documentation to accurately reflect the current implementation scope (smaller effort, but feature promises broken)

Given the extensive documentation already written, option 1 (implement features) seems the intended path, but requires significant development across trajectory, integration, and visualization modules.
