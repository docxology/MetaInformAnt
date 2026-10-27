# Single-Cell Module Rules (`metainformant.singlecell`)

## Purpose
Single-cell transcriptomic analysis: quality control, normalization, dimensionality reduction, clustering, and trajectory analysis.

## Dependencies
- **Required**: `core` (expression matrices are passed in; `singlecell` does not import `rna` or `ml`)
- **Optional**: `scanpy`, `anndata`

## Source Structure
```
src/metainformant/singlecell/
 analysis/
 clustering.py, dimensionality.py, trajectory.py
 data/
 integration.py, preprocessing.py
 visualization/
 visualization.py
```

## Package Management
- **ALWAYS use `uv`** for all Python package management and environment operations
- Install optional dependencies: `uv add scanpy anndata`
- Use `uv run` to execute commands: `uv run pytest`, `uv run metainformant singlecell --help`
- **NEVER use `pip` directly**

## Key Submodules

### Preprocessing
- Quality control metrics
- Cell and gene filtering
- Normalization (log-normalization, scaling)

### Dimensionality
- Principal Component Analysis (PCA)
- t-SNE
- UMAP

**Patterns**:
```python
from metainformant.singlecell.analysis import dimensionality
from metainformant.ml.features import dimensionality as ml_dimensionality

# Single-cell PCA on a SingleCellData object: returns the data with obsm["X_pca"] and PC columns in obs
data = dimensionality.pca_reduction(data, n_components=50)

# Shared dimensionality reduction on a plain cells x genes matrix: (transformed data, fitted PCA model)
transformed, pca_model = ml_dimensionality.pca_reduction(expression_matrix, n_components=50)

# Nonlinear embeddings: dimensionality.umap_reduction(data), dimensionality.tsne_reduction(data)
```

### Clustering
- Leiden algorithm
- Louvain algorithm
- Community detection

### Trajectory (`trajectory`)
- Pseudotime analysis
- Trajectory inference
- Branch point detection

**Patterns**:
```python
from metainformant.singlecell.analysis import trajectory

# Needs an embedding (run PCA first); root_cell is the index of the starting cell
data = trajectory.compute_diffusion_pseudotime(data, root_cell=0)
# Adds obs["dpt_pseudotime"] and uns["dpt"]

branches = trajectory.find_trajectory_branches(data, pseudotime_col="dpt_pseudotime")
# Returns: dict with n_branches, branch_sizes, branch_labels, branches
```

### Integration (`integration`)
- Multi-sample batch correction
- Canonical correlation analysis (CCA)
- Mutual nearest neighbors (MNN) correction
- Harmony batch correction
- Data integration quality assessment

**Patterns**:
```python
from metainformant.singlecell.data import integration

# Integrate multiple samples (adds a "batch" column to obs);
# integration_method is "scanorama" (default), "mnn", "harmony" or "bbknn"
integrated = integration.integrate_multiple_batches(
    [sample1_data, sample2_data, sample3_data],
    integration_method="scanorama",
)

# Harmony batch correction
corrected = integration.harmony_integration(integrated, batch_key="batch")

# Alternative: ComBat batch correction
corrected = integration.combat_integration(integrated, batch_key="batch")
```

### Visualization (`visualization`)
- QC metrics visualization
- Dimensionality reduction plots (UMAP, t-SNE)
- Gene expression overlays
- Trajectory visualizations
- Cluster annotation plots

**Patterns**:
```python-snippet
from metainformant.singlecell.visualization.visualization import (
    plot_qc_metrics,
    plot_umap,
    plot_tsne,
    plot_pca,
    plot_marker_expression,
    plot_cluster_comparison,
    plot_trajectory,
)

# Create visualizations (each returns a matplotlib Figure)
fig = plot_qc_metrics(data)
fig = plot_umap(data, color="leiden_cluster")  # needs umap_reduction first; color is an obs column
fig = plot_marker_expression(data, marker_genes=["GENE1"], method="dotplot")
fig = plot_cluster_comparison(data, cluster_cols=["leiden_cluster", "batch"])
fig = plot_trajectory(data, trajectory_key="dpt")  # key in data.uns, e.g. after compute_diffusion_pseudotime
```

## Patterns

### I/O Operations
- **CRITICAL**: Always use `metainformant.core.io` for JSON/CSV/TSV operations
- Never use direct `import json` or `import csv`

**Pattern**:
```python
from metainformant.core import io

# Save analysis results
io.dump_json(results, "output/singlecell/analysis/clustering.json")
io.write_csv(df, "output/singlecell/analysis/markers.csv")
```

### Path Handling
- Always use `metainformant.core.io.paths` utilities for path operations
- Always resolve and validate paths using `paths.expand_and_resolve()` and `paths.is_within()`

**Pattern**:
```python
from metainformant.core.io import paths

# Expand and resolve paths
output_path = paths.expand_and_resolve("output/singlecell/analysis/")
```

### Type Hints
- Comprehensive type hints throughout
- Use `from __future__ import annotations` for forward references
- Union types: `str | Path | None`

### scanpy Integration
- Use `scanpy` when available
- Fallback to scikit-learn implementations
- Handle `AnnData` objects when available

**Pattern**:
```python-snippet
try:
    import scanpy as sc
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False

def normalize_data(data):
    if SCANPY_AVAILABLE:
        return sc.pp.normalize_total(data)
    else:
        # Fallback implementation
        return normalize_manual(data)
```

### AnnData Objects
- Support `AnnData` format when available
- Convert to/from standard formats (pandas, numpy)
- Preserve metadata in AnnData objects

### Output Paths
- Base: `output/singlecell/<analysis>/`
- Preprocessing: `output/singlecell/<analysis>/preprocessing/`
- Clustering: `output/singlecell/<analysis>/clustering/`
- Trajectory: `output/singlecell/<analysis>/trajectory/`
- Integration: `output/singlecell/<analysis>/integration/`
- Visualizations: `output/singlecell/<analysis>/plots/`

### Dimensionality Reduction
- `singlecell.analysis.dimensionality` provides PCA, t-SNE and UMAP for `SingleCellData`; `ml.features.dimensionality` offers the same methods for plain matrices
- Support both scanpy and scikit-learn APIs
- Return consistent data structures

## Configuration

- No module-specific config (uses core config)
- Output paths default to `output/singlecell/`
- Environment prefix: `SC_` (e.g., `SC_THREADS`, `SC_WORK_DIR`)

## Integration

- **Uses**: `core`
- **Used by**: `multiomics` (single-cell omics integration), `information` (single-cell information analysis)
- **External**: `scanpy`, `anndata` when available

## Reference

See Global Project Rules (root `.cursorrules`) for:
- Common directory structure and path handling
- Configuration patterns with env overrides
- Testing policy (real-implementation policy)
- Import patterns and code style
- Documentation guidelines

## Testing

- **STRICTLY REAL IMPLEMENTATION**: Test real implementations only (see Global Project Rules (root `.cursorrules`) real-implementation policy)
- Skip tests if scanpy/anndata unavailable using `@pytest.mark.external_tool`
- Test with real expression matrices from `data/`
- Write test outputs to `output/singlecell/test/` using `tmp_path` fixture
- Real implementations only - use real computations


