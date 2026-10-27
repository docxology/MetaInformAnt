# Single-Cell Dimensionality Reduction

The dimensionality-reduction functions live in two modules: `metainformant.singlecell.analysis.pca_methods` (highly variable genes, PCA, ICA, factor analysis) and `metainformant.singlecell.analysis.nonlinear_methods` (neighbor graph, UMAP, t-SNE, diffusion maps, MDS, embedding metrics). `analysis.dimensionality` re-exports both for older imports.

## Conventions

- **Copies by default.** The functions that add results to a `SingleCellData` return a modified copy and leave their input untouched. `select_hvgs`, `compute_neighbors` and `compute_diffusion_map` also take a keyword-only `inplace` flag: with `inplace=True` they write into the argument and return that same object, which saves a copy of the expression matrix on large datasets.
- **Coordinates are stored twice.** PCA, UMAP and t-SNE write the coordinates to `data.obsm` (`X_pca`, `X_umap`, `X_tsne`) and as `obs` columns (`PC1`, `UMAP1`, `tSNE1`, ...). The plotting functions read the `obs` columns.
- **Seeds.** `random_state` is handed to scikit-learn or umap-learn by PCA, ICA, factor analysis, t-SNE, UMAP and MDS, which then seed themselves (pass an integer for reproducible coordinates). The one random draw this module makes itself, the cell subsample of `compute_dimensionality_metrics`, uses a local generator and never touches numpy's global random state. `compute_diffusion_map` starts ARPACK from a fixed vector, so repeated calls on the same data give identical coordinates, signs included.
- **Dense work.** HVG selection, PCA, UMAP, t-SNE and the neighbor graph convert `data.X` to a dense array internally.

## Highly Variable Genes (HVGs)

### select_hvgs()

Flag the genes with the highest biological variability:

```python
from metainformant.singlecell.analysis.pca_methods import select_hvgs

# Seurat-style selection (the default flavor)
data = select_hvgs(data, n_top_genes=2000)

# Other flavors; `method` is an alias of `flavor`
data = select_hvgs(data, n_top_genes=2000, flavor="variance")
data = select_hvgs(data, n_top_genes=2000, flavor="cell_ranger")

# Variance-stabilising selection of Stuart et al. (2019); needs raw counts in data.X
data = select_hvgs(data, n_top_genes=2000, flavor="seurat_v3")

# Write the annotations into `data` instead of returning a copy
select_hvgs(data, n_top_genes=2000, inplace=True)
```

**Parameters:**

- `n_top_genes`: number of genes to flag (default 2000)
- `flavor` (alias `method`): `"seurat"` (default), `"cell_ranger"`, `"seurat_v3"` or `"variance"`; anything else raises the built-in `ValueError` (this function does not use the package's `ValidationError`)
- `inplace`: keyword-only, default `False`

**Flavors:**

- `"seurat"`: fits a straight line of log10 variance against log10 mean over the genes with non-zero mean and ranks genes by variance divided by the fitted variance (with fewer than 10 such genes it ranks by variance)
- `"cell_ranger"`: ranks genes by the squared coefficient of variation (variance / mean²)
- `"seurat_v3"`: the variance-stabilising selection of Stuart et al. (2019) as scanpy defines it, for **raw counts** in `data.X`. The log10 variance of the genes with non-zero variance is regressed on their log10 mean by a local regression of degree 2 with span 0.3; the fitted value gives the expected standard deviation of a gene of that mean; the counts are clipped at `sd × √(number of cells) + mean`; and the genes are ranked by the variance of the clipped counts in units of the fitted standard deviation (ties keep gene order). That standardised variance is stored in `data.var['variances_norm']`
- `"variance"`: ranks genes by variance

**`seurat_v3` requirements.** It raises `ValueError` for fewer than two cells, for negative or non-finite values, for fewer than 17 genes with non-zero variance (too few for the local regression), and when `int(0.3 × n)` or more of those `n` genes have exactly the same mean (a neighbourhood of zero width: filter out rarely detected genes first). Non-integer values give a `UserWarning` and are used as they are. The exact local regression takes O(genes²) time. Scanpy gets its smoother from scikit-misc, which interpolates between kd-tree vertices, whereas this one is evaluated exactly at every gene, so genes right at the `n_top_genes` cut-off can rank differently from scanpy; the size of that difference has not been measured (scikit-misc is not installed here).

The function only flags genes; it does not subset the matrix. There are no `min_mean`, `max_mean` or `min_disp` thresholds.

**Adds to data:**

- `data.var['highly_variable']`: boolean mask of the selected genes
- `data.var['means']`, `data.var['variances']`: mean and variance per gene
- `data.var['variances_norm']`: depends on the flavor. For `"seurat_v3"` it is the standardised variance that ranks the genes; for the other flavors it is the z-score of the variance across genes
- `data.var['dispersions']`: squared coefficient of variation per gene
- `data.uns['hvg']`: `flavor`, `n_top_genes` and `n_selected`

`metainformant.singlecell.data.preprocessing.identify_highly_variable_genes(data, n_top_genes=2000, flavor="seurat")` is an alternative that also writes `data.var['highly_variable']` (flavors `"seurat"` and `"cell_ranger"`, other `var` columns and `uns['highly_variable_genes']`).

## Principal Component Analysis (PCA)

### compute_pca()

```python
from metainformant.singlecell.analysis.pca_methods import compute_pca

# PCA on all genes
data = compute_pca(data, n_components=50, random_state=42)

# PCA on the genes flagged by select_hvgs()
data = compute_pca(data, n_components=50, random_state=42, use_hvgs=True)

print(f"PC loadings shape: {data.varm['PCs'].shape}")           # genes x components
print(f"Explained variance: {data.uns['pca']['explained_variance_ratio'][:5]}")
```

**Parameters:**

- `n_components`: number of components (default 50); clamped to `min(n_obs - 1, genes used)`
- `random_state`: seed for the PCA solver
- `scale_data`: standardise every gene before the PCA (default `True`)
- `use_hvgs`: keyword-only; restrict the PCA to the genes flagged in `data.var['highly_variable']`. It raises `ValueError` when that column is missing, not boolean (or 0/1), has missing values or flags no gene. Scaling, the variance ratios and the component clamp then refer to those genes only, and `varm['PCs']` keeps one row per gene with exact zeros for the genes that were left out

`pca_reduction` is the same function, and `run_pca` is an alias that also accepts a raw `expression_matrix=` and then returns a dictionary.

**Adds to data:**

- `data.obsm['X_pca']` and the `obs` columns `PC1`, `PC2`, ...: cell coordinates
- `data.varm['PCs']`: gene loadings (genes × components)
- `data.uns['pca']`: `n_components`, `explained_variance_ratio`, `cumulative_explained_variance`, `singular_values`, `components_shape`, `random_state`, `scaled`, `use_hvgs` and `n_genes_used`. The keys `variance_ratio` and `variance` are aliases of `explained_variance_ratio` (both hold the ratio), and `cumulative_variance` is an alias of `cumulative_explained_variance`

### Choosing the Number of Components

```python
import matplotlib.pyplot as plt
import numpy as np

variance_ratio = data.uns["pca"]["explained_variance_ratio"]
cumulative = np.cumsum(variance_ratio)

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(range(1, len(variance_ratio) + 1), variance_ratio, "bo-")
plt.xlabel("Principal Component")
plt.ylabel("Explained Variance Ratio")
plt.title("Scree Plot")

plt.subplot(1, 2, 2)
plt.plot(range(1, len(cumulative) + 1), cumulative, "ro-")
plt.axhline(y=0.8, color="k", linestyle="--", alpha=0.5)
plt.xlabel("Principal Component")
plt.ylabel("Cumulative Explained Variance")
plt.show()
```

### Other linear reductions

| Function | Adds |
|----------|------|
| `ica_reduction(data, n_components=10, random_state=None, max_iter=1000)` | `obs` columns `IC1`, ..., `varm['ICs']`, `uns['ica']` |
| `factor_analysis_reduction(data, n_components=10, random_state=None)` | `obs` columns `FA1`, ..., `varm['FA_loadings']`, `uns['factor_analysis']` |

## Neighbor Graph Construction

### compute_neighbors()

Build the k-nearest-neighbor graph that diffusion maps and diffusion pseudotime use:

```python
from metainformant.singlecell.analysis.nonlinear_methods import compute_neighbors

data = compute_neighbors(data, n_neighbors=15, n_pcs=40)   # on obsm['X_pca'] when it exists, else on data.X

# Write the graph into `data` itself instead of a copy
compute_neighbors(data, n_neighbors=15, inplace=True)
```

**Parameters:**

- `n_neighbors`: neighbors per cell (default 15; the row of each cell also contains the cell itself, as its first entry)
- `metric`: distance metric passed to scikit-learn (default `"euclidean"`)
- `use_rep`: with `"X_pca"` (the default) the graph is built on `data.obsm['X_pca']` when `compute_pca` has been run, and on `data.X` otherwise. Any other value also uses `data.X`; no other `obsm` key is looked up
- `n_pcs`: keyword-only; use only the first `n_pcs` principal components
- `method`: keyword-only; only `"sklearn"` is supported, anything else raises `ValueError`
- `inplace`: keyword-only, default `False`

There is no `random_state` argument: the neighbor search is exact and deterministic.

**Adds to data:**

`data.uns['neighbors']`, a dictionary with `indices` and `distances` (cells × `n_neighbors` arrays), `connectivities` (a sparse cells × cells matrix of Gaussian-kernel affinities with the median distance as bandwidth, symmetrised by averaging) and `params`. Nothing is written to `data.obsp`.

The graph is indexed by cell. `filter_cells` does not update it, so recompute it after filtering the cells; `compute_diffusion_pseudotime` raises an error on a stored graph with the wrong number of cells.

## UMAP

### compute_umap()

Uniform Manifold Approximation and Projection:

```python
from metainformant.singlecell.analysis.nonlinear_methods import compute_umap

data = compute_umap(
    data,
    n_components=2,        # 1 to 3
    n_neighbors=15,        # 5 to 100
    min_dist=0.5,          # 0 to 1
    random_state=42,
)

umap_coords = data.obsm["X_umap"]    # also in data.obs as UMAP1, UMAP2
```

**Parameters:** `n_components` (1–3), `n_neighbors` (5–100), `min_dist` (0–1), `random_state`, `metric` (default `"euclidean"`) and the keyword-only `n_epochs`. Values outside the ranges raise `ValidationError`. There is no `spread` parameter.

UMAP is fitted on the expression matrix `data.X`, not on the PCA coordinates. When `data.uns['neighbors']` is missing the returned copy also gets a neighbor graph, which UMAP itself does not use. `data.uns['umap']` records the parameters. `umap_reduction` is the same function, and `run_umap` is a wrapper of it without the `n_epochs` argument.

**Note:** Requires the optional `umap-learn` package (installed by the `singlecell` extra); without it the call raises `ImportError`:

```bash
uv pip install umap-learn
```

## t-SNE

### compute_tsne()

t-distributed Stochastic Neighbor Embedding:

```python
from metainformant.singlecell.analysis.nonlinear_methods import compute_tsne

data = compute_tsne(
    data,
    n_components=2,
    perplexity=30,
    max_iter=1000,
    random_state=42,
)

tsne_coords = data.obsm["X_tsne"]    # also in data.obs as tSNE1, tSNE2
```

**Parameters:** `n_components`, `perplexity`, `random_state`, `learning_rate` (default 200), `max_iter` (default 1000) and the keyword-only `n_iter` (alias of `max_iter`).

t-SNE is fitted on `data.X`; with more than 5000 cells the matrix is first reduced to up to 50 principal components. The perplexity is lowered with a warning when it is not below `(n_cells - 1) / 3`, and `max_iter` is raised to at least 250. `data.uns['tsne']` records the values actually used (`perplexity`, `max_iter`) together with `kl_divergence` and `n_iter`. `tsne_reduction` is the same function without the `n_iter` alias, and `run_tsne` is a wrapper of it.

## Diffusion Maps

### compute_diffusion_map()

Diffusion map coordinates from the neighbor graph:

```python
from metainformant.singlecell.analysis.nonlinear_methods import compute_diffusion_map

data = compute_diffusion_map(data, n_components=15, alpha=1.0)

dc_coords = data.obsm["X_diffmap"]            # cells x n_components
eigenvalues = data.uns["diffmap"]["eigenvalues"]
```

**Parameters:**

- `n_components`: number of diffusion components (default 10); the first, constant eigenvector is skipped
- `n_neighbors`: neighbors for the graph, used only when `data.uns['neighbors']` does not exist yet (the graph that gets built is stored on the result)
- `alpha`: normalisation of the operator; `0` uses the unnormalised adjacency, `1` the symmetric normalisation, other values the generalised `D^-alpha A D^-alpha` form
- `inplace`: keyword-only, default `False`; writes `obsm['X_diffmap']`, `uns['diffmap']` (and the graph, when it had to be built) into `data`

The edge weights are `1 / (1 + distance)`. The leading eigenvectors come from ARPACK with a fixed start vector (a dense eigendecomposition is the fallback when ARPACK fails), so repeated calls give identical results.

**Adds to data:**

- `data.obsm['X_diffmap']`: diffusion coordinates
- `data.uns['diffmap']`: `eigenvalues` and `params` (`n_components`, `n_neighbors`, `alpha`)

### diffusion_map_reduction()

A separate dense implementation that does not use the neighbor graph: a Gaussian kernel over all pairwise distances with the median distance as bandwidth. `diffusion_map_reduction(data, n_components=10, n_neighbors=15, alpha=1.0)` validates `n_components` (2–50), `n_neighbors` (5–100) and `alpha` (0–1), writes the coordinates to `obs` as `DC1`, `DC2`, ... and the parameters to `uns['diffusion_map']`; `n_neighbors` is recorded but not used by the computation. A failure of the eigendecomposition propagates as `numpy.linalg.LinAlgError`.

## Other Embeddings and Metrics

### mds_reduction()

`mds_reduction(data, n_components=2, metric=True, random_state=None)` runs multidimensional scaling on the cell-by-cell Euclidean distances (`n_components` 2 or 3). It writes `obs` columns `MDS1`, ... and `uns['mds']` (including the stress).

### compute_dimensionality_metrics()

Quality metrics of an embedding stored as `obs` columns:

```python
from metainformant.singlecell.analysis.nonlinear_methods import compute_dimensionality_metrics

metrics = compute_dimensionality_metrics(data, ["UMAP1", "UMAP2"], random_state=0)
print(metrics["distance_correlation"])
# Also: embedding_dimensions, n_samples, n_features and, up to 1000 cells,
# trustworthiness, continuity and neighborhood_size
```

The distance correlation compares the original and embedded pairwise distances of at most 1000 cells. `random_state` seeds that subsample; it only matters above 1000 cells, and `None` (the default) draws a different subsample on each call.

## Parameter Guidelines

### HVG Selection

- **2000-3000 genes**: good balance for most datasets
- **More genes (5000+)**: for very heterogeneous datasets
- **Fewer genes (1000)**: for homogeneous populations

### PCA

- **50 components**: standard for most analyses
- **Use the scree plot** to determine the optimal number
- **80% variance**: common cutoff for cumulative variance

### UMAP Parameters

- **n_neighbors**: 15-50 for most datasets; higher for smoother embeddings, lower for more local structure
- **min_dist**: 0.1-0.5 for tight clusters, 0.5-1.0 for a spread-out visualisation

### t-SNE Parameters

- **perplexity**: 30-50 for most datasets; check several values for robustness

## Common Workflows

### Standard Dimensionality Reduction Pipeline

```python
from metainformant.singlecell.analysis.pca_methods import compute_pca, select_hvgs
from metainformant.singlecell.analysis.nonlinear_methods import (
    compute_neighbors,
    compute_tsne,
    compute_umap,
)

# 1. Flag highly variable genes
data = select_hvgs(data, n_top_genes=2000, flavor="seurat")
print(f"Selected {data.var['highly_variable'].sum()} HVGs")

# 2. PCA on the flagged genes
data = compute_pca(data, n_components=50, use_hvgs=True)
print(f"Explained variance (first 10 PCs): {sum(data.uns['pca']['explained_variance_ratio'][:10]):.3f}")

# 3. Neighbor graph on the first 40 PCs
data = compute_neighbors(data, n_neighbors=15, n_pcs=40)

# 4. UMAP and t-SNE embeddings (both are fitted on data.X)
data = compute_umap(data, n_components=2, min_dist=0.5)
data = compute_tsne(data, n_components=2, perplexity=30)

print(f"Available embeddings: {list(data.obsm.keys())}")
```

### Trajectory Analysis Preparation

Diffusion pseudotime reads the neighbor graph in `data.uns['neighbors']` (and builds one from `data.X` when there is none), so the graph is the only preparation it needs; the diffusion map is optional:

```python
from metainformant.singlecell.analysis.nonlinear_methods import compute_diffusion_map, compute_neighbors
from metainformant.singlecell.analysis.pca_methods import compute_pca
from metainformant.singlecell.analysis.trajectory import compute_diffusion_pseudotime

data = compute_pca(data, n_components=50)
data = compute_neighbors(data, n_neighbors=30, n_pcs=30)     # more neighbors for trajectories
data = compute_diffusion_map(data, n_components=15)          # optional: coordinates for plotting
data = compute_diffusion_pseudotime(data, root_cell=0)
```

## Quality Assessment

Compare the original and embedded geometry with `compute_dimensionality_metrics` (above). For a quick global check, correlate pairwise distances yourself:

```python
import numpy as np
from scipy.spatial.distance import pdist

def global_preservation(X_orig, X_embed, n=1000):
    """Correlation of pairwise distances on the first n cells."""
    return np.corrcoef(pdist(X_orig[:n]), pdist(X_embed[:n]))[0, 1]
```

## Error Handling

```python
from metainformant.singlecell.analysis.nonlinear_methods import compute_tsne, compute_umap

try:
    data = compute_umap(data)
except ImportError:
    print("umap-learn not installed. Install with: uv pip install umap-learn")
    data = compute_tsne(data)    # fall back to t-SNE
```

`compute_pca(..., use_hvgs=True)` raises `ValueError` when `select_hvgs` has not been run, and the range checks of `compute_umap` and `diffusion_map_reduction` raise `ValidationError` (from `metainformant.core.utils.errors`).

## Integration with Clustering

The clustering functions do not read PCA coordinates or the neighbor graph: each builds its own kNN graph from the expression matrix of the object it is given. To cluster on the principal components, pass them as a matrix (the function then returns the labels; see [Clustering](./clustering.md) for the Leiden requirements):

```python
from metainformant.singlecell.analysis.clustering import leiden_clustering

# Cluster the cells on their first 30 principal components
data.obs["leiden_pca"] = leiden_clustering(data.obsm["X_pca"][:, :30], resolution=0.5, random_state=0)

# Cluster on the expression matrix itself (a SingleCellData in, a SingleCellData out)
data = leiden_clustering(data, resolution=0.5)
```

## Testing

Dimensionality-reduction tests are in `tests/singlecell/test_singlecell_dimensionality.py`; the `inplace` contract is pinned by `tests/singlecell/test_singlecell_input_contracts.py` and the seeding by `tests/singlecell/test_singlecell_randomness.py`:

```bash
uv run pytest tests/singlecell/test_singlecell_dimensionality.py tests/singlecell/test_singlecell_input_contracts.py -v
```

## Related Documentation

- [Preprocessing](./preprocessing.md): For data preparation before dimensionality reduction
- [Clustering](./clustering.md): For clustering the cells
- [Visualization](./visualization.md): For plotting embeddings
- [Trajectory](./trajectory.md): For pseudotime analysis on the neighbor graph
