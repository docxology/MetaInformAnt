# Single-Cell Trajectory Analysis

`metainformant.singlecell.analysis.trajectory` orders cells along a developmental axis. It offers one graph-based method that follows the published algorithm (diffusion pseudotime) and three simplified ones (an embedding distance, a cluster-level graph and cluster-level lineages). The functions that compute a pseudotime return a modified copy and leave their input untouched.

| Function | Pseudotime column in `data.obs` | Results in `data.uns` | Method |
|----------|--------------------------------|-----------------------|--------|
| `compute_diffusion_pseudotime`, `dpt_trajectory` | `dpt_pseudotime` | `dpt` | Diffusion pseudotime (Haghverdi et al. 2016) on the neighbor graph |
| `compute_pseudotime_from_dimensionality_reduction` | `embedding_pseudotime` | `embedding_pseudotime` | Euclidean distance from a root cell in an embedding |
| `paga_trajectory` | `paga_pseudotime` | `paga` | Approximation: group connectivity from the row order of the cells |
| `slingshot_trajectory` | `slingshot_pseudotime` | `slingshot` | Approximation: straight lines between cluster centroids |

## Pseudotime Inference

### compute_diffusion_pseudotime()

Diffusion pseudotime (DPT) of Haghverdi et al. (2016, Nature Methods 13:845):

```python
from metainformant.singlecell.analysis.trajectory import compute_diffusion_pseudotime

data = compute_diffusion_pseudotime(
    data,
    root_cell=0,        # Row index of the root cell (a heuristic picks one when None)
    n_components=10,    # Number of non-trivial diffusion components in the distance
    n_neighbors=15,     # Only used when data.uns['neighbors'] does not exist
)

pseudotime = data.obs["dpt_pseudotime"]
print(data.uns["dpt"]["root_cell"], data.uns["dpt"]["n_unreachable"])
```

**How it works:**

1. **Graph.** The cells are connected by a k-nearest-neighbor graph with Gaussian affinities. The function uses `data.uns['neighbors']` when it exists (see `compute_neighbors` in [Dimensionality Reduction](./dimensionality.md); a graph built on the PCA coordinates is used as is) and otherwise builds one from `data.X` with `n_neighbors` neighbors and stores it in the returned object. It never reads `data.obsm['X_pca']` itself.
2. **Transition matrix.** The affinities are density-normalised into a random-walk transition matrix `T`.
3. **Distance.** The DPT distance between two cells is the Euclidean distance between their rows of the accumulated transition matrix `sum over t >= 1 of (T^t - steady state)`, evaluated on the leading `n_components` non-trivial eigenvectors of `T`. The eigenproblem is solved on the cells connected to the root: with a dense solver up to 1000 of them, and with ARPACK from a fixed start vector above that, so repeated calls give identical values.
4. **Scaling.** The pseudotime of a cell is its DPT distance from the root divided by the largest such distance: it spans [0, 1] with the root at 0.

Because the distance follows the graph, pseudotime keeps increasing along a curved trajectory, where the straight-line distance from the root would turn back. The function computes the pseudotime from one root only; the branching detection of the paper is not part of it.

**Parameters:**

- `root_cell`: integer row index (0 to `n_obs - 1`) of the starting cell. When `None`, a heuristic picks the cell with the largest gene-variance-weighted expression sum. DPT measures distance from the root, so the direction of the pseudotime is only as meaningful as this choice; the chosen cell is recorded in `uns['dpt']['root_cell']` and `uns['dpt']['root_selection']` (`"user"` or `"heuristic"`).
- `n_components`: number of non-trivial diffusion components (eigenvectors) that enter the distance (default 10).
- `n_neighbors`: keyword-only, default 15. It builds the graph only when `data.uns['neighbors']` is absent; an existing graph is used with its own number of neighbors.

**Cells the graph does not connect to the root** have no finite distance to it. They get `inf` in `dpt_pseudotime`, a warning is logged, and `uns['dpt']['n_unreachable']` counts them. Filter with `np.isfinite` before plotting or summarising the column.

**Errors.** `ValidationError` (from `metainformant.core.utils.errors`) is raised for an invalid `root_cell`, `n_components` or `n_neighbors`, for an argument that is not a `SingleCellData`, and for a stored neighbor graph that is stale: `uns['neighbors']` describing a different number of cells than `data` (for example after `filter_cells`, which does not update it), lacking `connectivities` or holding negative affinities. Recompute the graph with `compute_neighbors` in the first case.

**Adds to data:**

- `data.obs['dpt_pseudotime']`: pseudotime per cell
- `data.uns['dpt']`: `method`, `root_cell`, `root_selection`, `n_components` (eigenvectors actually used), `eigenvalues`, `n_neighbors` (of the graph that was used), `neighbors_source` (`"uns"` or `"X"`), `n_unreachable`, `min_pseudotime` and `max_pseudotime` (over the finite values)

### dpt_trajectory()

`dpt_trajectory(data, root_cell=None)` runs `compute_diffusion_pseudotime` with its defaults of 10 components and 15 neighbors.

### compute_pseudotime_from_dimensionality_reduction()

Pseudotime from the coordinates of an embedding stored as `obs` columns:

```python
from metainformant.singlecell.analysis.trajectory import compute_pseudotime_from_dimensionality_reduction

data = compute_pseudotime_from_dimensionality_reduction(
    data,
    dim_red_cols=["UMAP1", "UMAP2"],   # obs columns written by compute_umap
    root_cell=0,                       # default: the cell closest to the embedding origin
)
pseudotime = data.obs["embedding_pseudotime"]
```

The pseudotime is the straight-line (Euclidean) distance from the root cell in the embedding, scaled to [0, 1]. It is not a geodesic along the data manifold, so it turns back on a curved trajectory; use `compute_diffusion_pseudotime` for a graph-based pseudotime. `uns['embedding_pseudotime']` records the columns, the root cell and the range.

## Cluster-Level Trajectories (Approximations)

Both functions below work on cluster labels and are simplified approximations of PAGA and Slingshot, not implementations of them.

### paga_trajectory()

```python
from metainformant.singlecell.analysis.trajectory import paga_trajectory

data = paga_trajectory(data, groups="leiden_cluster")
connectivity = data.uns["paga"]["connectivity_matrix"]     # groups x groups, rows normalised
pseudotime = data.obs["paga_pseudotime"]                   # every cell takes its group's value
```

- The connectivity matrix counts how often a cell of one group is directly followed by a cell of another group **in the row order of the data**, normalised per row. It is not computed from the expression values or from the neighbor graph, so it depends on how the rows happen to be ordered: shuffling the cells gives a different graph.
- The group pseudotime is the number of hops from the group with the fewest incoming transitions, scaled to [0, 1].
- `uns['paga']` holds `groups`, `unique_groups`, `connectivity_matrix` and `n_groups`. A missing `groups` column raises `ValidationError`.

### slingshot_trajectory()

```python
from metainformant.singlecell.analysis.trajectory import slingshot_trajectory

data = slingshot_trajectory(
    data,
    start_cluster=0,          # Root cluster (same type as the labels in obs)
    end_clusters=[2, 4],      # Terminal clusters, one branch per end cluster
)
pseudotime = data.obs["slingshot_pseudotime"]
branches = data.uns["slingshot"]["trajectory_info"]["branches"]
```

- The cluster column is detected automatically: the first `obs` column whose name contains `cluster`, `leiden` or `louvain` (case-insensitive). None found, or a start or end cluster that is not among its values, raises `ValidationError`.
- Each branch is a straight line between the start and an end cluster: `path` is just `[start, end]` and `length` the Euclidean distance between their centroids in gene space. There is no minimum spanning tree, no principal curve and no smoothing.
- The pseudotime of a cluster is the distance of its centroid from the start centroid divided by the largest such distance, and every cell takes the value of its cluster.
- `uns['slingshot']` holds `start_cluster`, `end_clusters` and `trajectory_info` (`branches`, `cluster_pseudotime`, `centroids`).

## Lineage Analysis

### find_trajectory_branches()

Split cells that are ordered by an existing pseudotime column into consecutive groups:

```python
from metainformant.singlecell.analysis.trajectory import find_trajectory_branches

branches = find_trajectory_branches(
    data,
    pseudotime_col="dpt_pseudotime",
    min_branch_size=10,        # Minimum number of cells before a group may be closed
)

# Add the group of each cell to data (-1 = cell not assigned to a group)
data.obs["branch"] = branches["branch_labels"]
print(f"Found {branches['n_branches']} groups, sizes: {branches['branch_sizes']}")
```

This is a simple heuristic: it sorts the cells by pseudotime and closes a group at any gap to the next pseudotime value that exceeds the 90th percentile of all gaps, once the group holds at least `min_branch_size` cells. It detects gaps along the pseudotime axis, not bifurcations of the graph. Remove cells with `inf` pseudotime first (they make it return a single group).

### compute_trajectory_entropy()

Gene-expression entropy in a sliding window along the pseudotime axis:

```python
from metainformant.singlecell.analysis.trajectory import compute_trajectory_entropy

result = compute_trajectory_entropy(data, pseudotime_col="dpt_pseudotime", window_size=100)
# entropy_profile, pseudotime_windows, differentiation_points (local maxima), window_size, n_windows
```

`window_size` must be between 1 and the number of cells, otherwise `ValidationError` is raised.

## Visualization

### Trajectory Visualization

`plot_trajectory(data, trajectory_key)` colours the cells of an embedding by pseudotime; `trajectory_key` names the entry of `data.uns` to draw (`"dpt"` also marks the root cell):

```python
import matplotlib.pyplot as plt
from metainformant.singlecell.visualization.visualization import plot_trajectory

fig = plot_trajectory(data, "dpt")        # needs UMAP, t-SNE or PCA columns in data.obs
plt.show()
```

Plotting pseudotime on a stored embedding and a gene against pseudotime by hand:

```python
import matplotlib.pyplot as plt
import numpy as np

pseudotime = data.obs["dpt_pseudotime"].to_numpy()
finite = np.isfinite(pseudotime)

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
umap = data.obsm["X_umap"]
scatter = plt.scatter(umap[finite, 0], umap[finite, 1], c=pseudotime[finite], cmap="viridis", alpha=0.7)
plt.colorbar(scatter, label="Pseudotime")
plt.xlabel("UMAP 1")
plt.ylabel("UMAP 2")

plt.subplot(1, 2, 2)
gene = "GENE1"  # Replace with an actual gene
if gene in data.var.index:
    column = data.var.index.get_loc(gene)
    values = data.X[:, column]
    values = values.toarray().ravel() if hasattr(values, "toarray") else np.asarray(values).ravel()
    plt.scatter(pseudotime[finite], values[finite], alpha=0.5)
    plt.xlabel("Pseudotime")
    plt.ylabel(f"{gene} expression")

plt.tight_layout()
plt.show()
```

## Integration with Other Modules

### Preprocessing for Trajectory Analysis

```python
from metainformant.singlecell.data.preprocessing import (
    calculate_qc_metrics,
    filter_cells,
    log_transform,
    normalize_counts,
)
from metainformant.singlecell.analysis.pca_methods import compute_pca, select_hvgs
from metainformant.singlecell.analysis.nonlinear_methods import compute_diffusion_map, compute_neighbors
from metainformant.singlecell.analysis.trajectory import compute_diffusion_pseudotime

# 1. Standard preprocessing
data = calculate_qc_metrics(data)
data = filter_cells(data, min_genes=200, max_pct_mt=20)
data = normalize_counts(data, method="total_count")
data = log_transform(data)

# 2. Feature selection and dimensionality reduction
data = select_hvgs(data, n_top_genes=2000)
data = compute_pca(data, n_components=50, use_hvgs=True)

# 3. Neighbor graph (after filtering, so that it matches the remaining cells)
data = compute_neighbors(data, n_neighbors=30, n_pcs=30)      # more neighbors for trajectories
data = compute_diffusion_map(data, n_components=15)           # optional: coordinates for plotting

# 4. Trajectory inference on that graph
data = compute_diffusion_pseudotime(data, root_cell=0)
```

### Using several methods

Each method adds its own pseudotime column, so they can be run on the same object:

```python
from metainformant.singlecell.analysis.trajectory import (
    compute_diffusion_pseudotime,
    compute_pseudotime_from_dimensionality_reduction,
    paga_trajectory,
    slingshot_trajectory,
)

data = compute_diffusion_pseudotime(data, root_cell=0)                    # dpt_pseudotime
data = paga_trajectory(data, groups="leiden_cluster")                     # paga_pseudotime
data = slingshot_trajectory(data, start_cluster=0, end_clusters=[2])      # slingshot_pseudotime
data = compute_pseudotime_from_dimensionality_reduction(
    data, dim_red_cols=["UMAP1", "UMAP2"], root_cell=0                    # embedding_pseudotime
)
```

### Method Comparison

```python
import numpy as np
from scipy.stats import spearmanr

dpt = data.obs["dpt_pseudotime"].to_numpy()
embedding = data.obs["embedding_pseudotime"].to_numpy()
finite = np.isfinite(dpt)
print(spearmanr(dpt[finite], embedding[finite]))
```

## Choosing the Number of Components

`n_components` sets how many non-trivial diffusion components enter the DPT distance. The default of 10 suits most data; fewer components keep only the coarsest structure, more components let finer structure contribute. Check that the pseudotime still orders cells sensibly (for example against known marker genes) when you change it.

## Quality Control

### Trajectory Quality Assessment

```python
import matplotlib.pyplot as plt
import numpy as np

def assess_trajectory_quality(data, cluster_col="leiden_cluster", pseudotime_col="dpt_pseudotime"):
    """Plot the pseudotime distribution, overall and per cluster."""
    pseudotime = data.obs[pseudotime_col].to_numpy()
    finite = np.isfinite(pseudotime)

    plt.figure(figsize=(10, 4))

    plt.subplot(1, 2, 1)
    plt.hist(pseudotime[finite], bins=50, alpha=0.7)
    plt.xlabel("Pseudotime")
    plt.ylabel("Number of cells")
    plt.title("Pseudotime distribution")

    plt.subplot(1, 2, 2)
    for cluster in data.obs[cluster_col].unique():
        mask = finite & (data.obs[cluster_col] == cluster).to_numpy()
        plt.hist(pseudotime[mask], alpha=0.5, label=f"Cluster {cluster}")
    plt.xlabel("Pseudotime")
    plt.legend()
    plt.title("Pseudotime by cluster")

    plt.tight_layout()
    plt.show()

assess_trajectory_quality(data)
```

## Applications

### Developmental Biology

- **Embryonic development**: Track cell fate decisions
- **Organogenesis**: Understand tissue formation
- **Stem cell differentiation**: Map lineage relationships

### Disease Studies

- **Cancer progression**: Tumor evolution trajectories
- **Drug response**: Time-course treatment effects
- **Immune responses**: T cell activation dynamics

## Best Practices

1. **Quality control**: Ensure good quality cells before trajectory analysis
2. **Choose the root deliberately**: pass `root_cell` from biology (a known progenitor); the heuristic default is only a convenience
3. **Check for unreachable cells**: `uns['dpt']['n_unreachable']` should normally be 0; raise `n_neighbors` or remove isolated cells otherwise
4. **Match the graph to the cells**: recompute `compute_neighbors` after `filter_cells`
5. **Multiple methods**: compare results from different algorithms, remembering that the PAGA and Slingshot functions are approximations
6. **Biological validation**: confirm with known developmental markers

## Limitations

- **Requires prior knowledge**: you must name the root cell (or accept a heuristic one)
- **Assumes continuity**: may not work for discrete state transitions
- **Sample coverage**: sparse sampling can miss intermediate states
- **Technical noise**: sensitive to technical artifacts
- **Simplified cluster methods**: `paga_trajectory` depends on the row order of the cells and `slingshot_trajectory` draws straight lines between centroids

## Testing

Trajectory analysis is covered by `tests/singlecell/test_singlecell_diffusion_pseudotime.py` (a curved trajectory, scale invariance, a brute-force check of the DPT definition and the ARPACK path) and `tests/singlecell/test_singlecell_clustering_trajectory.py`:

```bash
uv run pytest tests/singlecell/test_singlecell_diffusion_pseudotime.py tests/singlecell/test_singlecell_clustering_trajectory.py -v
```

## Related Documentation

- [Dimensionality Reduction](./dimensionality.md): For the neighbor graph, PCA and diffusion maps
- [Clustering](./clustering.md): For identifying cell populations
- [Visualization](./visualization.md): For trajectory plotting
- [Preprocessing](./preprocessing.md): For data preparation
