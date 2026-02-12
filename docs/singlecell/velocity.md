# RNA Velocity

The velocity module estimates RNA velocity from spliced and unspliced count
matrices, projects velocity vectors into low-dimensional embeddings, derives
pseudotime orderings, evaluates confidence metrics, and fits full dynamical
kinetic models.

## Concepts

### RNA Velocity Theory

RNA velocity exploits the ratio of unspliced (nascent) to spliced (mature) mRNA
to infer each cell's future transcriptional state. Under steady-state:
`v_ig = u_ig - gamma_g * s_ig`, where `u` is unspliced count, `s` is spliced
count, and `gamma` is the degradation rate from regressing unspliced on spliced.

### Steady-State Model

`compute_velocity` fits a no-intercept linear regression per gene to estimate
gamma. Genes are filtered by `min_counts` and `r_squared_threshold`.

### Dynamical Model

`fit_dynamical_model` estimates alpha (transcription), beta (splicing), and
gamma (degradation) per gene via an iterative EM-like procedure initialized from
the steady-state gamma. Convergence is checked against a relative tolerance.

### Velocity Embedding

`velocity_embedding` projects velocity vectors into a 2D embedding (UMAP/tSNE).
Cosine similarity between each cell's velocity and expression displacements to
neighbors yields transition probabilities; the projected velocity is their
weighted sum in embedding space.

### Velocity Pseudotime

`velocity_pseudotime` derives temporal ordering from the transition matrix. Root
defaults to the cell with smallest outgoing velocity. Normalized to [0, 1].

### Confidence Metrics

`velocity_confidence` reports per-gene confidence (consistent velocity sign
fraction) and per-cell confidence (mean absolute velocity / expression level).

## Function Reference

### `compute_velocity(spliced, unspliced, gene_names, method, min_counts, r_squared_threshold) -> dict`

Estimate velocity with the steady-state linear model. Returns `velocity_matrix`,
`gamma`, `r_squared`, `velocity_genes`, and `n_velocity_genes`.

| Parameter             | Type        | Default          | Description                        |
|-----------------------|-------------|------------------|------------------------------------|
| `spliced`             | Any         | --               | Spliced counts (cells x genes)     |
| `unspliced`           | Any         | --               | Unspliced counts (cells x genes)   |
| `gene_names`          | `list[str]` | --               | Gene names for columns             |
| `method`              | `str`       | `"steady_state"` | Estimation method                  |
| `min_counts`          | `int`       | `10`             | Min total counts per gene          |
| `r_squared_threshold` | `float`     | `0.01`           | Min R-squared for velocity genes   |

### `velocity_embedding(velocity, embedding, n_neighbors) -> dict`

Project velocity into embedding space. Returns `velocity_embedding`,
`transition_matrix`, and `cell_velocities`.

| Parameter     | Type  | Default | Description                                   |
|---------------|-------|---------|-----------------------------------------------|
| `velocity`    | Any   | --      | Velocity matrix from `compute_velocity`       |
| `embedding`   | Any   | --      | 2D embedding coordinates (cells x 2)          |
| `n_neighbors` | `int` | `30`    | Neighbors for transition computation           |

### `velocity_pseudotime(velocity, embedding, root_cell, n_neighbors) -> dict`

Derive pseudotime from velocity transitions. Returns `pseudotime`, `root_cell`,
and `terminal_cells`.

| Parameter     | Type          | Default | Description                            |
|---------------|---------------|---------|----------------------------------------|
| `velocity`    | Any           | --      | Velocity matrix                        |
| `embedding`   | Any           | --      | Embedding coordinates                  |
| `root_cell`   | `int \| None` | `None`  | Root cell index (auto-selected if None)|
| `n_neighbors` | `int`         | `30`    | Neighbors for transitions              |

### `velocity_confidence(velocity, spliced) -> dict`

Evaluate velocity reliability. Returns `gene_confidence`, `cell_confidence`,
`overall_confidence`, `n_high_confidence_genes`, `n_high_confidence_cells`.

| Parameter  | Type | Default | Description                              |
|------------|------|---------|------------------------------------------|
| `velocity` | Any  | --      | Velocity matrix                          |
| `spliced`  | Any  | --      | Spliced counts for expression weighting  |

### `fit_dynamical_model(spliced, unspliced, gene_names, max_iter, tol, min_counts) -> dict`

Fit full kinetic parameters. Returns `alpha`, `beta`, `gamma`, `likelihood`,
`velocity_matrix`, `fitted_genes`, and `n_iterations`.

| Parameter    | Type        | Default | Description                             |
|--------------|-------------|---------|-----------------------------------------|
| `spliced`    | Any         | --      | Spliced counts (cells x genes)          |
| `unspliced`  | Any         | --      | Unspliced counts (cells x genes)        |
| `gene_names` | `list[str]` | --      | Gene names for columns                  |
| `max_iter`   | `int`       | `100`   | Maximum optimization iterations         |
| `tol`        | `float`     | `1e-4`  | Convergence tolerance                   |
| `min_counts` | `int`       | `10`    | Min total counts per gene               |

## Code Examples

```python
from metainformant.singlecell.velocity.rna_velocity import (
    compute_velocity, velocity_embedding, velocity_pseudotime,
    velocity_confidence, fit_dynamical_model,
)

# Count matrices (4 cells, 3 genes)
spliced = [[10, 5, 20], [8, 12, 15], [15, 3, 25], [6, 18, 10]]
unspliced = [[3, 2, 5], [2, 4, 4], [4, 1, 7], [1, 6, 3]]
genes = ["GAPDH", "TP53", "MYC"]

# Steady-state velocity
vel = compute_velocity(spliced, unspliced, genes)
print(f"Velocity genes: {vel['velocity_genes']}")

# Project onto UMAP and derive pseudotime
embedding = [[1.0, 2.0], [3.0, 1.5], [2.5, 4.0], [4.0, 3.0]]
emb = velocity_embedding(vel["velocity_matrix"], embedding, n_neighbors=3)
pt = velocity_pseudotime(vel["velocity_matrix"], embedding, n_neighbors=3)
print(f"Root: {pt['root_cell']}, terminals: {pt['terminal_cells']}")

# Confidence and dynamical model
conf = velocity_confidence(vel["velocity_matrix"], spliced)
print(f"Overall confidence: {conf['overall_confidence']:.3f}")
dyn = fit_dynamical_model(spliced, unspliced, genes, max_iter=50)
print(f"Fitted genes: {dyn['fitted_genes']}")
```

## Configuration

All parameters are function arguments. No external config files required.

## Import Path

```python
from metainformant.singlecell.velocity.rna_velocity import compute_velocity
from metainformant.singlecell.velocity.rna_velocity import velocity_embedding
from metainformant.singlecell.velocity.rna_velocity import velocity_pseudotime
from metainformant.singlecell.velocity.rna_velocity import velocity_confidence
from metainformant.singlecell.velocity.rna_velocity import fit_dynamical_model
```
