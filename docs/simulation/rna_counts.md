### Simulation: RNA Counts

The `metainformant.simulation.models.rna` module provides functions for simulating
RNA-seq count data using negative binomial distributions. It covers bulk RNA-seq,
single-cell RNA-seq (with dropout), differential expression, time-series expression,
spatial transcriptomics, and technical noise modeling.

---

## Negative Binomial Model

RNA-seq count data is modeled with the negative binomial (NB) distribution, which
accounts for both biological and technical overdispersion relative to a Poisson model.

The parameterization used:

- **Mean**: `mu` (average expression level per gene)
- **Dispersion**: `phi` (overdispersion parameter)
- **Variance**: `mu + mu^2 * phi`

Conversion to numpy's NB parameters:
- `n = 1 / phi` (number of successes)
- `p = 1 / (1 + mu * phi)` (success probability)

---

## Core Function

### `simulate_counts_negative_binomial`

The low-level function that all other simulators build upon. Requires explicit
mean and dispersion arrays.

```python
import numpy as np
from metainformant.simulation import simulate_counts_negative_binomial

means = np.array([100.0, 50.0, 200.0, 10.0, 500.0])
dispersions = np.array([0.1, 0.2, 0.05, 0.5, 0.1])

counts = simulate_counts_negative_binomial(
    n_samples=6,
    n_features=5,
    means=means,
    dispersions=dispersions,
)
# Shape: (6, 5) -- 6 samples x 5 genes
```

**Parameters:**

| Parameter     | Type         | Description                                           |
|---------------|--------------|-------------------------------------------------------|
| `n_samples`   | `int`        | Number of samples/cells to simulate (>= 1)            |
| `n_features`  | `int`        | Number of genes/features (>= 1)                       |
| `means`       | `np.ndarray` | Mean expression per feature, shape `(n_features,)`    |
| `dispersions` | `np.ndarray` | Dispersion per feature, shape `(n_features,)`, > 0    |
| `rng`         | `Random`     | Optional random number generator for reproducibility  |

---

## Convenience Wrapper

### `simulate_rnaseq_counts`

A user-friendly wrapper that auto-generates realistic mean and dispersion arrays
from simple parameters. Gene means follow a log-normal distribution.

```python
from metainformant.simulation import simulate_rnaseq_counts

counts = simulate_rnaseq_counts(n_genes=1000, n_samples=10)
# Shape: (10, 1000)

# Customize expression level and dispersion
counts = simulate_rnaseq_counts(
    n_genes=500,
    n_samples=6,
    mean_expression=200.0,
    dispersion=0.3,
)
```

**Parameters:**

| Parameter          | Type    | Default  | Description                          |
|--------------------|---------|----------|--------------------------------------|
| `n_genes`          | `int`   | `1000`   | Number of genes to simulate          |
| `n_samples`        | `int`   | `10`     | Number of samples                    |
| `mean_expression`  | `float` | `100.0`  | Average expression level             |
| `dispersion`       | `float` | `0.5`    | NB dispersion parameter              |
| `rng`              | `Random`| `None`   | Random number generator              |

---

## Differential Expression

### `simulate_differential_expression`

Simulate two-group experimental designs with specified fold changes for a subset
of genes. Samples are split evenly between conditions.

```python
import numpy as np
from metainformant.simulation import simulate_differential_expression

fold_changes = np.array([2.0, -1.5, 3.0, 0.5])
expression, labels = simulate_differential_expression(
    n_samples=20,
    n_features=5000,
    fold_changes=fold_changes,
)
# expression.shape: (20, 5000)
# labels: array of 0s and 1s indicating condition
# First 4 genes are differentially expressed
```

**Parameters:**

| Parameter      | Type         | Description                                          |
|----------------|--------------|------------------------------------------------------|
| `n_samples`    | `int`        | Total samples across both conditions (>= 2)          |
| `n_features`   | `int`        | Total genes                                          |
| `fold_changes` | `np.ndarray` | Fold change values for DE genes (positive=up, negative=down) |
| `rng`          | `Random`     | Optional random number generator                     |

**Returns:** Tuple of `(expression_matrix, group_labels)`

---

## Bulk RNA-seq

### `simulate_bulk_rnaseq`

Simulate bulk RNA-seq with realistic library sizes and multinomial count generation.

```python
from metainformant.simulation import simulate_bulk_rnaseq

counts = simulate_bulk_rnaseq(n_samples=12, n_genes=2000)
# Library sizes follow log-normal distribution (~1M reads)
# Gene means follow power-law distribution
```

**Parameters:**

| Parameter       | Type              | Default | Description                          |
|-----------------|-------------------|---------|--------------------------------------|
| `n_samples`     | `int`             | required| Number of samples                    |
| `n_genes`       | `int`             | required| Number of genes                      |
| `library_sizes` | `np.ndarray|None` | `None`  | Custom library sizes per sample      |
| `gene_means`    | `np.ndarray|None` | `None`  | Custom mean expression per gene      |
| `rng`           | `Random`          | `None`  | Random number generator              |

---

## Single-Cell RNA-seq

### `simulate_single_cell_rnaseq`

Simulate scRNA-seq data with cell type-specific expression patterns and dropout events.
Each cell type has marker genes with elevated expression. Dropout is modeled as a
Bernoulli process applied after count generation.

```python
from metainformant.simulation import simulate_single_cell_rnaseq

expression, cell_types = simulate_single_cell_rnaseq(
    n_cells=500,
    n_genes=2000,
    n_cell_types=5,
    dropout_rate=0.3,
)
# expression.shape: (500, 2000)
# cell_types: array of cell type labels (0-4)
```

**Parameters:**

| Parameter      | Type    | Default | Description                              |
|----------------|---------|---------|------------------------------------------|
| `n_cells`      | `int`   | required| Number of cells                          |
| `n_genes`      | `int`   | required| Number of genes                          |
| `n_cell_types` | `int`   | `5`     | Number of distinct cell types            |
| `dropout_rate` | `float` | `0.3`   | Probability of zero-inflation per gene   |
| `rng`          | `Random`| `None`  | Random number generator                  |

---

## Time-Series and Spatial Expression

### `simulate_time_series_expression`

Generate oscillatory gene expression over time points. Each gene has a random
frequency, phase, and amplitude. Final counts include Poisson noise.

```python
from metainformant.simulation import simulate_time_series_expression

counts = simulate_time_series_expression(n_timepoints=24, n_genes=100)
# Shape: (24, 100) -- 24 time points x 100 genes
```

### `simulate_spatial_expression`

Simulate spatially resolved expression data with configurable spatial patterns.

```python
from metainformant.simulation import simulate_spatial_expression

expression, coordinates = simulate_spatial_expression(
    n_spots=200, n_genes=500, spatial_patterns="gradient",
)
# expression.shape: (200, 500)
# coordinates.shape: (200, 2) -- (x, y) positions
```

Supported `spatial_patterns`: `"random"`, `"gradient"`, `"clusters"`.

---

## Technical Noise

### `add_technical_noise`

Layer PCR amplification bias and shot noise on top of an existing count matrix.

```python
from metainformant.simulation import simulate_rnaseq_counts, add_technical_noise

clean_counts = simulate_rnaseq_counts(n_genes=1000, n_samples=10)
noisy_counts = add_technical_noise(
    clean_counts,
    amplification_bias=0.1,
    sequencing_depth=1_000_000,
)
```

**Parameters:**

| Parameter            | Type         | Default | Description                            |
|----------------------|--------------|---------|----------------------------------------|
| `expression_matrix`  | `np.ndarray` | required| Input count matrix                     |
| `amplification_bias` | `float`      | `0.1`   | PCR bias strength (0-1)                |
| `sequencing_depth`   | `float|None` | `None`  | Target library size (rescales counts)  |
| `rng`                | `Random`     | `None`  | Random number generator                |

---

## See Also

- **[Sequence Simulation](sequences.md)** -- DNA/protein sequence generation
- **[Simulation Overview](index.md)** -- Full module architecture
- `metainformant.simulation.models.rna` -- Source module
