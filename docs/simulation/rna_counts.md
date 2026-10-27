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
from metainformant.simulation.models.rna import simulate_counts_negative_binomial

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
from metainformant.simulation.models.rna import simulate_rnaseq_counts

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
of genes. Samples are split between the two conditions (an odd total gives group 1
the extra sample).

Each gene starts from a baseline mean drawn from an exponential distribution with
mean 10; its dispersion is `1 / (mean + 1)` clipped to the range 0.1 to 10, and both
conditions use it. Group 0 keeps the baseline mean. A fold change rescales the
**group-1 mean** of the corresponding gene, and it is a *signed ratio*:

| `fold_changes` value | Group-1 mean                                                                 |
|----------------------|------------------------------------------------------------------------------|
| `fc > 0`             | multiplied by `fc` (`2.0` doubles it, `0.5` halves it)                       |
| `fc < 0`             | divided by `abs(fc)`: a fold *decrease* (`-4.0` gives about 0.25x, `-2.0` halves it) |
| `fc == 0`            | unchanged                                                                    |

A negative value needs `abs(fc) >= 1` because a decrease cannot be smaller than one
fold; `-1.0` leaves the mean unchanged. To express a decrease as a ratio, pass a
positive value below one (`0.25`) instead.

```python
import numpy as np
from metainformant.simulation.models.rna import simulate_differential_expression

fold_changes = np.array([2.0, -4.0, 3.0, 0.5])
expression, labels = simulate_differential_expression(
    n_samples=20,
    n_features=5000,
    fold_changes=fold_changes,
)
# expression.shape: (20, 5000)
# labels: array of 0s and 1s indicating condition
# The first 4 genes are differentially expressed in group 1:
# multiplied by 2, divided by 4 (about 0.25x), multiplied by 3 and halved
```

**Parameters:**

| Parameter      | Type         | Description                                          |
|----------------|--------------|------------------------------------------------------|
| `n_samples`    | `int`        | Total samples across both conditions (>= 2)          |
| `n_features`   | `int`        | Total genes                                          |
| `fold_changes` | `np.ndarray` | Signed fold changes for the first `len(fold_changes)` genes (see the table above); must be an array of finite values with fewer entries than `n_features` |
| `rng`          | `Random`     | Optional random number generator                     |

**Returns:** Tuple of `(expression_matrix, group_labels)`

**Raises:** `ValidationError` (from `metainformant.core.utils.errors`; it does not derive
from `ValueError`) when

- `n_samples` is below 2, `n_features` is below 1, or `fold_changes` is not an `np.ndarray`
- `len(fold_changes)` is `n_features` or more
- a fold change is NaN or infinite
- a fold change lies strictly between -1 and 0

---

## Bulk RNA-seq

### `simulate_bulk_rnaseq`

Simulate bulk RNA-seq with realistic library sizes and multinomial count generation.

```python
from metainformant.simulation.models.rna import simulate_bulk_rnaseq

counts = simulate_bulk_rnaseq(n_samples=12, n_genes=2000)
# Library sizes follow a log-normal distribution (median e^15, about 3.3 million reads)
# Gene means follow a power-function distribution (many low, few high)
```

**Parameters:**

| Parameter       | Type              | Default | Description                          |
|-----------------|-------------------|---------|--------------------------------------|
| `n_samples`     | `int`             | required| Number of samples                    |
| `n_genes`       | `int`             | required| Number of genes                      |
| `library_sizes` | `np.ndarray \| None` | `None`  | Custom library sizes per sample   |
| `gene_means`    | `np.ndarray \| None` | `None`  | Custom mean expression per gene   |
| `rng`           | `Random`          | `None`  | Random number generator              |

---

## Single-Cell RNA-seq

### `simulate_single_cell_rnaseq`

Simulate scRNA-seq data with cell type-specific expression patterns and dropout events.
Each cell type has marker genes with elevated expression. Dropout is modeled as a
Bernoulli process applied after count generation.

```python
from metainformant.simulation.models.rna import simulate_single_cell_rnaseq

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

Generate sine-based gene expression over time points. Each gene follows one smooth
profile over the sampling times `t = linspace(0, 4π, n_timepoints)`:
`mean = max(0, baseline + amplitude · sin(frequency · t + phase))`. The baseline
(uniform 5 to 20), amplitude (uniform 10 to 100) and phase (uniform 0 to 2π) are drawn
once per gene and stay fixed over time; the frequency is uniform 0.5 to 2.0 unless
`oscillation_freq` supplies one per gene. Each observation is an independent Poisson
draw around the gene's curve, so successive time points differ only by counting noise.

```python
from metainformant.simulation.models.rna import simulate_time_series_expression

counts = simulate_time_series_expression(n_timepoints=24, n_genes=100)
# Shape: (24, 100) -- 24 time points x 100 genes
```

### `simulate_spatial_expression`

Simulate spatially resolved expression data with configurable spatial patterns.

```python
from metainformant.simulation.models.rna import simulate_spatial_expression

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
from metainformant.simulation.models.rna import simulate_rnaseq_counts, add_technical_noise

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
| `sequencing_depth`   | `float \| None` | `None`  | Target library size (rescales counts)  |
| `rng`                | `Random`     | `None`  | Random number generator                |

---

## See Also

- **[Sequence Simulation](sequences.md)** -- DNA/protein sequence generation
- **[Simulation Overview](index.md)** -- Full module architecture
- `metainformant.simulation.models.rna` -- Source module
