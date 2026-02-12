# Math Core

Statistical utility functions and population genetics visualization for mathematical biology.

## Contents

| File | Purpose |
|------|---------|
| `utilities.py` | General statistical functions: correlation, regression, entropy |
| `visualization.py` | Plotting for allele frequencies, drift, coalescent trees, epidemics |

## Key Functions

| Function | Description |
|----------|-------------|
| `correlation()` | Pearson correlation coefficient between two variables |
| `linear_regression()` | Simple linear regression with slope, intercept, R-squared |
| `r_squared()` | Coefficient of determination |
| `fisher_exact_test()` | Fisher exact test on 2x2 contingency table |
| `covariance()` | Covariance between two variables |
| `shannon_entropy()` | Shannon entropy of a distribution |
| `jensen_shannon_divergence()` | Jensen-Shannon divergence between distributions |
| `plot_allele_frequency_spectrum()` | Allele frequency spectrum histogram |
| `plot_coalescent_tree()` | Coalescent tree visualization |
| `plot_genetic_drift_simulation()` | Genetic drift trajectory plots |
| `plot_epidemic_model_simulation()` | SIR/SEIR epidemic model curves |
| `plot_linkage_disequilibrium_decay()` | LD decay over physical distance |

## Usage

```python
from metainformant.math.core.utilities import correlation, shannon_entropy
from metainformant.math.core.visualization import plot_allele_frequency_spectrum

r = correlation([1.0, 2.0, 3.0], [1.1, 2.2, 2.9])
h = shannon_entropy([0.25, 0.25, 0.25, 0.25])
plot_allele_frequency_spectrum(freqs, output_path="output/afs.png")
```
