# Math Core

Statistical utility functions for mathematical biology.

## Contents

| File | Purpose |
|------|---------|
| `utilities.py` | General statistical functions: correlation, regression, entropy |

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

## Usage

```python
from metainformant.math.core.utilities import correlation, shannon_entropy

r = correlation([1.0, 2.0, 3.0], [1.1, 2.2, 2.9])
h = shannon_entropy([0.25, 0.25, 0.25, 0.25])
```

Population genetics and mathematical biology plotting now lives in
[`src/metainformant/visualization/analysis/math_plots.py`](../../../visualization/analysis/math_plots.py).
