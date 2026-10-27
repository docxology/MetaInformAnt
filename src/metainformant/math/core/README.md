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
| `fisher_exact_test()` | Exact two-sided Fisher test on a 2x2 table `(a, b, c, d)`; returns `(odds_ratio, p_value)` from `scipy.stats.fisher_exact` and raises `ImportError` when SciPy is missing (there is no approximate fallback) |
| `covariance()` | Covariance between two variables |
| `shannon_entropy()` | Shannon entropy of a distribution |
| `jensen_shannon_divergence()` | Jensen-Shannon divergence between distributions |

## Usage

```python
from metainformant.math.core.utilities import correlation, fisher_exact_test, shannon_entropy

r = correlation([1.0, 2.0, 3.0], [1.1, 2.2, 2.9])
h = shannon_entropy([0.25, 0.25, 0.25, 0.25])
odds_ratio, p_value = fisher_exact_test(8, 2, 1, 9)  # table [[8, 2], [1, 9]]: odds ratio 36.0, p about 0.0055
```

`fisher_exact_test` reports the sample odds ratio `a*d / (b*c)` (`inf` when `b*c` is 0) and the exact two-sided conditional p-value; a table with an empty row or column gives `(nan, 1.0)`.

