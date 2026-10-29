# Population Genetics Statistical Testing

Statistical testing helpers for population genetics analyses: bootstrap methods, permutation tests, outlier detection, group comparisons, and confidence intervals. They live in `metainformant.math.population_genetics.statistics`.

## Overview

The helpers implement mostly non-parametric methods that are useful for robust population genetics inference, particularly with small sample sizes or non-normal distributions. Functions that return p-values are inferential; see the statistical-claim boundary in `src/metainformant/math/AGENTS.md`.

## Functions

### Bootstrap Confidence Intervals

**`bootstrap_confidence_interval(values, statistic_func=None, confidence_level=0.95, n_bootstrap=None, n_bootstraps=1000, random_state=None)`**

Percentile bootstrap interval for any statistic.

**Parameters:**
- `values`: Sequence of observed values
- `statistic_func`: Function that computes the statistic from a sample (default: `numpy.mean`)
- `confidence_level`: Confidence level (default: 0.95 for 95% CI)
- `n_bootstraps`: Number of bootstrap resamples (default: 1000); `n_bootstrap` is an alias that takes precedence when given
- `random_state`: Integer seed for a private `numpy` generator. The same seed gives the same interval; `None` draws fresh entropy, so the interval varies between calls. numpy's global random state is not used

**Returns:**
Dictionary with `statistic` (the statistic of the original values), `ci_lower`, `ci_upper`, `confidence_level` and `n_bootstrap`. Empty `values` give NaN bounds and `n_bootstrap` 0.

**Example:**
```python-snippet
from metainformant.math.population_genetics.statistics import bootstrap_confidence_interval
import numpy as np

data = [0.01, 0.02, 0.03, 0.04, 0.05]
result = bootstrap_confidence_interval(data, np.mean, n_bootstraps=1000, random_state=1)
print(f"Mean: {result['statistic']:.4f}")
print(f"95% CI: [{result['ci_lower']:.4f}, {result['ci_upper']:.4f}]")
```

### Permutation Tests

**`permutation_test(values1, values2, n_permutations=1000, random_state=None)`**

Permutation test for a difference in means between two groups.

**Parameters:**
- `values1`: First group of observations
- `values2`: Second group of observations
- `n_permutations`: Number of random relabelings of the pooled values (default: 1000)
- `random_state`: Integer seed for a private `numpy` generator (`None` draws fresh entropy)

**Returns:**
Dictionary with `test_statistic` (also under `statistic`), the observed absolute difference in means; `p_value`, calculated as (1 + k) / (1 + n_permutations), where k is the number of permutations whose absolute difference is at least as large; `n_permutations`; and `significant` (`p_value < 0.05`). An empty group raises `ValueError`.

**Example:**
```python-snippet
from metainformant.math.population_genetics.statistics import permutation_test

pop1_diversity = [0.01, 0.02, 0.03]
pop2_diversity = [0.04, 0.05, 0.06]
result = permutation_test(pop1_diversity, pop2_diversity, random_state=7)
print(f"P-value: {result['p_value']:.4f}")
```

### Outlier Detection

**`detect_outliers(values, method="iqr", threshold=1.5)`**

Detect outliers in a sequence of values.

**Parameters:**
- `values`: Sequence of values to test
- `method`: `"iqr"` (default; values outside Q1 − threshold·IQR and Q3 + threshold·IQR) or `"zscore"` (absolute z-score above `threshold`, using the population standard deviation); any other value raises `ValueError`
- `threshold`: Fence multiplier for `"iqr"`, or the z-score cutoff for `"zscore"` (the default 1.5 suits the IQR method; use a larger value such as 3.0 for z-scores)

**Returns:**
Dictionary with `outliers` (also under `outlier_values`), `outlier_indices`, `method`, `threshold`, `total_values` and `outlier_percentage`. An empty input returns only empty `outliers` and `outlier_indices`.

**Example:**
```python-snippet
from metainformant.math.population_genetics.statistics import detect_outliers

tajimas_d_values = [-0.5, -0.3, 0.1, 5.0, 0.2]
result = detect_outliers(tajimas_d_values, threshold=2.0)
print(f"Outliers at indices: {result['outlier_indices']}")
```

### Statistical Comparisons

**`compare_statistics(stats1, stats2, test_type="ttest")`**

Compare two groups of values, or two dictionaries of named statistics.

**Parameters:**
- `stats1`, `stats2`: Two sequences of values, or two dictionaries of statistic names to values
- `test_type`: one of `COMPARISON_TESTS`: `"ttest"` (Welch's t-test, default), `"mannwhitney"` (two-sided Mann-Whitney U) or `"kruskal"` (Kruskal-Wallis H); any other value raises `ValueError`. Input that cannot be tested raises `ValueError` instead of returning a number: empty, non-finite or non-1-D groups, fewer than two observations per group for the t-test, constant groups, or all-identical values for the rank tests. Dictionaries of summary statistics support only `"ttest"`, because there is nothing to rank.

**Returns:**
For sequences, a dictionary with `test_statistic`, `p_value` and `test_type`. For dictionaries, a dictionary keyed by each statistic name present in both inputs, holding `statistic_name`, `population1_value`, `population2_value`, `difference`, `test_statistic`, `p_value` and `significant`; those p-values come from a rough normal approximation whose standard error is built from the `variance` entries (default 1.0) as if each population had 10 observations, so treat them as indicative only. A zero standard error raises `ValueError` instead of producing a p-value.

### Confidence Intervals

**`calculate_confidence_intervals(values, confidence_level=0.95, method="bootstrap")`**

Interval for a list of values: a confidence interval for the mean (`"bootstrap"`, `"normal"`) or the central range of the values themselves (`"percentile"`).

**Parameters:**
- `values`: List or 1-D array of observed values, or a dictionary mapping statistic names to lists of observations; empty input raises `ValueError`
- `confidence_level`: Confidence level, strictly between 0 and 1 (default: 0.95)
- `method`: `"bootstrap"` (default), `"normal"` or `"percentile"`; any other value raises `ValueError`
  - `"bootstrap"`: percentile interval of the means of `min(1000, 10 × len(values))` resamples. The resampling uses a private `numpy.random.RandomState(42)`, so the interval for a given list is always the same, and numpy's global random state is neither read nor changed
  - `"normal"`: mean ± z × standard error (z = 1.96 at 0.95, 2.576 at 0.99, otherwise the normal quantile)
  - `"percentile"`: the empirical percentiles of the values themselves

**Returns:**
A tuple `(lower, upper)` for a list or array. For a dictionary, `{name: {"ci_lower": ..., "ci_upper": ...}}`, each interval computed from that statistic's observations by the chosen method. A dictionary entry holding a single point value raises `ValueError`, because one value carries no sampling information.

**Example:**
```python-snippet
from metainformant.math.population_genetics.statistics import calculate_confidence_intervals

diversity = [0.1, 0.15, 0.12, 0.18, 0.09]
lower, upper = calculate_confidence_intervals(diversity)  # bootstrap, fixed seed
lower99, upper99 = calculate_confidence_intervals(diversity, confidence_level=0.99, method="normal")
```

## Use Cases

1. **Bootstrap Confidence Intervals**: Estimate uncertainty in diversity estimates (π, θ) when sample sizes are small
2. **Permutation Tests**: Test for differences between populations without assuming normal distributions
3. **Outlier Detection**: Identify loci with unusual Tajima's D or other statistics that may indicate selection
4. **Statistical Comparisons**: Compare diversity, Fst, or other statistics between populations with appropriate tests

## References

- Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.
- Manly, B. F. (2006). Randomization, bootstrap and Monte Carlo methods in biology. CRC press.
