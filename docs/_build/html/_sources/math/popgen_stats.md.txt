# Population Genetics Statistical Testing

This module provides statistical testing methods for population genetics analyses, including bootstrap methods, permutation tests, outlier detection, and confidence interval calculations.

## Overview

The `popgen_stats` module implements non-parametric statistical methods that are essential for robust population genetics inference, particularly when dealing with small sample sizes or non-normal distributions.

## Functions

### Bootstrap Confidence Intervals

**`bootstrap_confidence_interval(data, statistic_func, n_bootstrap=1000, confidence_level=0.95, method="percentile", random_state=None)`**

Calculate bootstrap confidence intervals for any statistic.

**Parameters:**
- `data`: Sequence of observed values
- `statistic_func`: Function that computes the statistic from data
- `n_bootstrap`: Number of bootstrap replicates (default: 1000)
- `confidence_level`: Confidence level (default: 0.95 for 95% CI)
- `method`: Method for CI calculation ("percentile" or "bca")
- `random_state`: Random seed for reproducibility

**Returns:**
Dictionary with statistic value, CI bounds, confidence level, and number of bootstrap replicates.

**Example:**
```python
from metainformant.math.popgen_stats import bootstrap_confidence_interval
import numpy as np

data = [0.01, 0.02, 0.03, 0.04, 0.05]
result = bootstrap_confidence_interval(data, np.mean, n_bootstrap=1000)
print(f"Mean: {result['statistic']:.4f}")
print(f"95% CI: [{result['ci_lower']:.4f}, {result['ci_upper']:.4f}]")
```

### Permutation Tests

**`permutation_test(group1, group2, statistic_func=None, n_permutations=10000, alternative="two-sided", random_state=None)`**

Perform permutation test to compare two groups.

**Parameters:**
- `group1`: First group of observations
- `group2`: Second group of observations
- `statistic_func`: Function to compute test statistic (default: difference of means)
- `n_permutations`: Number of permutations (default: 10000)
- `alternative`: Alternative hypothesis ("two-sided", "greater", "less")
- `random_state`: Random seed

**Returns:**
Dictionary with test statistic, p-value, number of permutations, and alternative hypothesis.

**Example:**
```python
from metainformant.math.popgen_stats import permutation_test

pop1_diversity = [0.01, 0.02, 0.03]
pop2_diversity = [0.04, 0.05, 0.06]
result = permutation_test(pop1_diversity, pop2_diversity)
print(f"P-value: {result['p_value']:.4f}")
```

### Outlier Detection

**`detect_outliers(values, method="zscore", threshold=3.0, fdr_correction=False)`**

Detect outliers in a sequence of values.

**Parameters:**
- `values`: Sequence of values to test
- `method`: Detection method ("zscore" or "iqr")
- `threshold`: Threshold for outlier detection (standard deviations for zscore)
- `fdr_correction`: Whether to apply FDR correction for multiple comparisons

**Returns:**
Dictionary with outlier indices, outlier values, z-scores, and p-values.

**Example:**
```python
from metainformant.math.popgen_stats import detect_outliers

tajimas_d_values = [-0.5, -0.3, 0.1, 5.0, 0.2]
result = detect_outliers(tajimas_d_values, threshold=2.0)
print(f"Outliers at indices: {result['outlier_indices']}")
```

### Statistical Comparisons

**`compare_statistics(stat1, stat2, test_type="mannwhitney")`**

Compare statistics between two groups using various statistical tests.

**Parameters:**
- `stat1`: Statistics from first group
- `stat2`: Statistics from second group
- `test_type`: Type of test ("mannwhitney", "ttest", "kruskal")

**Returns:**
Dictionary with test statistic, p-value, and test type.

### Confidence Intervals

**`calculate_confidence_intervals(statistics, data=None, method="normal", confidence_level=0.95)`**

Calculate confidence intervals for population genetics statistics.

**Parameters:**
- `statistics`: Dictionary of statistic names to values
- `data`: Optional raw data for bootstrap method
- `method`: Method for CI calculation ("normal" or "bootstrap")
- `confidence_level`: Confidence level

**Returns:**
Dictionary mapping statistic names to CI dictionaries.

## Use Cases

1. **Bootstrap Confidence Intervals**: Estimate uncertainty in diversity estimates (π, θ) when sample sizes are small
2. **Permutation Tests**: Test for differences between populations without assuming normal distributions
3. **Outlier Detection**: Identify loci with unusual Tajima's D or other statistics that may indicate selection
4. **Statistical Comparisons**: Compare diversity, Fst, or other statistics between populations with appropriate tests

## References

- Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.
- Manly, B. F. (2006). Randomization, bootstrap and Monte Carlo methods in biology. CRC press.


