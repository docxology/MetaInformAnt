### Math: Effective Population Size

Functions for estimating effective population size (Ne) from various parameters.

**Functions**: `harmonic_mean_effective_size`, `effective_size_sex_ratio`, `effective_size_from_family_size_variance`.

## Overview

Effective population size (Ne) is the size of an ideal population that would
experience the same amount of genetic drift as the actual population. Ne is
typically smaller than census population size due to factors like:
- Unequal sex ratios
- Variance in family size
- Population size fluctuations
- Overlapping generations

## Functions

### `harmonic_mean_effective_size(census_sizes)`

Calculate effective population size across multiple generations using harmonic mean.

**Parameters**:
- `census_sizes`: Iterable of census population sizes across generations

**Returns**: Harmonic mean effective size (float)

**Formula**: Ne = n / Σ(1/N_i) where n is number of generations

**Note**: The harmonic mean gives more weight to smaller population sizes,
reflecting their greater impact on genetic drift.

**Example**:
```python
from metainformant.math import harmonic_mean_effective_size

# Population sizes across 4 generations
sizes = [1000, 500, 2000, 800]
ne = harmonic_mean_effective_size(sizes)
# Returns: ~800 (harmonic mean)
```

### `effective_size_sex_ratio(num_males, num_females)`

Calculate effective population size with unequal sex ratio.

**Parameters**:
- `num_males`: Number of breeding males (Nm)
- `num_females`: Number of breeding females (Nf)

**Returns**: Effective population size (float)

**Formula**: Ne = 4 × Nm × Nf / (Nm + Nf)

**Example**:
```python
from metainformant.math import effective_size_sex_ratio

# Equal sex ratio
ne_equal = effective_size_sex_ratio(100, 100)
# Returns: 200.0

# Unequal sex ratio (more females)
ne_unequal = effective_size_sex_ratio(10, 90)
# Returns: 36.0 (much smaller than census size of 100)
```

### `effective_size_from_family_size_variance(census_size, variance_offspring_number)`

Calculate effective population size using Crow and Denniston approximation.

**Parameters**:
- `census_size`: Census population size (N)
- `variance_offspring_number`: Variance in family size (Vk)

**Returns**: Effective population size (float)

**Formula**: Ne = (4N - 2) / (Vk + 2)

**Example**:
```python
from metainformant.math import effective_size_from_family_size_variance

# Population with high variance in family size
ne = effective_size_from_family_size_variance(census_size=1000, variance_offspring_number=5)
# Higher variance reduces effective size
```

## Usage Examples

### Multiple Generations

```python
from metainformant.math import harmonic_mean_effective_size

# Track population size over time
generation_sizes = [1000, 800, 500, 1200, 900]
ne = harmonic_mean_effective_size(generation_sizes)
print(f"Effective size: {ne:.1f}")
```

### Sex Ratio Effects

```python
from metainformant.math import effective_size_sex_ratio

# Compare different sex ratios
ne1 = effective_size_sex_ratio(50, 50)   # Equal
ne2 = effective_size_sex_ratio(10, 90)   # Skewed
ne3 = effective_size_sex_ratio(90, 10)   # Reversed

print(f"Equal: {ne1}, Skewed: {ne2}, Reversed: {ne3}")
```

### Family Size Variance

```python
from metainformant.math import effective_size_from_family_size_variance

# Low variance (more equal reproduction)
ne_low = effective_size_from_family_size_variance(1000, variance_offspring_number=2)

# High variance (unequal reproduction)
ne_high = effective_size_from_family_size_variance(1000, variance_offspring_number=10)

print(f"Low variance Ne: {ne_low}, High variance Ne: {ne_high}")
```

## Integration with Other Modules

### With Population Genetics

```python
from metainformant.math import (
    effective_size_sex_ratio,
    effective_population_size_from_heterozygosity,
)

# Method 1: From sex ratio
ne_sex = effective_size_sex_ratio(num_males=100, num_females=100)

# Method 2: From heterozygosity (requires mutation rate)
ne_het = effective_population_size_from_heterozygosity(
    observed_heterozygosity=0.5,
    mutation_rate=1e-8,
    ploidy=2
)
```

## Notes

1. **Multiple Factors**: Effective size can be reduced by multiple factors
   simultaneously. The harmonic mean method accounts for population size
   fluctuations over time.

2. **Sex Ratio**: Unequal sex ratios significantly reduce Ne, even if
   census size is large.

3. **Family Size Variance**: Higher variance in reproductive success
   (some individuals have many offspring, others none) reduces Ne.

## References

- Crow, J. F., & Kimura, M. (1970). *An introduction to population genetics theory*.
  Harper & Row.
- Wright, S. (1931). Evolution in Mendelian populations. *Genetics*, 16(2), 97-159.

