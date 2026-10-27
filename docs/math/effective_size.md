### Math: Effective Population Size

Functions for estimating effective population size (Ne) from various parameters.

**Functions**: `harmonic_mean_effective_size`, `effective_size_sex_ratio`, `effective_population_size_from_heterozygosity` (all in `metainformant.math.population_genetics.effective_size`) and `effective_size_from_family_size_variance` (in `metainformant.math.population_genetics.statistics`).

## Overview

Effective population size (Ne) is the size of an ideal population that would
experience the same amount of genetic drift as the actual population. Ne is
typically smaller than census population size due to factors like:
- Unequal sex ratios
- Variance in family size
- Population size fluctuations
- Overlapping generations

## Functions

### `harmonic_mean_effective_size(population_sizes)`

Calculate effective population size across multiple generations using harmonic mean.

**Parameters**:
- `population_sizes`: List of census population sizes across generations (all must be positive; an empty list raises `ValueError`)

**Returns**: Harmonic mean effective size (float)

**Formula**: Ne = n / Σ(1/N_i) where n is number of generations

**Note**: The harmonic mean gives more weight to smaller population sizes,
reflecting their greater impact on genetic drift.

**Example**:
```python
from metainformant.math.population_genetics.effective_size import harmonic_mean_effective_size

# Population sizes across 4 generations
sizes = [1000, 500, 2000, 800]
ne = harmonic_mean_effective_size(sizes)
# Returns: ~842.1 (harmonic mean)
```

### `effective_size_sex_ratio(male_size, female_size, sex_ratio=0.5)`

Calculate effective population size with unequal sex ratio.

**Parameters**:
- `male_size`: Number of breeding males (Nm)
- `female_size`: Number of breeding females (Nf)
- `sex_ratio`: Proportion of the population that is male (default 0.5; must be strictly between 0 and 1)

**Returns**: Effective population size (float)

**Formula**: Ne = 4 × Nm × Nf / (Nm + Nf) at the default `sex_ratio=0.5`. Any other `sex_ratio` r switches to a simplified adjustment, 2·Nm·Nf / (Nm + Nf) / (r(1 − r)), which is not continuous with the default formula, so leave `sex_ratio` at its default unless that approximation is what you want.

**Example**:
```python
from metainformant.math.population_genetics.effective_size import effective_size_sex_ratio

# Equal sex ratio
ne_equal = effective_size_sex_ratio(100, 100)
# Returns: 200.0

# Unequal sex ratio (more females)
ne_unequal = effective_size_sex_ratio(10, 90)
# Returns: 36.0 (much smaller than census size of 100)
```

### `effective_population_size_from_heterozygosity(heterozygosity, mutation_rate=1e-8)`

Estimate effective population size from observed heterozygosity, using the diploid relation H = 4·Ne·μ / (1 + 4·Ne·μ).

**Parameters**:
- `heterozygosity`: Observed heterozygosity H (strictly between 0 and 1)
- `mutation_rate`: Mutation rate per site per generation (default 1e-8; must be positive)

**Returns**: Effective population size (float)

**Formula**: Ne = H / (4μ(1 − H))

**Example**:
```python
from metainformant.math.population_genetics.effective_size import effective_population_size_from_heterozygosity

ne = effective_population_size_from_heterozygosity(0.5, mutation_rate=1e-8)
# Returns: 25000000.0
```

### `effective_size_from_family_size_variance(family_sizes, Vk=None)`

Calculate effective population size from variance in family size. The function has two forms:

- With `Vk`, `family_sizes` is the census size N and must be a **float** (an `int` fails an internal assertion). It uses the Crow and Denniston approximation.
- Without `Vk`, `family_sizes` is a list of offspring counts, one per parent.

**Parameters**:
- `family_sizes`: Census population size as a float (when `Vk` is given), or a list of offspring counts per parent
- `Vk`: Variance in family size (Vk)

**Returns**: Effective population size (float)

**Formula**: With `Vk`, Ne = (4N - 2) / (Vk + 2). With a list, Ne = k / (V - 1), where k is the mean family size and V the population variance (denominator n) of the list; the result is `inf` when V is at most 1 and 0.0 for an empty list.

**Example**:
```python
from metainformant.math.population_genetics.statistics import effective_size_from_family_size_variance

# Population with high variance in family size
ne = effective_size_from_family_size_variance(1000.0, Vk=5)
# Returns: ~571.1; higher variance reduces effective size
```

## Usage Examples

### Multiple Generations

```python
from metainformant.math.population_genetics.effective_size import harmonic_mean_effective_size

# Track population size over time
generation_sizes = [1000, 800, 500, 1200, 900]
ne = harmonic_mean_effective_size(generation_sizes)
print(f"Effective size: {ne:.1f}")
```

### Sex Ratio Effects

```python
from metainformant.math.population_genetics.effective_size import effective_size_sex_ratio

# Compare different sex ratios
ne1 = effective_size_sex_ratio(50, 50)   # Equal
ne2 = effective_size_sex_ratio(10, 90)   # Skewed
ne3 = effective_size_sex_ratio(90, 10)   # Reversed

print(f"Equal: {ne1}, Skewed: {ne2}, Reversed: {ne3}")
```

### Family Size Variance

```python
from metainformant.math.population_genetics.statistics import effective_size_from_family_size_variance

# Low variance (more equal reproduction)
ne_low = effective_size_from_family_size_variance(1000.0, Vk=2)

# High variance (unequal reproduction)
ne_high = effective_size_from_family_size_variance(1000.0, Vk=10)

print(f"Low variance Ne: {ne_low}, High variance Ne: {ne_high}")
```

## Integration with Other Modules

### With Population Genetics

```python
from metainformant.math.population_genetics.effective_size import (
    effective_size_sex_ratio,
    effective_population_size_from_heterozygosity,
)

# Method 1: From sex ratio
ne_sex = effective_size_sex_ratio(male_size=100, female_size=100)

# Method 2: From heterozygosity (requires mutation rate)
ne_het = effective_population_size_from_heterozygosity(
    heterozygosity=0.5,
    mutation_rate=1e-8
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
