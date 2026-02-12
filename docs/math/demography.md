# Math: Demographic Models

Demographic models for population genetics analysis.

**Functions**: `exponential_growth_effective_size`, `bottleneck_effective_size`, `two_epoch_effective_size`.

## Overview

Demographic history (population size changes) affects genetic diversity and the distribution of polymorphisms. This module provides basic demographic models to estimate effective population size under different demographic scenarios.

## Functions

### `exponential_growth_effective_size(current_size, growth_rate, generations)`

Calculate effective population size under exponential growth.

**Parameters**:
- `current_size`: Current population size (N_t)
- `growth_rate`: Per-generation growth rate (r). Positive for growth, negative for decline
- `generations`: Number of generations over which to calculate harmonic mean

**Returns**: Effective population size (harmonic mean)

**Formula**: Under exponential growth N(t) = N₀ × e^(r×t), the effective size is the harmonic mean: Ne = t / Σ(1/N_i)

**Example**:
```python
from metainformant.math.demography import exponential_growth_effective_size

# Population growing from ~1000 to 10000 over 10 generations
ne = exponential_growth_effective_size(
    current_size=10000,
    growth_rate=0.23,  # ~23% per generation
    generations=10
)
# Returns: ~4342 (harmonic mean of sizes over time)
```

### `bottleneck_effective_size(pre_bottleneck_size, bottleneck_size, bottleneck_duration, recovery_generations=0)`

Calculate effective population size through a bottleneck.

**Parameters**:
- `pre_bottleneck_size`: Population size before bottleneck
- `bottleneck_size`: Population size during bottleneck
- `bottleneck_duration`: Number of generations at bottleneck size
- `recovery_generations`: Number of generations after bottleneck (default: 0)

**Returns**: Effective population size (harmonic mean)

**Example**:
```python
from metainformant.math.demography import bottleneck_effective_size

# Severe bottleneck: 10000 → 100 for 5 generations
ne = bottleneck_effective_size(
    pre_bottleneck_size=10000,
    bottleneck_size=100,
    bottleneck_duration=5
)
# Returns: ~109 (heavily weighted by bottleneck)

# With recovery to original size over 10 generations
ne_recovery = bottleneck_effective_size(
    pre_bottleneck_size=10000,
    bottleneck_size=100,
    bottleneck_duration=5,
    recovery_generations=10
)
# Returns: ~200 (recovery improves effective size)
```

### `two_epoch_effective_size(ancient_size, current_size, time_since_change)`

Calculate effective population size for two-epoch model.

**Parameters**:
- `ancient_size`: Population size before change
- `current_size`: Population size after change
- `time_since_change`: Number of generations since size change

**Returns**: Effective population size (harmonic mean)

**Example**:
```python
from metainformant.math.demography import two_epoch_effective_size

# Population expansion: 1000 → 10000, 50 generations ago
ne = two_epoch_effective_size(
    ancient_size=1000,
    current_size=10000,
    time_since_change=50
)
# Returns: ~1960 (harmonic mean of 1000 and 10000)
```

## Use Cases

### 1. Interpreting Genetic Diversity

Effective population size affects genetic diversity:

```python
from metainformant.math.demography import exponential_growth_effective_size
from metainformant.dna.population import nucleotide_diversity

# Calculate diversity
sequences = ["AAAA", "AAAT", "AATT"]
pi = nucleotide_diversity(sequences)

# Estimate effective size from growth model
ne = exponential_growth_effective_size(
    current_size=10000,
    growth_rate=0.1,
    generations=100
)

# Compare observed diversity to expected under model
# (requires mutation rate)
```

### 2. Modeling Population History

Combine multiple demographic events:

```python
from metainformant.math.demography import (
    bottleneck_effective_size,
    exponential_growth_effective_size,
)

# Model: Bottleneck, then recovery and growth
# Step 1: Bottleneck
ne_bottleneck = bottleneck_effective_size(
    pre_bottleneck_size=10000,
    bottleneck_size=100,
    bottleneck_duration=5,
    recovery_generations=10
)

# Step 2: Growth after recovery
ne_growth = exponential_growth_effective_size(
    current_size=10000,
    growth_rate=0.05,
    generations=50
)

# Overall effective size is dominated by bottleneck
```

### 3. Comparing Demographic Models

Compare different demographic scenarios:

```python
from metainformant.math.demography import (
    exponential_growth_effective_size,
    two_epoch_effective_size,
)

# Model 1: Exponential growth
ne_exp = exponential_growth_effective_size(10000, 0.1, 50)

# Model 2: Two-epoch (sudden change)
ne_epoch = two_epoch_effective_size(1000, 10000, 50)

# Compare effective sizes
print(f"Exponential growth Ne: {ne_exp:.2f}")
print(f"Two-epoch Ne: {ne_epoch:.2f}")
```

## Theory

### Harmonic Mean

Effective population size under changing population size is the harmonic mean:

Ne = t / (1/N₁ + 1/N₂ + ... + 1/Nₜ)

This reflects that coalescence rates are inversely proportional to population size, so periods of small population size dominate the effective size.

### Exponential Growth

Under exponential growth: N(t) = N₀ × e^(r×t)

The effective size is less than the current size because the harmonic mean is weighted by the smaller historical sizes.

### Bottlenecks

Bottlenecks dramatically reduce effective population size because:
1. Small population sizes during bottleneck dominate harmonic mean
2. Genetic drift is much stronger during bottleneck
3. Recovery period may not fully compensate

## Limitations

1. **Deterministic models**: These models assume deterministic population size changes. Stochastic fluctuations are not modeled.

2. **Simple scenarios**: Complex demographic histories (multiple bottlenecks, migration, etc.) require more sophisticated models.

3. **No migration**: Models assume isolated populations. Migration would require additional parameters.

4. **No selection**: Models assume neutral evolution. Selection would affect genetic diversity differently.

## Integration with Other Modules

### With Coalescent Theory

```python
from metainformant.math.demography import exponential_growth_effective_size
from metainformant.math.population_genetics.coalescent import expected_segregating_sites

# Calculate expected diversity under demographic model
ne = exponential_growth_effective_size(10000, 0.1, 100)
theta = 4 * ne * 0.0001  # Assuming mutation rate μ = 0.0001
S_expected = expected_segregating_sites(
    sample_size=10,
    theta=theta,
    sequence_length=1000
)
```

### With Effective Size Estimators

```python
from metainformant.math.demography import bottleneck_effective_size
from metainformant.math.popgen import effective_population_size_from_heterozygosity

# Estimate Ne from observed heterozygosity
observed_H = 0.3
mutation_rate = 0.0001
ne_estimated = effective_population_size_from_heterozygosity(
    observed_heterozygosity=observed_H,
    mutation_rate=mutation_rate,
    ploidy=2
)

# Compare to demographic model
ne_model = bottleneck_effective_size(10000, 100, 5)
```

## References

- **Nei, M., Maruyama, T., & Chakraborty, R. (1975)**. The bottleneck effect and genetic variability in populations. *Evolution*, 29(1), 1-10.

- **Maruyama, T., & Fuerst, P. A. (1984)**. Population bottlenecks and nonequilibrium models in population genetics. I. Allele numbers when populations evolve from zero variability. *Genetics*, 108(3), 745-763.

- **Slatkin, M., & Hudson, R. R. (1991)**. Pairwise comparisons of mitochondrial DNA sequences in stable and exponentially growing populations. *Genetics*, 129(2), 555-562.

## See Also

- [`math.effective_size`](effective_size.md) - Effective population size estimators
- [`math.population_genetics.coalescent`](coalescent.md) - Coalescent theory and expected diversity
- [`math.popgen`](../math/popgen.md) - Population genetics models

