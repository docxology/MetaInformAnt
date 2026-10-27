# Math: Demographic Models

Demographic models for population genetics analysis.

**Functions**: `exponential_growth_effective_size`, `bottleneck_effective_size`, `two_epoch_effective_size`.

## Overview

Demographic history (population size changes) affects genetic diversity and the distribution of polymorphisms. This module provides basic demographic models to estimate effective population size under different demographic scenarios.

## Functions

### `exponential_growth_effective_size(initial_size, growth_rate, generations)`

Calculate effective population size under exponential growth or decline.

**Parameters**:
- `initial_size`: Population size at generation 0 (N₀)
- `growth_rate`: Per-generation growth rate (r > −1). Positive for growth, negative for decline
- `generations`: Number of generations (t)

**Returns**: Effective population size, the harmonic mean of the census sizes over generations 0..t−1

**Formula**: The census size follows `exponential_growth_model`, Nᵢ = N₀ (1 + r)ⁱ, and Ne = t / Σᵢ 1/Nᵢ (Wright 1938).

**Example**:
```python
from metainformant.math.population_genetics.demography import exponential_growth_effective_size

# Population growing from 1000 by 23% per generation for 10 generations (to ~6400)
ne = exponential_growth_effective_size(initial_size=1000, growth_rate=0.23, generations=10)
# Returns: ~2140, weighted toward the early, small generations
```

### `bottleneck_effective_size(pre_bottleneck_size, bottleneck_size, post_bottleneck_size, recovery_generations=0)`

Calculate effective population size through a bottleneck.

**Parameters**:
- `pre_bottleneck_size`: Population size before bottleneck
- `bottleneck_size`: Population size during bottleneck
- `post_bottleneck_size`: Population size after bottleneck
- `recovery_generations`: Number of generations of recovery after bottleneck (default: 0)

**Returns**: Effective population size, the harmonic mean over one generation at each of the three sizes plus `recovery_generations` further generations at the post-bottleneck size: Ne = (3 + k) / (1/N_pre + 1/N_bottleneck + (1 + k)/N_post)

**Example**:
```python
from metainformant.math.population_genetics.demography import bottleneck_effective_size

# Severe bottleneck: 10000 → 100 → 10000
ne = bottleneck_effective_size(
    pre_bottleneck_size=10000,
    bottleneck_size=100,
    post_bottleneck_size=10000
)
# Returns: ~294 (the bottleneck generation dominates)

# With 10 further generations at the recovered size
ne_recovery = bottleneck_effective_size(
    pre_bottleneck_size=10000,
    bottleneck_size=100,
    post_bottleneck_size=10000,
    recovery_generations=10
)
# Returns: ~1161
```

`math.population_genetics.statistics.bottleneck_effective_size(initial_size, bottleneck_size, bottleneck_duration, final_size)` is the same harmonic mean for a bottleneck lasting `bottleneck_duration` generations between one generation at each of the other sizes.

### `two_epoch_effective_size(ancient_size, recent_size, split_time)`

Calculate the coalescent effective size of a population whose size changed once.

**Parameters**:
- `ancient_size`: Population size before the change
- `recent_size`: Population size since the change
- `split_time`: Generations since the size changed (T)

**Returns**: The constant population size with the same expected pairwise coalescence time, and hence the same expected nucleotide diversity. For a diploid population Ne = N_recent + (N_ancient − N_recent) e^(−T / 2N_recent).

**Example**:
```python
from metainformant.math.population_genetics.demography import two_epoch_effective_size

# Population expansion: 1000 → 10000, 50 generations ago
ne = two_epoch_effective_size(ancient_size=1000, recent_size=10000, split_time=50)
# Returns: ~1022 (diversity has barely begun to rise)

# The reverse contraction is not symmetric
ne = two_epoch_effective_size(ancient_size=10000, recent_size=1000, split_time=50)
# Returns: ~9778 (diversity has barely begun to fall)
```

## Use Cases

### 1. Interpreting Genetic Diversity

Effective population size affects genetic diversity:

```python-snippet
from metainformant.math.population_genetics.demography import exponential_growth_effective_size
from metainformant.dna.population import nucleotide_diversity

# Calculate diversity
sequences = ["AAAA", "AAAT", "AATT"]
pi = nucleotide_diversity(sequences)

# Estimate effective size from growth model
ne = exponential_growth_effective_size(
    initial_size=10000,
    growth_rate=0.1,
    generations=100
)

# Compare observed diversity to expected under model
# (requires mutation rate)
```

### 2. Modeling Population History

Combine multiple demographic events:

```python-snippet
from metainformant.math.population_genetics.demography import (
    bottleneck_effective_size,
    exponential_growth_effective_size,
)

# Model: Bottleneck, then recovery and growth
# Step 1: Bottleneck
ne_bottleneck = bottleneck_effective_size(
    pre_bottleneck_size=10000,
    bottleneck_size=100,
    post_bottleneck_size=10000,
    recovery_generations=10
)

# Step 2: Growth after recovery
ne_growth = exponential_growth_effective_size(
    initial_size=10000,
    growth_rate=0.05,
    generations=50
)

# Overall effective size is dominated by bottleneck
```

### 3. Comparing Demographic Models

Compare different demographic scenarios:

```python-snippet
from metainformant.math.population_genetics.demography import (
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

Over t generations of changing census size, the effective population size is the harmonic mean (Wright 1938):

Ne = t / (1/N₁ + 1/N₂ + ... + 1/Nₜ)

Drift is strongest when the population is small, so small generations dominate the effective size.

### Exponential Growth

Under discrete exponential growth Nᵢ = N₀ (1 + r)ⁱ, the effective size lies between the initial and current sizes and below their arithmetic mean, because the harmonic mean is weighted by the smaller early sizes.

### Two Epochs

Going back in time, two lineages coalesce at rate 1/(2N) per generation. With rate 1/(2N_recent) for the T generations since a size change and 1/(2N_ancient) before it, the expected coalescence time is 2[N_recent + (N_ancient − N_recent) e^(−T / 2N_recent)], twice the effective size `two_epoch_effective_size` reports.

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

```python-snippet
from metainformant.math.population_genetics.demography import exponential_growth_effective_size
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

```python-snippet
from metainformant.math.population_genetics.demography import bottleneck_effective_size
from metainformant.math.population_genetics.effective_size import effective_population_size_from_heterozygosity

# Estimate Ne from observed heterozygosity
observed_H = 0.3
mutation_rate = 0.0001
ne_estimated = effective_population_size_from_heterozygosity(
    heterozygosity=observed_H,
    mutation_rate=mutation_rate
)

# Compare to demographic model
ne_model = bottleneck_effective_size(10000, 100, 10000)
```

## References

- **Nei, M., Maruyama, T., & Chakraborty, R. (1975)**. The bottleneck effect and genetic variability in populations. *Evolution*, 29(1), 1-10.

- **Maruyama, T., & Fuerst, P. A. (1984)**. Population bottlenecks and nonequilibrium models in population genetics. I. Allele numbers when populations evolve from zero variability. *Genetics*, 108(3), 745-763.

- **Slatkin, M., & Hudson, R. R. (1991)**. Pairwise comparisons of mitochondrial DNA sequences in stable and exponentially growing populations. *Genetics*, 129(2), 555-562.

## See Also

- [`math.effective_size`](effective_size.md) - Effective population size estimators
- [`math.population_genetics.coalescent`](coalescent.md) - Coalescent theory and expected diversity
- [`math.popgen`](../math/popgen.md) - Population genetics models

