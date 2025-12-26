# Simulation: Population Genetics

Synthetic data generation for population genetics analysis.

**Functions**: `generate_population_sequences`, `generate_two_populations`, `generate_genotype_matrix`, `simulate_bottleneck_population`, `simulate_population_expansion`, `generate_site_frequency_spectrum`, `generate_linkage_disequilibrium_data`.

## Overview

This module provides comprehensive synthetic data generation for population genetics analysis. You can generate:

- **Sequences with specified diversity** (π, θ_W)
- **Two populations with specified Fst**
- **Genotype matrices** with Hardy-Weinberg equilibrium
- **Demographic scenarios** (bottlenecks, expansions)
- **Site frequency spectra** under neutral models
- **Linkage disequilibrium** data

## Functions

### `generate_population_sequences(n_sequences, sequence_length, ...)`

Generate a population of sequences with specified diversity.

**Parameters**:
- `n_sequences`: Number of sequences to generate
- `sequence_length`: Length of each sequence
- `nucleotide_diversity`: Target π (average pairwise differences per site)
- `wattersons_theta`: Target θ_W (Watterson's estimator)
- `reference_sequence`: Starting sequence (optional)
- `mutation_rate`: Per-site mutation rate
- `gc_content`: GC content for reference sequence

**Returns**: List of sequences

**Example**:
```python
from metainformant.simulation.popgen import generate_population_sequences

# Generate 10 sequences with π = 0.01
sequences = generate_population_sequences(
    n_sequences=10,
    sequence_length=1000,
    nucleotide_diversity=0.01
)

# Analyze the generated data
from metainformant.dna.population import nucleotide_diversity
pi = nucleotide_diversity(sequences)
print(f"Observed π: {pi:.4f}")
```

### `generate_two_populations(n_pop1, n_pop2, sequence_length, fst=0.1, ...)`

Generate two populations with specified Fst.

**Parameters**:
- `n_pop1`: Number of sequences in population 1
- `n_pop2`: Number of sequences in population 2
- `sequence_length`: Length of each sequence
- `fst`: Target Fst value (0-1, higher = more differentiation)
- `within_pop_diversity`: Target nucleotide diversity within each population

**Returns**: Tuple of (pop1_sequences, pop2_sequences)

**Example**:
```python
from metainformant.simulation.popgen import generate_two_populations
from metainformant.dna.population import hudson_fst

# Generate two populations with Fst = 0.2
pop1, pop2 = generate_two_populations(
    n_pop1=10,
    n_pop2=10,
    sequence_length=1000,
    fst=0.2,
    within_pop_diversity=0.01
)

# Verify Fst
fst_observed = hudson_fst(pop1, pop2)
print(f"Observed Fst: {fst_observed:.4f}")
```

### `generate_genotype_matrix(n_individuals, n_sites, ...)`

Generate a genotype matrix with specified allele frequencies.

**Parameters**:
- `n_individuals`: Number of individuals (rows)
- `n_sites`: Number of sites (columns)
- `allele_frequencies`: Optional list of allele frequencies per site
- `min_maf`: Minimum minor allele frequency
- `max_maf`: Maximum minor allele frequency
- `hwe`: If True, genotypes follow Hardy-Weinberg equilibrium
- `ploidy`: Ploidy level (2 for diploid, 1 for haploid)

**Returns**: Genotype matrix (list of lists). For diploid: 0=AA, 1=AB, 2=BB.

**Example**:
```python
from metainformant.simulation.popgen import generate_genotype_matrix
from metainformant.dna.population import allele_frequencies

# Generate genotypes with specified frequencies
genotypes = generate_genotype_matrix(
    n_individuals=100,
    n_sites=10,
    allele_frequencies=[0.2, 0.3, 0.4, 0.1, 0.5, 0.3, 0.2, 0.4, 0.1, 0.6],
    hwe=True
)

# Verify frequencies
freqs = allele_frequencies(genotypes)
print(f"Allele frequencies: {freqs}")
```

### `simulate_bottleneck_population(n_sequences, sequence_length, ...)`

Simulate a population that went through a bottleneck.

**Parameters**:
- `n_sequences`: Number of sequences to generate
- `sequence_length`: Length of each sequence
- `pre_bottleneck_diversity`: Nucleotide diversity before bottleneck
- `bottleneck_size`: Effective population size during bottleneck
- `bottleneck_duration`: Number of generations at bottleneck size
- `recovery_generations`: Number of generations since bottleneck

**Returns**: List of sequences reflecting bottleneck signature

**Example**:
```python
from metainformant.simulation.popgen import simulate_bottleneck_population
from metainformant.dna.population import nucleotide_diversity, tajimas_d

# Simulate bottleneck: 10000 → 5 → recovery
sequences = simulate_bottleneck_population(
    n_sequences=20,
    sequence_length=1000,
    pre_bottleneck_diversity=0.01,
    bottleneck_size=5,
    bottleneck_duration=10,
    recovery_generations=20
)

# Bottleneck creates negative Tajima's D
pi = nucleotide_diversity(sequences)
d = tajimas_d(sequences)
print(f"π: {pi:.4f}, Tajima's D: {d:.4f}")  # D should be negative
```

### `simulate_population_expansion(n_sequences, sequence_length, ...)`

Simulate a population that underwent recent expansion.

**Parameters**:
- `n_sequences`: Number of sequences to generate
- `sequence_length`: Length of each sequence
- `initial_diversity`: Nucleotide diversity in small ancestral population
- `expansion_factor`: How much population grew (e.g., 10x)
- `growth_rate`: Per-generation growth rate

**Returns**: List of sequences reflecting expansion signature

**Example**:
```python
from metainformant.simulation.popgen import simulate_population_expansion

# Simulate 10x expansion
sequences = simulate_population_expansion(
    n_sequences=20,
    sequence_length=1000,
    expansion_factor=10.0,
    growth_rate=0.1
)

# Expansion also creates negative Tajima's D
from metainformant.dna.population import tajimas_d
d = tajimas_d(sequences)
print(f"Tajima's D: {d:.4f}")  # Should be negative
```

### `generate_site_frequency_spectrum(sample_size, n_sites, ...)`

Generate a site frequency spectrum with specified properties.

**Parameters**:
- `sample_size`: Number of sampled sequences (n)
- `n_sites`: Number of polymorphic sites to generate
- `theta`: Population mutation parameter θ
- `folded`: If True, return folded SFS (minor allele frequencies only)

**Returns**: List of counts per frequency bin

**Example**:
```python
from metainformant.simulation.popgen import generate_site_frequency_spectrum

# Generate folded SFS for 10 samples, 100 sites
sfs = generate_site_frequency_spectrum(
    sample_size=10,
    n_sites=100,
    theta=0.01,
    folded=True
)

# SFS should have more rare alleles (first bins)
print(f"SFS: {sfs}")  # Should show excess of rare alleles
```

### `generate_linkage_disequilibrium_data(n_individuals, n_sites, ...)`

Generate genotype data with specified linkage disequilibrium.

**Parameters**:
- `n_individuals`: Number of individuals
- `n_sites`: Number of sites
- `r_squared_target`: Target r² value (LD measure)
- `recombination_rate`: Recombination rate between sites
- `allele_frequencies`: Optional allele frequencies per site

**Returns**: Genotype matrix with specified LD patterns

**Example**:
```python
from metainformant.simulation.popgen import generate_linkage_disequilibrium_data

# Generate genotypes with LD
genotypes = generate_linkage_disequilibrium_data(
    n_individuals=100,
    n_sites=10,
    r_squared_target=0.5,
    recombination_rate=0.01
)

# Analyze LD
from metainformant.math.ld import r_squared
# Calculate LD between adjacent sites
ld = r_squared([row[0] for row in genotypes], [row[1] for row in genotypes])
print(f"Observed r²: {ld:.4f}")
```

## Use Cases

### 1. Testing Population Genetics Methods

Generate synthetic data to test analysis methods:

```python
from metainformant.simulation.popgen import generate_population_sequences
from metainformant.dna.population_analysis import calculate_summary_statistics

# Generate data with known properties
sequences = generate_population_sequences(
    n_sequences=20,
    sequence_length=1000,
    nucleotide_diversity=0.01
)

# Test analysis methods
stats = calculate_summary_statistics(sequences=sequences)
print(f"π: {stats['nucleotide_diversity']:.4f}")
print(f"S: {stats['segregating_sites']}")
print(f"θ_W: {stats['wattersons_theta']:.4f}")
```

### 2. Comparing Demographic Models

Compare different demographic scenarios:

```python
from metainformant.simulation.popgen import (
    simulate_bottleneck_population,
    simulate_population_expansion,
)
from metainformant.dna.population_analysis import neutrality_test_suite

# Bottleneck scenario
bottleneck_seqs = simulate_bottleneck_population(
    n_sequences=20,
    sequence_length=1000,
    bottleneck_size=5,
    bottleneck_duration=10
)

# Expansion scenario
expansion_seqs = simulate_population_expansion(
    n_sequences=20,
    sequence_length=1000,
    expansion_factor=10.0
)

# Compare neutrality tests
bottleneck_results = neutrality_test_suite(bottleneck_seqs)
expansion_results = neutrality_test_suite(expansion_seqs)

print(f"Bottleneck D: {bottleneck_results['tajimas_d']:.4f}")
print(f"Expansion D: {expansion_results['tajimas_d']:.4f}")
# Both should be negative
```

### 3. Generating Training Data

Create datasets for machine learning:

```python
from metainformant.simulation.popgen import generate_genotype_matrix

# Generate large genotype matrix
genotypes = generate_genotype_matrix(
    n_individuals=1000,
    n_sites=10000,
    min_maf=0.05,
    max_maf=0.5,
    hwe=True
)

# Use for GWAS simulation, PCA, etc.
from metainformant.gwas.structure import compute_pca
pca_result = compute_pca(genotypes, n_components=10)
```

### 4. Simulating Fst Scenarios

Generate populations with different levels of differentiation:

```python
from metainformant.simulation.popgen import generate_two_populations
from metainformant.dna.population import hudson_fst

# Low differentiation
pop1_low, pop2_low = generate_two_populations(
    n_pop1=10, n_pop2=10, sequence_length=1000, fst=0.05
)
fst_low = hudson_fst(pop1_low, pop2_low)

# High differentiation
pop1_high, pop2_high = generate_two_populations(
    n_pop1=10, n_pop2=10, sequence_length=1000, fst=0.5
)
fst_high = hudson_fst(pop1_high, pop2_high)

print(f"Low Fst: {fst_low:.4f}")
print(f"High Fst: {fst_high:.4f}")
```

## Integration with Analysis

### Full Workflow Example

```python
from metainformant.simulation.popgen import generate_two_populations
from metainformant.dna.population_analysis import compare_populations
from metainformant.math.demography import bottleneck_effective_size

# Generate populations
pop1, pop2 = generate_two_populations(
    n_pop1=20,
    n_pop2=20,
    sequence_length=1000,
    fst=0.2,
    within_pop_diversity=0.01
)

# Analyze
comparison = compare_populations(
    pop1_sequences=pop1,
    pop2_sequences=pop2
)

print(f"Fst: {comparison['fst']:.4f}")
print(f"Differentiation: {comparison['differentiation']}")

# Compare to demographic model
ne_bottleneck = bottleneck_effective_size(
    pre_bottleneck_size=10000,
    bottleneck_size=100,
    bottleneck_duration=5
)
print(f"Effective size under bottleneck: {ne_bottleneck:.2f}")
```

## Best Practices

### 1. Use Seeds for Reproducibility

```python
import random

rng = random.Random(42)

sequences = generate_population_sequences(
    n_sequences=10,
    sequence_length=1000,
    mutation_rate=0.01,
    rng=rng
)
```

### 2. Validate Generated Data

Always check that generated data matches expectations:

```python
from metainformant.simulation.popgen import generate_population_sequences
from metainformant.dna.population import nucleotide_diversity

sequences = generate_population_sequences(
    n_sequences=20,
    sequence_length=1000,
    nucleotide_diversity=0.01
)

# Verify
pi_observed = nucleotide_diversity(sequences)
print(f"Target π: 0.01, Observed π: {pi_observed:.4f}")
assert abs(pi_observed - 0.01) < 0.005  # Allow small variance
```

### 3. Combine Multiple Methods

Combine different simulation methods for complex scenarios:

```python
# Generate ancestral population
ancestral = generate_population_sequences(
    n_sequences=10,
    sequence_length=1000,
    nucleotide_diversity=0.01
)

# Create bottleneck
bottleneck = simulate_bottleneck_population(
    n_sequences=20,
    sequence_length=1000,
    reference_sequence=ancestral[0],
    bottleneck_size=5
)
```

## Limitations

1. **Simplified models**: These simulations use simplified models. For publication-quality analysis, consider more sophisticated simulators (e.g., msprime, simuPOP).

2. **Coalescent approximation**: The sequence generation uses approximate coalescent models. Exact coalescent simulations may give different results.

3. **Linkage**: LD simulation is simplified. For realistic LD patterns, use specialized tools.

4. **Selection**: Current implementation focuses on neutral evolution. Selection scenarios require additional parameters.

## References

- **Hudson, R. R. (2002)**. Generating samples under a Wright-Fisher neutral model of genetic variation. *Bioinformatics*, 18(2), 337-338.

- **Hernandez, R. D. (2008)**. A flexible forward simulator for populations subject to selection and demography. *Bioinformatics*, 24(23), 2786-2787.

- **Kelleher, J., et al. (2016)**. Efficient coalescent simulation and genealogical analysis for large sample sizes. *PLoS Computational Biology*, 12(5), e1004842.

## See Also

- [`dna.population`](../dna/population.md) - Population genetics analysis
- [`dna.population_analysis`](../dna/population_workflows.md) - Analysis workflows
- [`simulation.sequences`](sequences.md) - Basic sequence generation
- [`math.demography`](../math/demography.md) - Demographic models

