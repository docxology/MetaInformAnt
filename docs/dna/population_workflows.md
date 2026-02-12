# DNA: Population Genetics Workflows

This guide demonstrates common workflows and orchestration patterns for population genetics analysis using METAINFORMANT.

## Overview

The `dna.population_analysis` module provides high-level orchestrators that combine multiple population genetics statistics into convenient workflows. These functions are designed to simplify common analysis tasks while maintaining flexibility.

## Quick Start

```python
from metainformant.dna.population_analysis import (
    calculate_summary_statistics,
    compare_populations,
    neutrality_test_suite,
)
```

## Workflows

### 1. Comprehensive Summary Statistics

Calculate multiple population genetics statistics in a single call:

```python
from metainformant.dna.population_analysis import calculate_summary_statistics

# Using sequence data
sequences = ["AAAA", "AAAT", "AATT", "ATTT"]
stats = calculate_summary_statistics(sequences=sequences)

print(f"Nucleotide diversity (π): {stats['nucleotide_diversity']:.4f}")
print(f"Segregating sites: {stats['segregating_sites']}")
print(f"Watterson's theta: {stats['wattersons_theta']:.4f}")
print(f"Tajima's D: {stats['tajimas_d']:.4f}")
```

**Output includes:**
- `nucleotide_diversity`: Average pairwise differences per site (π)
- `segregating_sites`: Number of polymorphic sites
- `wattersons_theta`: Watterson's theta estimator
- `tajimas_d`: Tajima's D statistic
- `sample_size`: Number of sequences
- `sequence_length`: Length of sequences

**Using genotype matrices:**

```python
genotypes = [
    [0, 1, 0],  # Individual 1
    [0, 1, 1],  # Individual 2
    [1, 0, 1],  # Individual 3
]

stats = calculate_summary_statistics(genotype_matrix=genotypes)
print(f"Allele frequencies: {stats['allele_frequencies']}")
print(f"Observed heterozygosity: {stats['observed_heterozygosity']:.4f}")
```

### 2. Population Comparison

Compare two populations using multiple statistics:

```python
from metainformant.dna.population_analysis import compare_populations

# Define two populations
pop1_sequences = ["AAAA", "AAAA", "AAAT"]
pop2_sequences = ["TTTT", "TTTT", "TTTA"]

# Compare populations
comparison = compare_populations(
    pop1_sequences=pop1_sequences,
    pop2_sequences=pop2_sequences
)

print(f"Fst: {comparison['fst']:.4f}")
print(f"Differentiation: {comparison['differentiation']}")

# Access individual population statistics
print(f"Pop1 diversity: {comparison['pop1_stats']['nucleotide_diversity']:.4f}")
print(f"Pop2 diversity: {comparison['pop2_stats']['nucleotide_diversity']:.4f}")
```

**Differentiation levels:**
- `"none"`: Fst < 0.05
- `"low"`: 0.05 ≤ Fst < 0.15
- `"moderate"`: 0.15 ≤ Fst < 0.25
- `"high"`: Fst ≥ 0.25

### 3. Neutrality Test Suite

Run a comprehensive suite of neutrality tests:

```python
from metainformant.dna.population_analysis import neutrality_test_suite

sequences = ["AAAA", "AAAT", "AATT", "ATTT", "TTTT"]
results = neutrality_test_suite(sequences)

print(f"Tajima's D: {results['tajimas_d']:.4f}")
print(f"π/θ ratio: {results['pi_theta_ratio']:.4f}")
print(f"Interpretation: {results['interpretation']}")
```

**Interpretation values:**
- `"strong_negative_d"`: D < -2.0 (suggests population expansion or purifying selection)
- `"negative_d"`: -2.0 ≤ D < -1.0
- `"neutral"`: -1.0 ≤ D ≤ 1.0 (consistent with neutral evolution)
- `"positive_d"`: 1.0 < D ≤ 2.0
- `"strong_positive_d"`: D > 2.0 (suggests balancing selection or population contraction)

## Integration Examples

### Combining with Math Module

```python
from metainformant.dna.population_analysis import calculate_summary_statistics
from metainformant.math.population_genetics.coalescent import tajimas_D, wattersons_theta

# Calculate basic statistics
sequences = ["AAAA", "AAAT", "AATT"]
stats = calculate_summary_statistics(sequences=sequences)

# Use math module for more advanced calculations
theta_w = wattersons_theta(
    num_segregating_sites=stats["segregating_sites"],
    sample_size=stats["sample_size"]
)

# Full Tajima's D with variance (for publication-quality analysis)
tajima_d_full = tajimas_D(
    sequences=sequences,
    sample_size=stats["sample_size"]
)
```

### Combining with GWAS Structure Analysis

```python
from metainformant.dna.population_analysis import compare_populations
from metainformant.gwas.structure import compute_pca, compute_kinship_matrix

# Compare populations
pop1 = ["AAAA", "AAAA", "AAAT"]
pop2 = ["TTTT", "TTTT", "TTTA"]
comparison = compare_populations(
    pop1_sequences=pop1,
    pop2_sequences=pop2
)

# Convert sequences to genotype matrix for PCA
# (In practice, you'd use real genotype data)
genotype_matrix = [
    [0, 0, 0, 0],  # Convert sequences to 0/1/2 encoding
    [0, 0, 0, 0],
    [0, 0, 0, 1],
    [1, 1, 1, 1],
    [1, 1, 1, 1],
    [1, 1, 1, 0],
]

# PCA analysis
pca_result = compute_pca(genotype_matrix, n_components=2)

# Kinship analysis
kinship_result = compute_kinship_matrix(genotype_matrix, method="vanraden")
```

## Best Practices

### 1. Sequence Alignment

Always ensure sequences are aligned before analysis:

```python
from metainformant.dna.alignment import align_sequences

# Align sequences first
unaligned = ["AAAA", "AAAT", "AATT"]
aligned = align_sequences(unaligned)

# Then calculate statistics
stats = calculate_summary_statistics(sequences=aligned)
```

### 2. Handling Missing Data

The orchestrators handle sequence length mismatches automatically (with warnings):

```python
import warnings

# Sequences with different lengths
sequences = ["AAAA", "AAAT", "AAT"]  # Last sequence is shorter

# Warning will be issued, but analysis proceeds
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    stats = calculate_summary_statistics(sequences=sequences)
```

### 3. Sample Size Considerations

Many statistics require sufficient sample size:

```python
sequences = ["AAAA", "AAAT"]  # Only 2 sequences

stats = calculate_summary_statistics(sequences=sequences)
# Some statistics may be less reliable with small samples
```

### 4. Quality Control

Check data quality before analysis:

```python
def validate_sequences(sequences):
    """Validate sequences before analysis."""
    if len(sequences) < 2:
        raise ValueError("Need at least 2 sequences")
    
    lengths = [len(s) for s in sequences]
    if len(set(lengths)) > 1:
        print(f"Warning: Sequence length mismatch (min={min(lengths)}, max={max(lengths)})")
    
    return True

# Validate before analysis
if validate_sequences(sequences):
    stats = calculate_summary_statistics(sequences=sequences)
```

## Advanced Usage

### Custom Workflows

Combine multiple functions for custom analyses:

```python
from metainformant.dna.population import (
    nucleotide_diversity,
    segregating_sites,
    hudson_fst,
)
from metainformant.math.population_genetics.coalescent import expected_segregating_sites

def custom_diversity_analysis(sequences):
    """Custom workflow combining multiple statistics."""
    pi = nucleotide_diversity(sequences)
    S = segregating_sites(sequences)
    n = len(sequences)
    
    # Calculate expected vs observed
    theta = 0.01  # Assumed mutation parameter
    S_expected = expected_segregating_sites(
        sample_size=n,
        theta=theta,
        sequence_length=len(sequences[0])
    )
    
    return {
        "observed_segregating": S,
        "expected_segregating": S_expected,
        "nucleotide_diversity": pi,
        "ratio": S / S_expected if S_expected > 0 else 0,
    }
```

### Batch Processing

Process multiple populations:

```python
def analyze_multiple_populations(population_dict):
    """Analyze multiple populations."""
    results = {}
    
    for pop_name, sequences in population_dict.items():
        results[pop_name] = calculate_summary_statistics(sequences=sequences)
    
    return results

# Example usage
populations = {
    "pop1": ["AAAA", "AAAT", "AATT"],
    "pop2": ["TTTT", "TTTA", "TTAA"],
    "pop3": ["GGGG", "GGGA", "GGAA"],
}

all_results = analyze_multiple_populations(populations)
```

## References

- **Tajima, F. (1989)**. Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics*, 123(3), 585-595.

- **Watterson, G. A. (1975)**. On the number of segregating sites in genetical models without recombination. *Theoretical Population Biology*, 7(2), 256-276.

- **Hudson, R. R., Slatkin, M., & Maddison, W. P. (1992)**. Estimation of levels of gene flow from DNA sequence data. *Genetics*, 132(2), 583-589.

## See Also

- [`dna.population`](population.md) - Core population genetics functions
- [`math.population_genetics.coalescent`](../math/coalescent.md) - Coalescent theory and advanced statistics
- [`math.fst`](../math/fst.md) - F-statistics and population differentiation
- [`gwas.structure`](../gwas/structure.md) - Population structure analysis

