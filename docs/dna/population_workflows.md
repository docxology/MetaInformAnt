# DNA: Population Genetics Workflows

This guide demonstrates common workflows and orchestration patterns for population genetics analysis using METAINFORMANT.

## Overview

The `dna.population.analysis` module provides high-level orchestrators that combine multiple population genetics statistics into convenient workflows. These functions are designed to simplify common analysis tasks while maintaining flexibility.

## Quick Start

```python-snippet
from metainformant.dna.population.analysis import (
    calculate_summary_statistics,
    compare_populations,
    neutrality_test_suite,
)
```

## Workflows

### 1. Comprehensive Summary Statistics

Calculate multiple population genetics statistics in a single call:

```python-snippet
from metainformant.dna.population.analysis import calculate_summary_statistics

# Using sequence data
sequences = ["AAAA", "AAAT", "AATT", "ATTT"]
stats = calculate_summary_statistics(sequences=sequences)

print(f"Nucleotide diversity (π): {stats['nucleotide_diversity']:.4f}")
print(f"Segregating sites: {stats['segregating_sites']}")
print(f"Watterson's theta: {stats['watterson_theta']:.4f}")
print(f"Tajima's D: {stats['tajima_d']:.4f}")
```

**Output includes:**
- `nucleotide_diversity`: Average pairwise differences per site (π)
- `segregating_sites`: Number of polymorphic sites
- `watterson_theta`: Watterson's theta estimator
- `tajima_d`: Tajima's D statistic
- `tajima_d_p_value`: p-value for Tajima's D
- `fu_li_d`, `fu_li_d_p_value`: Fu and Li's D statistic and its p-value
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
print(f"Individuals: {stats['num_individuals']}, loci: {stats['num_loci']}")
print(f"Mean observed heterozygosity: {stats['mean_heterozygosity']:.4f}")
```

**Genotype-matrix output:** `num_individuals`, `num_loci`, `mean_heterozygosity` (mean of the per-locus observed heterozygosity), `heterozygosity_variance`.

### 2. Population Comparison

Compare two populations using multiple statistics:

```python-snippet
from metainformant.dna.population import hudson_fst
from metainformant.dna.population.analysis import calculate_summary_statistics, compare_populations

# Define two populations
pop1_sequences = ["AAAA", "AAAA", "AAAT"]
pop2_sequences = ["TTTT", "TTTT", "TTTA"]

# Summarize each population, then compare the two summaries
pop1_stats = calculate_summary_statistics(sequences=pop1_sequences)
pop2_stats = calculate_summary_statistics(sequences=pop2_sequences)
comparison = compare_populations(pop1_stats, pop2_stats)

# Fst is computed from the sequences directly
fst = hudson_fst(pop1_sequences, pop2_sequences)
print(f"Fst: {fst:.4f}")
print(f"Diversity ratio (pop2/pop1): {comparison['nucleotide_diversity_ratio']:.4f}")
print(f"Tajima's D difference: {comparison['tajima_d_difference']:.4f}")
print(f"Pop1 selection pattern: {comparison['pop1_selection_pattern']}")

# Access individual population statistics
print(f"Pop1 diversity: {pop1_stats['nucleotide_diversity']:.4f}")
print(f"Pop2 diversity: {pop2_stats['nucleotide_diversity']:.4f}")
```

**Comparison output:**
- `nucleotide_diversity_ratio`, `nucleotide_diversity_difference`: π of population 2 relative to population 1
- `tajima_d_difference`: Tajima's D of population 2 minus population 1
- `pop1_selection_pattern`, `pop2_selection_pattern`: `"balancing_selection"` (D > 2), `"positive_selection"` (D < -2), or `"neutral"`
- `fst_difference`, `fst_ratio`: only present when both summaries contain an `fst` entry

### 3. Neutrality Test Suite

Run a comprehensive suite of neutrality tests:

```python-snippet
from metainformant.dna.population.analysis import neutrality_test_suite

sequences = ["AAAA", "AAAT", "AATT", "ATTT", "TTTT"]
results = neutrality_test_suite(sequences)

print(f"Tajima's D: {results['tajima_d']:.4f}")
print(f"Fu and Li's D*: {results['fu_li_d_star']:.4f}")
print(f"Interpretation: {results['neutrality_summary']['tajima_d']}")
```

**Output includes:**
- `tajima_d`, `fu_li_d_star`, `fu_li_f_star`, `fay_wu_h`: test statistics, each with a matching `<name>_p_value` (a statistic is `None` if it could not be computed)
- `num_sequences`, `sequence_length`: Number and length of sequences
- `neutrality_summary`: Interpretation strings for `tajima_d` and `fu_li_d_star`

**Interpretation values** (based on the sign of the statistic):
- `neutrality_summary["tajima_d"]`: `"balancing_selection_or_population_bottleneck"` (D > 0) or `"positive_selection_or_population_expansion"` (D < 0) when the p-value is below `alpha` (default 0.05); otherwise `"neutral_evolution"` (consistent with neutrality); `"calculation_failed"` when D is `None`
- `neutrality_summary["fu_li_d_star"]`: `"balancing_selection"` (D* > 0) or `"positive_selection_or_population_expansion"` (D* < 0) under the same rule

## Integration Examples

### Combining with Math Module

```python-snippet
from metainformant.dna.population.analysis import calculate_summary_statistics
from metainformant.math.population_genetics.coalescent import tajimas_D, watterson_theta

# Calculate basic statistics
sequences = ["AAAA", "AAAT", "AATT"]
stats = calculate_summary_statistics(sequences=sequences)

# Use math module for more advanced calculations
theta_w = watterson_theta(
    num_segregating_sites=stats["segregating_sites"],
    sample_size=stats["sample_size"]
)

# Full Tajima's D with variance (for publication-quality analysis)
tajima_d_full = tajimas_D(
    pi=stats["nucleotide_diversity"],
    s=stats["segregating_sites"],
    n=stats["sample_size"]
)
```

### Combining with GWAS Structure Analysis

```python-snippet
from metainformant.dna.population.analysis import calculate_summary_statistics, compare_populations
from metainformant.gwas.analysis.structure import compute_pca, compute_kinship_matrix

# Compare populations
pop1 = ["AAAA", "AAAA", "AAAT"]
pop2 = ["TTTT", "TTTT", "TTTA"]
comparison = compare_populations(
    calculate_summary_statistics(sequences=pop1),
    calculate_summary_statistics(sequences=pop2)
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

```python-snippet
from metainformant.dna.alignment.msa import align_msa

# Align sequences first (align_msa takes and returns {id: sequence} dicts)
unaligned = {"seq1": "AAAA", "seq2": "AAAT", "seq3": "AATT"}
aligned = align_msa(unaligned)

# Then calculate statistics
stats = calculate_summary_statistics(sequences=list(aligned.values()))
```

### 2. Handling Missing Data

Sequence statistics need an alignment: sequences of different lengths raise
`ValueError`. Within an alignment, gaps and ambiguous bases are skipped:
nucleotide diversity compares each pair over the sites both sequences cover,
and the neutrality tests use only columns where every sequence has A, C, G or T.

```python
sequences = ["AAAA", "AAAT", "AAT"]  # Last sequence is shorter

try:
    calculate_summary_statistics(sequences=sequences)
except ValueError as exc:
    print(exc)  # Sequences must be aligned (same length)
```

### 3. Sample Size Considerations

Many statistics require sufficient sample size:

```python
sequences = ["AAAA", "AAAT"]  # Only 2 sequences

stats = calculate_summary_statistics(sequences=sequences)
# Neutrality tests need at least 4 sequences; below that they report 0.0 with p-value 1.0
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
        n_samples=n,
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

