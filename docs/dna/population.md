### DNA: Population Genetics

Population genetic analysis and diversity metrics for DNA sequences.

**Functions**: `allele_frequencies`, `observed_heterozygosity`, `nucleotide_diversity`, `tajimas_d`, `hudson_fst`, `segregating_sites`, `wattersons_theta`.

```mermaid
graph TD
  A[Sequences/Genotypes] --> B[Summary Stats]
  B --> C[pi]
  B --> D[Tajima's D]
  B --> E[Fst]
  B --> F[S, theta_w]
```

## Overview

This module provides population genetics statistics for DNA sequence data, including:
- **Diversity metrics**: Nucleotide diversity (π), segregating sites
- **Neutrality tests**: Tajima's D (simplified version)
- **Population differentiation**: F-statistics (Hudson's Fst)
- **Allele frequencies**: From genotype matrices
- **Heterozygosity**: Observed heterozygosity from diploid genotypes

## Functions

### `allele_frequencies(genotype_matrix)`

Compute frequency of allele '1' per site from a genotype matrix.

**Parameters**:
- `genotype_matrix`: List of lists where rows are individuals, columns are sites; values 0/1

**Returns**: List of frequencies per site (float values in [0, 1])

**Example**:
```python
from metainformant.dna import population

genotypes = [
    [0, 1, 0],  # Individual 1
    [0, 1, 1],  # Individual 2
    [1, 0, 1],  # Individual 3
]
freqs = population.allele_frequencies(genotypes)
# Returns: [0.333..., 0.666..., 0.666...]
```

### `observed_heterozygosity(genotypes)`

Calculate the proportion of heterozygous individuals among diploid genotypes.

**Parameters**:
- `genotypes`: Iterable of tuples (a1, a2) with alleles 0/1

**Returns**: Float in [0, 1] representing proportion of heterozygotes

**Example**:
```python
genotypes = [(0, 0), (0, 1), (1, 1), (1, 0)]
h_obs = population.observed_heterozygosity(genotypes)
# Returns: 0.5 (2 of 4 are heterozygous)
```

### `nucleotide_diversity(seqs)`

Calculate average pairwise nucleotide difference per site (π).

**Parameters**:
- `seqs`: Sequence of DNA sequences (strings)

**Returns**: Average pairwise diversity per site (float)

**Note**: If sequences have different lengths, truncates to shortest length.

**Example**:
```python
seqs = ["AAAA", "AAAT", "AATT"]
pi = population.nucleotide_diversity(seqs)
# Returns average pairwise differences per site
```

### `tajimas_d(seqs)`

Calculate simplified Tajima's D statistic for departure from neutral equilibrium.

**Parameters**:
- `seqs`: Sequence of DNA sequences (strings)

**Returns**: Simplified Tajima's D statistic (float)

**Note**: This is a simplified version. For publication-quality analysis with proper variance calculation, use `metainformant.math.coalescent.tajimas_D()`.

**Interpretation**:
- D > 0: Suggests balancing selection or population contraction
- D < 0: Suggests directional selection or population expansion
- D ≈ 0: Consistent with neutral equilibrium

**Example**:
```python
seqs = ["AAAA", "AAAT", "AATT", "ATTT"]
d = population.tajimas_d(seqs)
```

### `hudson_fst(pop1, pop2)`

Calculate Hudson's Fst estimator for two populations.

**Parameters**:
- `pop1`: Sequence of DNA sequences from population 1
- `pop2`: Sequence of DNA sequences from population 2

**Returns**: Fst value in [0, 1] where:
- 0 = no differentiation
- 1 = complete differentiation (fixed differences)

**Example**:
```python
pop1 = ["AAAA", "AAAA", "AAAT"]
pop2 = ["TTTT", "TTTT", "TTTA"]
fst = population.hudson_fst(pop1, pop2)
# Returns: 1.0 for fixed differences
```

### `segregating_sites(seqs)`

Count the number of sites with more than one allele among sequences.

**Parameters**:
- `seqs`: Sequence of DNA sequences (strings)

**Returns**: Integer count of segregating (polymorphic) sites

**Example**:
```python
seqs = ["AAAA", "AAAT", "AATT"]
S = population.segregating_sites(seqs)
# Returns: 2 (sites at positions 2 and 3 are polymorphic)
```

### `wattersons_theta(seqs)`

Calculate Watterson's theta per site: θ_W = S / (a₁ × L).

**Parameters**:
- `seqs`: Sequence of DNA sequences (strings)

**Returns**: Watterson's theta estimate per site (float)

**Formula**: θ_W = S / (a₁ × L) where:
- S = number of segregating sites
- a₁ = Σᵢ₌₁ⁿ⁻¹ 1/i (harmonic sum)
- L = sequence length

**Example**:
```python
seqs = ["AAAA", "AAAT", "AATT"]
theta_w = population.wattersons_theta(seqs)
```

## Usage Examples

### Basic Population Analysis

```python
from metainformant.dna import population

# Load sequences (from FASTA, etc.)
sequences = ["ATCGATCG", "ATCGTTCG", "ATCGATCG", "ATCGTTCG"]

# Calculate diversity metrics
pi = population.nucleotide_diversity(sequences)
print(f"Nucleotide diversity (π): {pi:.6f}")

# Count segregating sites
S = population.segregating_sites(sequences)
print(f"Segregating sites: {S}")

# Calculate Watterson's theta
theta_w = population.wattersons_theta(sequences)
print(f"Watterson's theta: {theta_w:.6f}")

# Calculate Tajima's D
d = population.tajimas_d(sequences)
print(f"Tajima's D: {d:.6f}")
```

### Population Differentiation

```python
# Compare two populations
pop1 = ["AAAA", "AAAT", "AAAA"]
pop2 = ["TTTT", "TTTA", "TTTT"]

fst = population.hudson_fst(pop1, pop2)
print(f"Fst between populations: {fst:.4f}")
```

### Genotype Analysis

```python
# From genotype matrix (0/1 encoding)
genotype_matrix = [
    [0, 1, 0],  # Individual 1
    [0, 1, 1],  # Individual 2
    [1, 0, 1],  # Individual 3
    [1, 1, 0],  # Individual 4
]

# Calculate allele frequencies
freqs = population.allele_frequencies(genotype_matrix)
print(f"Allele frequencies: {freqs}")

# Calculate observed heterozygosity
genotypes = [(0, 0), (0, 1), (1, 1), (1, 0)]
h_obs = population.observed_heterozygosity(genotypes)
print(f"Observed heterozygosity: {h_obs:.4f}")
```

## Integration with Other Modules

### With Math Module

For more advanced population genetics calculations, combine with `math.popgen`:

```python
from metainformant.dna import population
from metainformant.math.coalescent import tajimas_D

# Calculate observed statistics
seqs = ["AAAA", "AAAT", "AATT"]
pi = population.nucleotide_diversity(seqs)
S = population.segregating_sites(seqs)

# Use full Tajima's D implementation
d_full = tajimas_D(S, pi, len(seqs))
```

### With GWAS Module

Population genetics statistics complement GWAS analysis:

```python
from metainformant.dna import population
from metainformant.gwas import parse_vcf_full

# Parse VCF and extract sequences
vcf_data = parse_vcf_full("variants.vcf")
# ... convert to sequences ...

# Calculate population diversity
pi = population.nucleotide_diversity(sequences)
```

## Notes and Limitations

1. **Sequence Length**: Functions that process multiple sequences truncate to the shortest length if sequences differ.

2. **Tajima's D**: The `tajimas_d()` function is a simplified version. For statistical testing, use `metainformant.math.coalescent.tajimas_D()`.

3. **Genotype Encoding**: Allele frequency functions assume binary encoding (0/1). Ensure your data is properly encoded.

4. **Input Validation**: Functions handle edge cases (empty inputs, single sequences) but may not validate all input types.

## References

- Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics*, 123(3), 585-595.
- Hudson, R. R., Slatkin, M., & Maddison, W. P. (1992). Estimation of levels of gene flow from DNA sequence data. *Genetics*, 132(2), 583-589.
- Watterson, G. A. (1975). On the number of segregating sites in genetical models without recombination. *Theoretical Population Biology*, 7(2), 256-276.
