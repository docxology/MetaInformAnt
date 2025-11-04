### Math: F-statistics

F-statistics for measuring population differentiation and inbreeding.

**Functions**: `fst_from_heterozygosity`, `fst_from_allele_freqs`.

## Overview

F-statistics (Fst) measure the proportion of genetic diversity that is between
subpopulations rather than within them. High Fst indicates strong population
differentiation, while low Fst suggests gene flow between populations.

## Functions

### `fst_from_heterozygosity(Hs, Ht)`

Calculate Wright's Fst from heterozygosity values.

**Parameters**:
- `Hs`: Average heterozygosity within subpopulations
- `Ht`: Total heterozygosity in the entire population

**Returns**: Fst value in [0, 1]. Returns 0.0 if Ht <= 0.

**Formula**: Fst = (Ht - Hs) / Ht

**Example**:
```python
from metainformant.math import fst_from_heterozygosity

# Within-population heterozygosity
Hs = 0.3
# Total population heterozygosity
Ht = 0.4
fst = fst_from_heterozygosity(Hs, Ht)
# Returns: 0.25 (25% of diversity is between populations)
```

### `fst_from_allele_freqs(subpop_allele_freqs)`

Compute Fst for a biallelic locus from subpopulation allele frequencies.

**Parameters**:
- `subpop_allele_freqs`: Iterable of allele A frequencies (one per subpopulation)

**Returns**: Fst value in [0, 1]

**Formula**: 
- Hs = mean(2p(1-p)) across subpopulations
- Ht = 2p̄(1-p̄) where p̄ is the mean allele frequency
- Fst = (Ht - Hs) / Ht

**Note**: This function assumes biallelic loci.

**Example**:
```python
from metainformant.math import fst_from_allele_freqs

# Allele frequencies in 4 subpopulations
subpop_freqs = [0.2, 0.3, 0.25, 0.28]
fst = fst_from_allele_freqs(subpop_freqs)
```

## Interpretation

- **Fst = 0**: No differentiation (populations panmictic)
- **Fst = 0.01-0.05**: Low differentiation (moderate gene flow)
- **Fst = 0.05-0.15**: Moderate differentiation
- **Fst = 0.15-0.25**: Large differentiation
- **Fst > 0.25**: Very large differentiation (little gene flow)
- **Fst = 1**: Complete differentiation (fixed differences)

## Usage Examples

### Basic Fst Calculation

```python
from metainformant.math import fst_from_heterozygosity, fst_from_allele_freqs

# Method 1: From heterozygosity values
fst1 = fst_from_heterozygosity(Hs=0.3, Ht=0.4)

# Method 2: From allele frequencies
subpop_freqs = [0.2, 0.3, 0.25, 0.28]
fst2 = fst_from_allele_freqs(subpop_freqs)
```

### Integration with DNA Module

```python
from metainformant.dna import population
from metainformant.math import fst_from_heterozygosity

# Calculate heterozygosity within populations
pop1 = ["AAAA", "AAAT", "AAAA"]
pop2 = ["TTTT", "TTTA", "TTTT"]

# Calculate Fst using DNA module
fst_dna = population.hudson_fst(pop1, pop2)

# Or calculate from heterozygosity values
# (requires calculating heterozygosity separately)
```

## References

- Wright, S. (1949). The genetical structure of populations. *Annals of Eugenics*, 15(1), 323-354.
- Weir, B. S., & Cockerham, C. C. (1984). Estimating F-statistics for the analysis of population structure. *Evolution*, 38(6), 1358-1370.

