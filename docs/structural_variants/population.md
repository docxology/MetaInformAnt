# Structural Variants: Population Analysis

Population-scale SV genotyping, allele frequency, association testing, and structure analysis.

## Concept Overview

Moving from individual variant detection to population-level analysis requires
genotyping SVs across many samples, computing allele frequencies, and testing
for associations with phenotypes. This submodule provides:

- **Population Genotyping** -- Genotype known SVs across a sample cohort
- **Allele Frequency** -- Compute SV allele frequencies with optional stratification
- **Association Testing** -- Test SV-phenotype associations (chi-square, logistic regression)
- **Population Structure** -- PCA-based structure analysis from SV genotypes
- **LD Analysis** -- Linkage disequilibrium between SVs and SNPs
- **Callset Merging** -- Combine SV calls across multiple samples

## Function Reference

### genotype_sv_population

Genotype a set of known SVs across multiple samples using depth or split-read evidence.

```python
from metainformant.structural_variants.population.sv_population import (
    genotype_sv_population,
)

genotype_matrix = genotype_sv_population(
    sv_calls=known_variants,       # List of known SVs
    samples=sample_data,           # Dict[str, sample_data]
    method="depth"                 # "depth" or "split_read"
)
# genotype_matrix: np.ndarray (n_variants x n_samples)
# Values: 0 (hom-ref), 1 (het), 2 (hom-alt)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sv_calls` | `List` | required | Known structural variants to genotype |
| `samples` | `Dict` | required | Sample ID to alignment/depth data mapping |
| `method` | `str` | `"depth"` | Genotyping method: `"depth"` or `"split_read"` |

### sv_allele_frequency

Compute allele frequencies from a genotype matrix, optionally stratified by population.

```python
from metainformant.structural_variants.population.sv_population import (
    sv_allele_frequency,
)

freq = sv_allele_frequency(
    genotype_matrix=gt_matrix,
    sample_labels=["EUR", "AFR", "EAS", "EUR", "AFR"]  # Optional population labels
)
# freq["overall_af"]       -> np.ndarray of per-variant AFs
# freq["population_af"]    -> Dict[pop, np.ndarray] per-population AFs
# freq["het_count"]        -> np.ndarray heterozygote counts
# freq["hom_alt_count"]    -> np.ndarray homozygous alt counts
```

### sv_association_test

Test for association between SV genotypes and a phenotype.

```python
from metainformant.structural_variants.population.sv_population import (
    sv_association_test,
)

results = sv_association_test(
    genotypes=gt_matrix,           # (n_variants x n_samples)
    phenotypes=phenotype_vector,   # np.ndarray (n_samples,)
    covariates=covariate_matrix    # Optional (n_samples x n_covariates)
)
# results["p_values"]       -> np.ndarray per-variant p-values
# results["effect_sizes"]   -> np.ndarray (odds ratios or betas)
# results["test_statistics"] -> np.ndarray
# results["method"]         -> str ("chi_square" or "logistic")
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `genotypes` | `np.ndarray` | required | Genotype matrix |
| `phenotypes` | `np.ndarray` | required | Binary or continuous phenotype |
| `covariates` | `np.ndarray` | `None` | Adjustment covariates |

### sv_population_structure

PCA-based population structure analysis from SV genotypes.

```python
from metainformant.structural_variants.population.sv_population import (
    sv_population_structure,
)

pca_result = sv_population_structure(
    genotype_matrix=gt_matrix,
    n_components=10
)
# pca_result["components"]          -> np.ndarray (n_samples x n_components)
# pca_result["explained_variance"]  -> np.ndarray variance per PC
# pca_result["loadings"]            -> np.ndarray variant loadings
```

### sv_ld_analysis

Compute linkage disequilibrium between SVs and nearby SNPs.

```python
from metainformant.structural_variants.population.sv_population import (
    sv_ld_analysis,
)

ld = sv_ld_analysis(
    genotype_matrix_sv=sv_gt,          # SV genotypes (n_sv x n_samples)
    genotype_matrix_snp=snp_gt,        # SNP genotypes (n_snp x n_samples)
    sv_positions=sv_pos,               # List[Tuple[str, int]] (chrom, pos)
    snp_positions=snp_pos              # List[Tuple[str, int]] (chrom, pos)
)
# ld["r2_matrix"]      -> np.ndarray (n_sv x n_snp) pairwise r-squared
# ld["d_prime_matrix"] -> np.ndarray (n_sv x n_snp) D-prime
# ld["tag_snps"]       -> Dict[sv_idx, snp_idx] best tagging SNP per SV
```

### merge_sv_callsets

Merge SV calls from multiple samples into a unified callset.

```python
from metainformant.structural_variants.population.sv_population import (
    merge_sv_callsets,
)

merged = merge_sv_callsets(
    callsets={"sample1": calls1, "sample2": calls2, "sample3": calls3},
    min_overlap=0.5,
    max_breakpoint_distance=500
)
```

## End-to-End Population Analysis

```python
from metainformant.structural_variants.population.sv_population import (
    genotype_sv_population,
    sv_allele_frequency,
    sv_association_test,
    sv_population_structure,
)

# Step 1: Genotype known SVs across cohort
gt_matrix = genotype_sv_population(known_svs, sample_data, method="depth")

# Step 2: Compute allele frequencies
freq = sv_allele_frequency(gt_matrix, sample_labels=populations)
print(f"Mean AF: {freq['overall_af'].mean():.4f}")

# Step 3: Population structure
pca = sv_population_structure(gt_matrix, n_components=5)
print(f"Variance explained (PC1-3): {pca['explained_variance'][:3]}")

# Step 4: Association testing
assoc = sv_association_test(
    gt_matrix,
    phenotype_vector,
    covariates=pca["components"][:, :3]  # Adjust for population structure
)

# Report significant associations
import numpy as np
significant = np.where(assoc["p_values"] < 5e-8)[0]
for idx in significant:
    print(f"SV {idx}: p={assoc['p_values'][idx]:.2e}, "
          f"OR={assoc['effect_sizes'][idx]:.3f}")
```

## Configuration

Population analysis parameters use the `SV_` environment prefix:

| Parameter | Default | Env Variable |
|-----------|---------|-------------|
| `n_components` | `10` | `SV_PCA_COMPONENTS` |
| `min_overlap` | `0.5` | `SV_MIN_OVERLAP` |
| `max_breakpoint_distance` | `500` | `SV_MAX_BP_DISTANCE` |

## See Also

- [Detection](detection.md) -- Generating per-sample SV calls
- [Annotation](annotation.md) -- Functional annotation of population SVs
- [Filtering](filtering.md) -- Filtering before population analysis
