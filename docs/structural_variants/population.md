# Structural Variants: Population Analysis

Population-scale SV genotyping, allele frequency, association testing, and structure analysis.

## Concept Overview

Moving from individual variant detection to population-level analysis requires
genotyping SVs across many samples, computing allele frequencies, and testing
for associations with phenotypes. This submodule provides:

- **Population Genotyping** -- Genotype known SVs across a sample cohort
- **Allele Frequency** -- Compute SV allele frequencies with optional stratification
- **Association Testing** -- Test SV-phenotype associations (linear or logistic regression)
- **Population Structure** -- PCA-based structure analysis from SV genotypes
- **LD Analysis** -- Linkage disequilibrium between SVs and SNPs
- **Callset Merging** -- Combine SV calls across multiple samples

Genotype matrices are `n_svs x n_samples` and hold the number of alternate alleles (0, 1 or 2).
The functions work on in-memory lists and NumPy arrays and open no BAM files.

## Function Reference

### genotype_sv_population

Genotype a set of known SVs across multiple samples using depth or split-read evidence.

```python
from metainformant.structural_variants.population.sv_population import (
    genotype_sv_population,
)

result = genotype_sv_population(
    sv_calls=known_variants,       # SV dicts carrying per-sample evidence (below)
    samples=sample_ids,            # List of sample names
    method="depth"                 # "depth" or "split"
)
genotype_matrix = result["genotype_matrix"]   # list of lists, n_svs x n_samples
# Values: 0 (hom-ref), 1 (het), 2 (hom-alt)
# result["quality_scores"]      -> n_svs x n_samples genotype qualities (0-100)
# result["allele_frequencies"]  -> one allele frequency per SV
```

Each SV dictionary has `chrom`, `start`, `end`, `sv_type` and a `samples` mapping from sample
name to its evidence: `depth_ratio` (observed/expected depth, a linear ratio) for `"depth"`;
`split_reads`, `discordant_pairs` and `total_reads` (default 30) for `"split"`. A sample
without evidence is genotyped as reference.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sv_calls` | `List[Dict]` | required | Known structural variants with per-sample evidence under `samples` |
| `samples` | `List[str]` | required | Sample names, in the column order of the result |
| `method` | `str` | `"depth"` | Genotyping method: `"depth"` or `"split"` (any other value raises `ValueError`) |

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
# freq["frequencies"]            -> allele frequency per SV
# freq["maf"]                    -> minor allele frequency per SV
# freq["n_polymorphic"]          -> number of SVs with MAF > 0
# freq["sfs"]                    -> site frequency spectrum (counts per alternate-allele count)
# freq["population_frequencies"] -> Dict[pop, per-SV frequencies] (only with sample_labels)
```

### sv_association_test

Test for association between the genotypes of one SV and a phenotype: linear regression for a
continuous phenotype, logistic regression for a binary 0/1 phenotype.

```python
from metainformant.structural_variants.population.sv_population import (
    sv_association_test,
)

result = sv_association_test(
    genotypes=gt_matrix[0],        # One SV: one genotype (0, 1, 2) per sample
    phenotypes=phenotype_vector,   # One value per sample
    covariates=covariate_matrix    # Optional (n_samples x n_covariates)
)
# result["beta"], result["se"], result["p_value"]
# result["odds_ratio"] -> exp(beta) for a binary phenotype, None otherwise
# result["n_samples"]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `genotypes` | `List[int]` | required | Genotypes of one SV, one per sample |
| `phenotypes` | `List[float]` | required | Binary (0/1) or continuous phenotype, one per sample |
| `covariates` | `List[List[float]]` | `None` | Adjustment covariates (n_samples x n_covariates) |

Test every SV with a loop over the rows of the genotype matrix.

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
# pca_result["pcs"]               -> (n_samples x n_components) PC coordinates
# pca_result["eigenvalues"]       -> eigenvalue per returned component
# pca_result["variance_explained"] -> each component's share of the total variance
# pca_result["loadings"]          -> (n_svs x n_components) SV loadings
```

`variance_explained[k]` is `eigenvalue[k]` divided by the total variance, the sum of all
eigenvalues of the covariance matrix (rounding noise below zero counts as zero), and not by the
sum of the returned ones. The returned values therefore sum to less than 1 when `n_components`
is below the number of components of the data, and to 1 when all of them are returned
(`n_components` is capped at the smaller of the number of samples and of SVs).

### sv_ld_analysis

Compute linkage disequilibrium (r²) between each SV and every SNP you pass.

```python
from metainformant.structural_variants.population.sv_population import (
    sv_ld_analysis,
)

ld = sv_ld_analysis(
    genotype_matrix_sv=sv_gt,          # SV genotypes (n_sv x n_samples)
    genotype_matrix_snp=snp_gt,        # SNP genotypes (n_snp x n_samples)
    sv_positions=sv_pos,               # List[int], one position per SV
    snp_positions=snp_pos              # List[int], one position per SNP
)
# ld["ld_matrix"]      -> (n_sv x n_snp) r-squared values as lists
# ld["tag_snps"]       -> one dict per SV: snp_index, snp_position, r2 (the best tagging SNP)
# ld["max_r2_per_sv"]  -> highest r-squared per SV
```

There is no window argument: r² is computed against every SNP in `genotype_matrix_snp`.

### merge_sv_callsets

Merge SV calls from multiple samples into a unified callset.

```python
from metainformant.structural_variants.population.sv_population import (
    merge_sv_callsets,
)

merged = merge_sv_callsets(
    callsets=[calls1, calls2, calls3],   # one list of SV dicts per sample or caller
    min_overlap=0.5,
    max_breakpoint_distance=500
)
# Each merged dict has chrom, start, end, sv_type (the consensus call) plus
# n_samples, samples, support_count and merged_calls
```

A call's `sample` key names its sample; calls without one are named `sample_<index of its callset>`.

## End-to-End Population Analysis

```python
import numpy as np

from metainformant.structural_variants.population.sv_population import (
    genotype_sv_population,
    sv_allele_frequency,
    sv_association_test,
    sv_population_structure,
)

# Step 1: Genotype known SVs across the cohort
genotypes = genotype_sv_population(known_svs, sample_ids, method="depth")
gt_matrix = np.asarray(genotypes["genotype_matrix"])

# Step 2: Compute allele frequencies
freq = sv_allele_frequency(gt_matrix, sample_labels=populations)
print(f"Mean AF: {np.mean(freq['frequencies']):.4f}")

# Step 3: Population structure
pca = sv_population_structure(gt_matrix, n_components=5)
print(f"Variance explained (PC1-3): {pca['variance_explained'][:3]}")

# Step 4: Association testing, one SV at a time, adjusted for population structure
results = [
    sv_association_test(row, phenotype_vector, covariates=pca["pcs"][:, :3])
    for row in gt_matrix
]

# Report significant associations
for idx, result in enumerate(results):
    if result["p_value"] < 5e-8:
        print(f"SV {idx}: p={result['p_value']:.2e}, beta={result['beta']:.3f}")
```

## Configuration

Population analysis parameters are function arguments; none has an environment variable:

| Parameter | Default | Function |
|-----------|---------|----------|
| `n_components` | `10` | `sv_population_structure` |
| `min_overlap` | `0.5` | `merge_sv_callsets` |
| `max_breakpoint_distance` | `500` | `merge_sv_callsets` |

See [CONFIGURATION.md](CONFIGURATION.md) for every argument.

## See Also

- [Detection](detection.md) -- Generating per-sample SV calls
- [Annotation](annotation.md) -- Functional annotation of population SVs
- [Filtering](filtering.md) -- Filtering before population analysis
