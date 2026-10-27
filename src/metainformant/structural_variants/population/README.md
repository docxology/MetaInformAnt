# Population

Population-scale structural variant analysis providing cohort genotyping, allele frequency computation, association testing, PCA-based population structure, and SV-SNP linkage disequilibrium analysis.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports sv_population submodule |
| `sv_population.py` | Population genotyping, allele frequency, association, PCA, LD |

## Key Functions

| Function | Description |
|----------|-------------|
| `sv_population.genotype_sv_population()` | Population-scale SV genotyping across samples |
| `sv_population.sv_allele_frequency()` | Compute SV allele frequencies in a cohort |
| `sv_population.sv_association_test()` | Test SV association with phenotypes |
| `sv_population.sv_population_structure()` | PCA-based population structure from SV genotypes; each `variance_explained` value is that component's share of the total variance (the sum of all eigenvalues) |
| `sv_population.sv_ld_analysis()` | Compute LD between SVs and nearby SNPs |
| `sv_population.merge_sv_callsets()` | Merge multi-sample SV callsets for population analysis |

## Usage

```python
from metainformant.structural_variants.population import sv_population

genotypes = sv_population.genotype_sv_population(sv_calls, samples)
gt_matrix = genotypes["genotype_matrix"]                       # n_svs x n_samples
freqs = sv_population.sv_allele_frequency(gt_matrix)
assoc = sv_population.sv_association_test(gt_matrix[0], phenotypes)  # one SV at a time
pca = sv_population.sv_population_structure(gt_matrix, n_components=10)
ld = sv_population.sv_ld_analysis(sv_gt_matrix, snp_gt_matrix, sv_positions, snp_positions)
```
