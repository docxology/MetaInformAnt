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
| `sv_population.sv_population_structure()` | PCA-based population structure from SV genotypes |
| `sv_population.sv_ld_analysis()` | Compute LD between SVs and nearby SNPs |
| `sv_population.merge_sv_callsets()` | Merge multi-sample SV callsets for population analysis |

## Usage

```python
from metainformant.structural_variants.population import sv_population

genotypes = sv_population.genotype_sv_population(sv_calls, samples)
freqs = sv_population.sv_allele_frequency(genotypes)
assoc = sv_population.sv_association_test(genotypes, phenotypes)
pca = sv_population.sv_population_structure(genotypes, n_components=10)
ld = sv_population.sv_ld_analysis(sv_genotypes, snp_genotypes)
```
