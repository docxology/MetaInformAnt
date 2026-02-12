# GWAS Heritability Estimation

SNP heritability estimation from GWAS summary statistics using LD Score Regression, partitioned heritability, genetic correlation, and related methods.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `estimation` module |
| `estimation.py` | LDSC, partitioned h2, genetic correlation, HE regression, GREML, liability scale |

## Key Functions

| Function | Description |
|----------|-------------|
| `estimate_h2_ldsc()` | Estimate SNP heritability using LD Score Regression |
| `partitioned_h2()` | Partition heritability by functional annotation categories |
| `genetic_correlation()` | Estimate cross-trait genetic correlation via bivariate LDSC |
| `haseman_elston_regression()` | Heritability estimation via Haseman-Elston regression |
| `greml_simple()` | Simplified Genomic REML heritability estimation |
| `compute_liability_h2()` | Transform observed-scale h2 to liability scale for case-control |

## Usage

```python
from metainformant.gwas.heritability import estimation

h2 = estimation.estimate_h2_ldsc(chi2_stats, ld_scores, n_samples=50000)
partitioned = estimation.partitioned_h2(chi2_stats, ld_scores, annotations)
rg = estimation.genetic_correlation(chi2_trait1, chi2_trait2, ld_scores)
h2_liability = estimation.compute_liability_h2(h2_observed, prevalence=0.01)
```
