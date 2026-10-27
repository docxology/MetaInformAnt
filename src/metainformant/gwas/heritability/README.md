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
partitioned = estimation.partitioned_h2(chi2_stats, ld_scores_by_category, n_samples=50000)
rg = estimation.genetic_correlation(z1, z2, ld_scores, n1=50000, n2=40000)
h2_liability = estimation.compute_liability_h2(h2_observed, prevalence=0.01, sample_prevalence=0.5)
```

- `partitioned_h2` takes `ld_scores_by_category`, a dictionary from category name (for example `"promoter"` or `"coding"`) to that category's LD scores, one per SNP.
- `genetic_correlation` takes the Z-scores `z1` and `z2` of the two traits and their sample sizes `n1` and `n2`; `n_overlap` (default 0) is the number of samples the two studies share.
- `compute_liability_h2` takes the population `prevalence` and the proportion of cases in the study, `sample_prevalence`.
