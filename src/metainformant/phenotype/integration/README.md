# Phenotype Integration

Cross-omic integration connecting phenotype data with DNA variants, gene expression, and environmental factors. Enables genotype-phenotype associations, trait-expression correlations, and gene-environment interaction analysis.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `cross_omic` module |
| `cross_omic.py` | Phenotype-genotype, trait-expression, and GxE integration functions |

## Key Functions

| Function | Description |
|----------|-------------|
| `phenotype_genotype_association()` | Test associations between phenotypes and genetic variants (`method` `"linear_regression"` or `"correlation"`) |
| `trait_expression_correlation()` | Correlate trait values with gene expression levels (`method` `"pearson"` or `"spearman"`) |
| `multi_phenotype_integration()` | Integrate and analyze multiple phenotype domains together |
| `phenotype_environment_interaction()` | Analyze gene-by-environment interactions on phenotypes |

## Usage

```python
from metainformant.phenotype.integration import cross_omic

results = cross_omic.phenotype_genotype_association(phenotypes, genotypes)
correlations = cross_omic.trait_expression_correlation(traits, expression)
```

## Statistics

The p-values come from `scipy.stats` and are two-sided:

| Test | Statistic and distribution |
|------|----------------------------|
| Pearson correlation (`"correlation"`, `"pearson"`, the interaction term of `phenotype_environment_interaction()`) | t = `r * sqrt((n - 2) / (1 - r**2))` with `n - 2` degrees of freedom, as in `scipy.stats.pearsonr`, whatever the sample size |
| Linear regression (`"linear_regression"`) | F = `R**2 / ((1 - R**2) / (n - 2))` with 1 and `n - 2` degrees of freedom, equal to the t-test of the slope (`scipy.stats.linregress`) |
| Spearman correlation (`"spearman"`) | Pearson's test on the ranks; tied values share the mean of the ranks they span, so the result equals `scipy.stats.spearmanr` on data with ties and does not depend on the order of the samples |

A variant or gene needs at least 3 complete samples (5 for `phenotype_environment_interaction()`) and is skipped with fewer. A constant variable gives a correlation (or slope) of 0 and a p-value of 1 rather than an error.
