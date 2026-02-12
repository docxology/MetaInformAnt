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
| `phenotype_genotype_association()` | Test associations between phenotypes and genetic variants |
| `trait_expression_correlation()` | Correlate trait values with gene expression levels |
| `multi_phenotype_integration()` | Integrate and analyze multiple phenotype domains together |
| `phenotype_environment_interaction()` | Analyze gene-by-environment interactions on phenotypes |

## Usage

```python
from metainformant.phenotype.integration import cross_omic

results = cross_omic.phenotype_genotype_association(phenotypes, genotypes)
correlations = cross_omic.trait_expression_correlation(traits, expression)
```
