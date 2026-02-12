# Phenotype GWAS Integration

Phenome-wide association study (PheWAS) and GWAS-phenotype integration tools. Tests genetic variants against multiple phenotypes simultaneously, computes genetic risk scores, and screens phenotypes for heritability.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `phewas` module |
| `phewas.py` | PheWAS scans, correlation, risk scores, and heritability screening |

## Key Functions

| Function | Description |
|----------|-------------|
| `run_phewas()` | Run phenome-wide association scan for a variant |
| `phenotype_correlation_matrix()` | Compute pairwise correlation matrix across phenotypes |
| `genetic_risk_score()` | Calculate polygenic/genetic risk scores |
| `phenotype_heritability_screen()` | Screen phenotypes for heritability estimates |
| `categorize_phenotypes()` | Classify phenotypes into functional categories |

## Usage

```python
from metainformant.phenotype.gwas_integration import phewas

results = phewas.run_phewas(variant_genotypes, phenotype_data)
grs = phewas.genetic_risk_score(weights, genotypes)
```
