# Specification: integration

## 🎯 Scope
Phenotype integration with other -omics modules.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 2 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `cross_omic.py`

## Behavioral Contracts
- P-values are two-sided and come from `scipy.stats`: Student's t with `n - 2` degrees of freedom for Pearson correlations, the F distribution with 1 and `n - 2` degrees of freedom for the linear-regression slope, and Pearson's test on average ranks (ties share the mean rank) for Spearman correlations, as `scipy.stats.spearmanr` does.
- A variant or gene with fewer than 3 complete samples (5 for `phenotype_environment_interaction`) is skipped; a constant variable gives a correlation (or slope) of 0 and a p-value of 1.
