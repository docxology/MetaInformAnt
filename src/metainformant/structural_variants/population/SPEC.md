# Specification: population

## 🎯 Scope
Population-scale structural variant analysis submodule.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 2 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `sv_population.py`

### Function contracts (behavior pinned by tests)

- `sv_population.sv_population_structure(genotype_matrix, n_components)` returns
  `variance_explained[k] = eigenvalue[k] / total variance`, where the total variance is the sum of
  all eigenvalues of the covariance matrix (rounding noise below zero counts as zero), not the
  sum of the returned ones. The leading components therefore sum to less than 1 when
  `n_components` is below the rank of the data, and to 1 when all of them are returned.
