# Specification: network_info

## 🎯 Scope
Network information flow subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 2 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `information_flow.py`

### Behavior Contracts

#### `information_flow.py`
- `granger_causality(source, target, max_lag=5)` — nested OLS F-test on residual
  SS of autoregressive target models with and without lagged source. Candidate
  lags are scored by the BIC of the unrestricted model
  (`n_obs * log(rss/n_obs) + (2*lag + 1) * log(n_obs)`); raw residual SS is
  monotone in lag, so the unpenalised minimum is never used. p-values use
  `scipy.stats.f.sf` when scipy is available. Raises `ValueError` when no
  candidate lag is fittable.
- `_ols_rss` solves `y ~ 1 + x` with `numpy.linalg.lstsq` (NumPy) or exact
  Gaussian elimination with partial pivoting (pure Python); no iterative
  Gauss-Seidel.
