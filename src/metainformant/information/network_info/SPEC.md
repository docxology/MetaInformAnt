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
- `network_entropy(adjacency_matrix)` — Von Neumann entropy of the density matrix `rho = L / trace(L)` with `L = D - W`,
  where `W` holds the absolute values of the off-diagonal edge weights (self-loops are ignored) and `D` the weighted
  degrees. `S = -sum(lambda * log2(lambda))` over the eigenvalues of `rho` from a full eigendecomposition
  (`numpy.linalg.eigvalsh` with NumPy, a cyclic Jacobi solver in pure Python otherwise; eigenvalues that are zero or
  rounding noise contribute 0). A complete graph on `n` nodes gives `log2(n - 1)`; a matrix without edges gives 0.0.
  Returns `entropy`, `normalized_entropy` (divided by `log2(n_nodes)`), `n_nodes` and `density`. Raises `ValueError` for a
  matrix that is not square or not symmetric (tolerance `1e-9 * max(1, largest absolute entry)`).
- `transfer_entropy(source, target, lag=1, n_bins=10, random_state=None)`,
  `information_flow_network(time_series, method="transfer_entropy", threshold=0.05, random_state=None)` and
  `mutual_information_network(data_matrix, variable_names, n_bins=10, threshold=0.05, random_state=None)` draw their
  permutations (100 for transfer entropy, 50 per pair for the MI network) from a generator built from `random_state`
  (an integer seed or a `random.Random`). The same seed gives the same p-values and edges, `None` draws fresh entropy,
  and the global `random` state is neither seeded nor consumed. `information_flow_network` shares one generator across
  all node pairs; with `method="granger_causality"` (a deterministic F-test) the seed is unused.
- `information_flow.network_entropy` is not `integration.networks.network_entropy(graph, attribute=None,
  method="shannon")`: that function takes a NetworkX graph or a NumPy array and returns a float.
