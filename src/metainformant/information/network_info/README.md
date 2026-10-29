# Network Information Flow

Information-theoretic analysis for network and time series data, including transfer entropy, Granger causality, and the Von Neumann entropy of a network, using discretization-based estimation.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `information_flow` module |
| `information_flow.py` | Transfer entropy, Granger causality, network entropy, information flow networks |

## Key Functions

| Function | Description |
|----------|-------------|
| `transfer_entropy()` | Compute transfer entropy between two time series, with a permutation p-value (`random_state` seeds the permutations) |
| `granger_causality()` | Test Granger causality between time series variables |
| `network_entropy()` | Von Neumann entropy of a weighted adjacency matrix: the entropy of the eigenvalues of rho = L / trace(L) |
| `information_flow_network()` | Construct directed information flow network from a dict of time series (`random_state` seeds the transfer-entropy permutation tests) |
| `mutual_information_network()` | Build undirected MI-based network from a data matrix and variable names (`random_state` seeds the permutation tests) |

## Usage

```python
from metainformant.information.network_info import information_flow

te = information_flow.transfer_entropy(source_ts, target_ts, lag=1, random_state=0)
gc = information_flow.granger_causality(series_x, series_y, max_lag=5)
ne = information_flow.network_entropy(adjacency_matrix)

# series_by_name: {"A": [...], "B": [...]}, equal-length series
flow_net = information_flow.information_flow_network(series_by_name, method="transfer_entropy", random_state=0)

# data_matrix: observations x variables; variable_names: one name per column
mi_net = information_flow.mutual_information_network(data_matrix, variable_names, random_state=0)
```

## Network entropy

`network_entropy(adjacency_matrix)` takes a symmetric weighted adjacency matrix (a list of lists) and returns a dictionary
with `entropy` (bits), `normalized_entropy` (the entropy divided by `log2(n_nodes)`), `n_nodes` and `density` (the
fraction of node pairs with a non-zero entry).

- The entropy is the Von Neumann entropy of the density matrix `rho = L / trace(L)` (Braunstein, Ghosh and Severini,
  2006), where `L = D - W` is the Laplacian, `W` holds the absolute values of the edge weights, self-loops are ignored
  and `D` holds the weighted degrees: `S = -sum(lambda_i * log2(lambda_i))` over the eigenvalues of `rho`.
- The eigenvalues come from an eigendecomposition: `numpy.linalg.eigvalsh` when NumPy is available, otherwise a
  pure-Python Jacobi solver. The result depends on the structure of the network, not only on its degree distribution.
- A complete graph on `n` nodes has `S = log2(n - 1)`. A network without edges has no density matrix and gets 0.0, and
  so does the empty matrix.
- A matrix that is not square, or not symmetric (beyond a tolerance of `1e-9 * max(1, largest absolute entry)`),
  raises `ValueError`.

## Seeds

The permutation tests of `transfer_entropy`, `information_flow_network` and `mutual_information_network` take
`random_state`: an integer seed or a `random.Random`. The same seed gives the same p-values, effective transfer entropy
and edges; `None` draws fresh entropy; the global `random` state is neither seeded nor consumed.
`information_flow_network` shares one generator across all node pairs, so the whole network is reproducible. Granger
causality is a deterministic F-test: `granger_causality` takes no seed, and `information_flow_network` leaves
`random_state` unused when `method="granger_causality"`.

## Not to be confused with

`metainformant.information.integration.networks.network_entropy(graph, attribute=None, method="shannon")` is a different
function. It takes a NetworkX graph (or a NumPy adjacency array) and returns a single float: by default the Shannon
entropy of the degree distribution (or of a node attribute), or, with `method="von_neumann"`, the same Von Neumann
entropy of the Laplacian density matrix `rho = L / trace(L)` that the function documented here computes (it calls this
function; a complete graph on `n` nodes gives `log2(n - 1)`). The `von_neumann` method is defined for undirected graphs
and raises `ValueError` for a directed one. The function documented here lives in `information_flow` and returns the
dictionary described above. Import them by module path (`network_info.information_flow` or `integration.networks`) so
the two names do not collide.
