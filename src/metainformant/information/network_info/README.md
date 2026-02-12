# Network Information Flow

Information-theoretic analysis for network and time series data, including transfer entropy, Granger causality, and network entropy computation using discretization-based estimation.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `information_flow` module |
| `information_flow.py` | Transfer entropy, Granger causality, network entropy, information flow networks |

## Key Functions

| Function | Description |
|----------|-------------|
| `transfer_entropy()` | Compute transfer entropy between two time series |
| `granger_causality()` | Test Granger causality between time series variables |
| `network_entropy()` | Compute Von Neumann entropy of a network adjacency matrix |
| `information_flow_network()` | Construct directed information flow network from multivariate series |
| `mutual_information_network()` | Build undirected MI-based network from multivariate data |

## Usage

```python
from metainformant.information.network_info import information_flow

te = information_flow.transfer_entropy(source_ts, target_ts, lag=1)
gc = information_flow.granger_causality(series_x, series_y, max_lag=5)
ne = information_flow.network_entropy(adjacency_matrix)
flow_net = information_flow.information_flow_network(multivariate_data)
mi_net = information_flow.mutual_information_network(multivariate_data)
```
