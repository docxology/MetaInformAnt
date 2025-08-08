### Simulation: RNA Counts

Function: `simulate_counts_negative_binomial`

```mermaid
flowchart LR
  A[NB mean, dispersion] --> B[simulate_counts_negative_binomial]
  B --> C[genes x samples matrix]
```

Example

```python
from metainformant.simulation import simulate_counts_negative_binomial

matrix = simulate_counts_negative_binomial(100, 6, mean_expression=100.0, dispersion=0.1)
```


