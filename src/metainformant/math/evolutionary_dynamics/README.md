# Evolutionary Dynamics

Population dynamics models and evolutionary game theory including logistic maps, Lotka-Volterra systems, and replicator dynamics.

## Contents

| File | Purpose |
|------|---------|
| `core.py` | Population dynamics: logistic map, Lotka-Volterra, replicator dynamics |
| `egt.py` | Evolutionary game theory: replicator equation and strategy evolution |

## Key Functions

| Function | Description |
|----------|-------------|
| `logistic_map()` | Generate logistic map sequence for chaos analysis |
| `lotka_volterra_step()` | Single time step of predator-prey dynamics |
| `replicator_derivative()` | Derivative of replicator dynamics for strategy frequencies |
| `replicator_step()` | Discrete step of replicator dynamics evolution |

## Usage

```python
from metainformant.math.evolutionary_dynamics.core import logistic_map, lotka_volterra_step
from metainformant.math.evolutionary_dynamics.egt import replicator_derivative

trajectory = logistic_map(r=3.5, x0=0.5, n_iterations=100)
prey, pred = lotka_volterra_step(prey=100, predator=20)
derivs = replicator_derivative(fitnesses=[1.1, 0.9], frequencies=[0.6, 0.4])
```
