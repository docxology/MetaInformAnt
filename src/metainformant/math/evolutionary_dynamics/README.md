# Evolutionary Dynamics

This submodule models population dynamics and ecological interactions.

## Components

- **core.py**: Discrete and continuous time dynamics (Logistic map, Lotka-Volterra).
- **egt.py**: Evolutionary Game Theory (Replicator dynamics).

## Usage

```python
from metainformant.math.evolutionary_dynamics import lotka_volterra_step

# Simulate receptor-prey step
prey, pred = lotka_volterra_step(prey=100, predator=10)
```
