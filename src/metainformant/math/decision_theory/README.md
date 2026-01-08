# Decision Theory

This submodule implements mathematical models of decision making.

## Components

- **ddm.py**: Drift-Diffusion Models for reaction times and accuracy.

## Usage

```python
from metainformant.math.decision_theory import ddm_analytic_accuracy

acc = ddm_analytic_accuracy(drift_rate=0.5, boundary=1.0)
```
