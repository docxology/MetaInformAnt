# Decision Theory

Drift-diffusion models (DDM) for decision making and cognitive processes under uncertainty.

## Contents

| File | Purpose |
|------|---------|
| `ddm.py` | Drift-diffusion model fitting, timing, accuracy, and likelihood |

## Key Functions

| Function | Description |
|----------|-------------|
| `ddm_mean_decision_time()` | Mean decision time given drift rate and boundary separation |
| `ddm_analytic_accuracy()` | Analytic accuracy for drift-diffusion model |
| `ddm_log_likelihood()` | Log-likelihood of observed data under DDM parameters |
| `fit_ddm_parameters()` | Fit DDM parameters to observed reaction time data |
| `island_model_update()` | Island model allele frequency update with migration |

## Usage

```python
from metainformant.math.decision_theory.ddm import (
    ddm_mean_decision_time,
    ddm_analytic_accuracy,
    fit_ddm_parameters,
)

mean_rt = ddm_mean_decision_time(drift_rate=0.3, boundary=1.5)
accuracy = ddm_analytic_accuracy(drift_rate=0.3, boundary=1.5)
params = fit_ddm_parameters(reaction_times, choices)
```
