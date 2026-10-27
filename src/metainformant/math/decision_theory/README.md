# Decision Theory

Drift-diffusion models (DDM) for decision making and cognitive processes under uncertainty.

## Contents

| File | Purpose |
|------|---------|
| `ddm.py` | Drift-diffusion model: closed-form timing and accuracy, an approximate trial likelihood and fitting |

## Key Functions

| Function | Description |
|----------|-------------|
| `ddm_mean_decision_time()` | Closed-form mean decision time for an unbiased start, E[T] = (a / 2v) tanh(v a / 2s²); uses `noise_sd` and treats drift of either sign by symmetry |
| `ddm_analytic_accuracy()` | Closed-form probability of the upper-boundary response, 1 / (1 + exp(−v a / s²)) |
| `ddm_log_likelihood()` | Log-likelihood of one trial under a normal approximation of the first-passage time (not the exact density) |
| `fit_ddm_parameters()` | Fit drift, boundary separation, start point and non-decision time to that approximate likelihood |
| `island_model_update()` | Island model allele frequency update with migration |

## Usage

```python
from metainformant.math.decision_theory.ddm import (
    ddm_mean_decision_time,
    ddm_analytic_accuracy,
    fit_ddm_parameters,
)

mean_rt = ddm_mean_decision_time(drift_rate=0.3, boundary=1.5)  # 0.553
accuracy = ddm_analytic_accuracy(drift_rate=0.3, boundary=1.5, noise_sd=1.0)  # 0.611
params = fit_ddm_parameters(reaction_times, choices)
```

## Notes

- The two closed forms (Palmer, Huk and Shadlen 2005; Bogacz et al. 2006) describe dx = v dt + s dW between absorbing boundaries a apart, started midway. They exclude non-decision time and raise `ValueError` unless `boundary` and `noise_sd` are positive and finite and `drift_rate` is finite.
- `ddm_log_likelihood` and `fit_ddm_parameters` use a normal approximation of the first-passage time: no noise parameter, the opposite boundary is ignored, and a response against the sign of the drift has a log-likelihood of −inf. Their time scale therefore differs from the closed forms, and fitted values are approximate.
- See [Drift–Diffusion Model](../../../../docs/math/ddm.md) for the formulas and the details of the approximation.
