# Decision Theory

Drift-diffusion models (DDM) for decision making and cognitive processes under uncertainty.

## Contents

| File | Purpose |
|------|---------|
| `ddm.py` | Drift-diffusion model: closed-form timing and accuracy, the exact first-passage-time density, trial likelihood and maximum-likelihood fitting |

## Key Functions

| Function | Description |
|----------|-------------|
| `ddm_mean_decision_time()` | Closed-form mean decision time for an unbiased start, E[T] = (a / 2v) tanh(v a / 2s²); uses `noise_sd` and treats drift of either sign by symmetry |
| `ddm_analytic_accuracy()` | Closed-form probability of the upper-boundary response, 1 / (1 + exp(−v a / s²)) |
| `ddm_first_passage_density()` | Density (or log-density) of a response through either boundary at a given time, from the Navarro and Fuss (2009) small-time and large-time series; scalar or array input |
| `ddm_log_likelihood()` | Exact log-likelihood of one trial: the log of that density at the observed boundary and time |
| `fit_ddm_parameters()` | Maximum-likelihood fit of drift, boundary separation, start point and non-decision time, with the noise fixed by `noise_sd` |
| `island_model_update()` | Island model allele frequency update with migration (no drift); frequencies and rate must lie in [0, 1] |

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
- `ddm_first_passage_density`, `ddm_log_likelihood` and `fit_ddm_parameters` use the exact first-passage density of dx = v dt + s dW between boundaries at 0 and `a`, started at `starting_point * a`, for both boundaries (choice 1 is the upper one, 2 the lower). It is zero for `rt <= non_decision_time`, integrates over time to the choice probability (`expit(v a / s^2)` for an unbiased start and the upper boundary, the same value as `ddm_analytic_accuracy`) and is consistent with the closed forms above. A response against the sign of the drift has a finite likelihood.
- The series are truncated for a relative error below 1e-12 (see `_series_term_counts` for the bound); the density is evaluated through its logarithm, so slow or fast responses far in the tails keep a finite log-likelihood. `noise_sd` (default 1) is the diffusion standard deviation; drift, boundary and noise are identified only up to a common factor, so `fit_ddm_parameters` fixes the noise and reports drift and boundary in its units.
- `fit_ddm_parameters` is deterministic (L-BFGS-B from the three best points of a fixed grid), accepts partial `bounds` (a missing parameter keeps its default; unknown names raise), lowers the upper bound of `non_decision_time` to just below the fastest response, and raises `ValueError` for invalid data and `RuntimeError` when it does not converge. It never returns default parameters.
- See [Drift–Diffusion Model](../../../../docs/math/ddm.md) for the closed forms.
