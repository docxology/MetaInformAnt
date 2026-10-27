# Life Events: Survival Analysis

Time-to-event analysis for life event data.

## Concept Overview

Survival analysis estimates the probability of an event occurring over time,
accounting for censored observations (subjects who leave the study or have not
yet experienced the event). This module provides five core methods:

- **Kaplan-Meier** -- Non-parametric survival curve estimation
- **Cox Proportional Hazards** -- Semi-parametric regression with covariates
- **Competing Risks** -- Multiple mutually exclusive event types
- **Recurrent Events** -- Repeated occurrences per subject
- **Time-Varying Covariates** -- Covariates that change during follow-up

The Kaplan-Meier, competing-risks, recurrent-event, and time-varying-covariate routines are pure Python
step-function estimators; the Cox fit (`cox_ph_model`) runs on NumPy. No external survival analysis package is required.
For life event sequences, `SurvivalPredictor` (see [Models](models.md)) fits Cox, exponential, or Kaplan-Meier models on
event-count covariates.

## Function Reference

### kaplan_meier_estimator

Non-parametric survival probability estimation with optional group comparison.

```python
from metainformant.life_events.survival.time_to_event import kaplan_meier_estimator

result = kaplan_meier_estimator(
    times=[5.0, 8.0, 12.0, 15.0, 20.0, 25.0],
    events=[1, 0, 1, 1, 0, 1],         # 1=event, 0=censored
    groups=["A", "A", "A", "B", "B", "B"]  # Optional stratification
)
```

**Returns:**

| Key | Type | Description |
|-----|------|-------------|
| `time_points` | `List[float]` | Sorted unique times at which an event occurred (times with only censored subjects are not listed) |
| `survival_probability` | `List[float]` | S(t) at each time point |
| `confidence_interval` | `List[Tuple[float, float]]` | (lower, upper) 95% bounds per time point: S(t) ± 1.96 × the Greenwood standard error, clipped to [0, 1] |
| `n_at_risk` | `List[int]` | Number at risk at each time point |
| `median_survival` | `float` or `None` | First event time with S(t) ≤ 0.5; `None` when it is not reached |
| `n_events` | `int` | Total observed events |
| `n_censored` | `int` | Total censored observations |
| `groups` | `Dict` | Per-group results, each a complete result of this table; present only when `groups` is given |

With `groups`, the top-level `survival_probability`, `confidence_interval`, and `n_at_risk` are empty lists,
`median_survival` is `None`, `time_points` is the union of the groups' event times, and the counts are totals; the curves
are in `result["groups"][label]`. `ValueError` is raised when `times` and `events` differ in length.

### cox_ph_model

Semi-parametric regression modeling the effect of covariates on hazard rates.

```python
from metainformant.life_events.survival.time_to_event import cox_ph_model

result = cox_ph_model(
    times=[3.0, 7.0, 10.0, 14.0, 18.0],
    events=[1, 1, 0, 1, 0],
    covariates=[[1.2, 0], [0.8, 1], [1.5, 0], [0.9, 1], [1.1, 0]],
    covariate_names=["risk_score", "treatment"]
)
for name, hr, se, p in zip(result["covariate_names"], result["hazard_ratios"], result["se"], result["p_values"]):
    print(f"{name}: HR={hr:.2f}, SE(log HR)={se:.2f}, p={p:.3f}")
```

The coefficients maximise the Breslow partial likelihood by Newton-Raphson with the exact gradient and information matrix
and step halving (at most 50 iterations); tied event times share one risk set, and right-censoring is supported. The standard
errors are the square roots of the diagonal of the inverse information matrix, so correlated covariates are accounted for, and the
p-values are two-sided Wald tests. The same fit with `penalizer=0` is available on event-count covariates as
`SurvivalPredictor(method="cox", penalizer=0)`.

**Returns:** a dictionary whose list values are plain Python lists, one entry per covariate in the order given.

| Key | Type | Description |
|-----|------|-------------|
| `coefficients` | `List[float]` | Beta coefficients (log hazard ratios) |
| `hazard_ratios` | `List[float]` | exp(beta) per covariate |
| `se` | `List[float]` | Standard errors of the coefficients; `inf` for a covariate that is not identified |
| `p_values` | `List[float]` | Wald test p-values; 1 for a covariate that is not identified |
| `concordance` | `float` | Harrell's C-index (discrimination); pairs with equal risk scores are not counted, and 0.5 is returned when no pair is comparable |
| `covariate_names` | `List[str]` | `covariate_names` as given, or `X0`, `X1`, ... |
| `n_subjects` | `int` | Number of subjects |
| `n_events` | `int` | Number of events |
| `log_likelihood` | `float` | Maximised partial log-likelihood |
| `converged` | `bool` | Whether Newton-Raphson converged |
| `n_iterations` | `int` | Number of Newton-Raphson iterations |

A covariate that is not identified (constant, or a linear combination of other covariates) gets an infinite standard error and
a p-value of 1. A covariate that perfectly orders the event times makes the likelihood monotone: its coefficient grows until the
likelihood stops improving and its standard error is very large. A fit that does not converge within 50 iterations emits a
`RuntimeWarning`. `ValueError` is raised when `times`, `events`, and `covariates` differ in length, there is no covariate, a value is
not a finite number, or `events` holds anything other than 0 and 1 (booleans are accepted).

### competing_risks

Cumulative incidence estimation when multiple event types compete.

```python
from metainformant.life_events.survival.time_to_event import competing_risks

result = competing_risks(
    times=[2.0, 5.0, 8.0, 10.0, 12.0],
    events=[1, 1, 1, 0, 1],
    event_types=[1, 2, 1, 0, 2]     # event type of each subject with an event (1 = relapse, 2 = death)
)
```

**Returns:**

| Key | Type | Description |
|-----|------|-------------|
| `cumulative_incidence_per_type` | `Dict[type, List[Tuple[float, float]]]` | `(time, cumulative incidence)` pairs per event type |
| `cause_specific_hazards` | `Dict[type, List[Tuple[float, float]]]` | `(time, hazard)` pairs per event type, at the times that event type occurs |
| `overall_survival` | `List[Tuple[float, float]]` | `(time, survival)` pairs of the combined survival curve |

The event type of a censored subject (`events` value 0) is ignored.

### recurrent_events

Analysis of events that can happen multiple times per subject. `times` and `events` hold one list per subject.

```python
from metainformant.life_events.survival.time_to_event import recurrent_events

result = recurrent_events(
    times=[[1.0, 3.0, 5.0], [2.0, 6.0]],      # event and censoring times of S1 and S2
    events=[[1, 1, 1], [1, 0]],               # 1 = event, 0 = censoring time
    subject_ids=["S1", "S2"]
)
```

**Returns:**

| Key | Type | Description |
|-----|------|-------------|
| `mean_cumulative_function` | `List[Tuple[float, float]]` | `(time, expected cumulative events)` at each event time |
| `rate_per_subject` | `Dict[str, float]` | Each subject's number of events divided by their longest recorded time |
| `gap_time_distribution` | `Dict[str, float]` | `mean`, `median`, `std`, `min`, `max`, and `n_gaps` of the times between consecutive events; empty when no subject has two events |
| `n_subjects` | `int` | Number of subjects |
| `total_events` | `int` | Events across all subjects |

### time_varying_covariates

Expand data for covariates that change during follow-up.

```python
from metainformant.life_events.survival.time_to_event import time_varying_covariates

intervals = [
    {"subject_id": "S1", "start": 0, "stop": 5, "event": 0, "covariates": {"treatment": 0}},
    {"subject_id": "S1", "start": 5, "stop": 10, "event": 1, "covariates": {"treatment": 1}},
    {"subject_id": "S2", "start": 0, "stop": 8, "event": 0, "covariates": {"treatment": 1}},
]

result = time_varying_covariates(intervals)
```

**Returns:**

| Key | Type | Description |
|-----|------|-------------|
| `expanded_data` | `List[Dict]` | Interval records sorted by subject and start, with `subject_id`, `start`, `stop`, `event`, and one key per covariate |
| `interval_counts` | `Dict` | Intervals per subject |
| `n_subjects` | `int` | Unique subjects |
| `n_intervals` | `int` | Total intervals |
| `covariate_names` | `List[str]` | Sorted covariate names |

## End-to-End Example

```python
import numpy as np

from metainformant.life_events.core.utils import generate_cohort_sequences
from metainformant.life_events.survival.time_to_event import cox_ph_model, kaplan_meier_estimator

sequences = generate_cohort_sequences(n_sequences=200, avg_events_per_sequence=10, random_state=1)

# Simulated follow-up: the hazard rises with the number of distinct life domains a person has events in
rng = np.random.default_rng(1)
n_events = np.array([len(s.events) for s in sequences])
n_domains = np.array([len(s.get_domains()) for s in sequences])
latent = rng.exponential(scale=20.0 / np.exp(0.4 * n_domains))
times = np.minimum(latent, 15.0)                  # administrative censoring at t = 15
events = (latent <= 15.0).astype(int)             # 1 = event observed, 0 = censored
groups = ["many_domains" if d >= 4 else "few_domains" for d in n_domains]

# Kaplan-Meier by group
km = kaplan_meier_estimator(times.tolist(), events.tolist(), groups)
for group, result in km["groups"].items():
    print(f"  {group}: events={result['n_events']}, median = {result['median_survival']}")

# Cox model with covariates
covariates = np.column_stack([n_events, n_domains]).tolist()
cox = cox_ph_model(times.tolist(), events.tolist(), covariates, ["n_events", "n_domains"])
for name, hr, p in zip(cox["covariate_names"], cox["hazard_ratios"], cox["p_values"]):
    print(f"  {name}: HR={hr:.3f}, p={p:.4f}")
```

## Configuration

The functions above take no configuration object. `SurvivalPredictor` (see [Models](models.md)) exposes `penalizer`, `max_iter`,
and `tol` for its Cox and exponential fits.

## See Also

- [Events](events.md) -- Input event data structures
- [Models](models.md) -- SurvivalPredictor for ML-based survival
- [Visualization](visualization.md) -- Plotting survival curves
- [Workflow](workflow.md) -- Integrated analysis pipelines
