# Life Events: Survival Analysis

Pure-Python time-to-event analysis for life event data.

## Concept Overview

Survival analysis estimates the probability of an event occurring over time,
accounting for censored observations (subjects who leave the study or have not
yet experienced the event). This module provides five core methods:

- **Kaplan-Meier** -- Non-parametric survival curve estimation
- **Cox Proportional Hazards** -- Semi-parametric regression with covariates
- **Competing Risks** -- Multiple mutually exclusive event types
- **Recurrent Events** -- Repeated occurrences per subject
- **Time-Varying Covariates** -- Covariates that change during follow-up

All implementations are pure NumPy with no external survival analysis dependencies.

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
| `time_points` | `np.ndarray` | Unique event times |
| `survival_probability` | `np.ndarray` | S(t) at each time point |
| `confidence_interval` | `Dict` | Upper and lower 95% CI bounds |
| `n_at_risk` | `np.ndarray` | Number at risk at each time |
| `median_survival` | `float` | Time at S(t) = 0.5 |
| `n_events` | `int` | Total observed events |
| `n_censored` | `int` | Total censored observations |
| `groups` | `Dict` | Per-group results (if groups provided) |

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
```

**Returns:**

| Key | Type | Description |
|-----|------|-------------|
| `coefficients` | `Dict[str, float]` | Beta coefficients per covariate |
| `hazard_ratios` | `Dict[str, float]` | exp(beta) per covariate |
| `se` | `Dict[str, float]` | Standard errors |
| `p_values` | `Dict[str, float]` | Wald test p-values |
| `concordance` | `float` | C-index (discrimination) |

### competing_risks

Cumulative incidence estimation when multiple event types compete.

```python
from metainformant.life_events.survival.time_to_event import competing_risks

result = competing_risks(
    times=[2.0, 5.0, 8.0, 10.0, 12.0],
    events=[1, 1, 1, 0, 1],
    event_types=["relapse", "death", "relapse", "censored", "death"]
)
```

**Returns:**

| Key | Type | Description |
|-----|------|-------------|
| `cumulative_incidence` | `Dict[str, np.ndarray]` | CIF per event type |
| `cause_specific_hazards` | `Dict[str, np.ndarray]` | Hazard per cause |
| `overall_survival` | `np.ndarray` | Combined survival curve |

### recurrent_events

Analysis of events that can happen multiple times per subject.

```python
from metainformant.life_events.survival.time_to_event import recurrent_events

result = recurrent_events(
    times=[1.0, 3.0, 5.0, 2.0, 6.0],
    events=[1, 1, 1, 1, 0],
    subject_ids=["S1", "S1", "S1", "S2", "S2"]
)
```

**Returns:**

| Key | Type | Description |
|-----|------|-------------|
| `mean_cumulative_function` | `np.ndarray` | Expected cumulative events over time |
| `rate_per_subject` | `Dict[str, float]` | Event rate per subject |
| `gap_time_distribution` | `np.ndarray` | Distribution of inter-event times |

### time_varying_covariates

Expand data for covariates that change during follow-up.

```python
from metainformant.life_events.survival.time_to_event import time_varying_covariates

intervals = [
    {"subject": "S1", "start": 0, "stop": 5, "event": 0, "treatment": 0},
    {"subject": "S1", "start": 5, "stop": 10, "event": 1, "treatment": 1},
    {"subject": "S2", "start": 0, "stop": 8, "event": 0, "treatment": 1},
]

result = time_varying_covariates(intervals)
```

**Returns:**

| Key | Type | Description |
|-----|------|-------------|
| `expanded_data` | `List[Dict]` | Interval-format records |
| `interval_counts` | `Dict` | Intervals per subject |
| `n_subjects` | `int` | Unique subjects |
| `n_intervals` | `int` | Total intervals |
| `covariate_names` | `List[str]` | Detected covariate names |

## End-to-End Example

```python
from metainformant.life_events.core.utils import generate_cohort_sequences
from metainformant.life_events.survival.time_to_event import (
    kaplan_meier_estimator,
    cox_ph_model,
)

# Generate synthetic cohort
sequences = generate_cohort_sequences(n_sequences=200, avg_events_per_sequence=10)

# Extract time-to-event data
times = [float(len(s.events)) for s in sequences]
events = [1 if len(s.events) > 8 else 0 for s in sequences]
groups = ["high" if len(s.events) > 12 else "low" for s in sequences]

# Kaplan-Meier by group
km = kaplan_meier_estimator(times, events, groups)
print(f"Median survival: {km['median_survival']:.1f}")
for group, data in km["groups"].items():
    print(f"  {group}: median = {data['median_survival']:.1f}")

# Cox model with covariates
covariates = [[len(s.events), len(s.get_domains())] for s in sequences]
cox = cox_ph_model(times, events, covariates, ["n_events", "n_domains"])
for name, hr in cox["hazard_ratios"].items():
    print(f"  {name}: HR={hr:.3f}, p={cox['p_values'][name]:.4f}")
```

## Configuration

Survival analysis parameters can be set via workflow config:

```yaml
work_dir: output/life_events
random_seed: 42
```

Environment prefix: `LE_` (e.g., `LE_RANDOM_SEED=123`).

## See Also

- [Events](events.md) -- Input event data structures
- [Models](models.md) -- SurvivalPredictor for ML-based survival
- [Visualization](visualization.md) -- Plotting survival curves
- [Workflow](workflow.md) -- Integrated analysis pipelines
