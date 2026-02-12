# Survival Analysis

Time-to-event analysis methods including Kaplan-Meier estimation, Cox proportional hazards modelling, competing risks, recurrent events, and time-varying covariates.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `time_to_event` module |
| `time_to_event.py` | KM estimator, Cox PH, competing risks, recurrent events, time-varying covariates |

## Key Functions

| Function | Description |
|----------|-------------|
| `kaplan_meier_estimator()` | Kaplan-Meier survival curve estimation with Greenwood variance |
| `cox_ph_model()` | Cox proportional hazards model via Newton-Raphson optimization |
| `competing_risks()` | Cumulative incidence estimation for competing risk events |
| `recurrent_events()` | Analyze recurrent event processes |
| `time_varying_covariates()` | Prepare time-varying covariate data for survival models |

## Usage

```python
from metainformant.life_events.survival import time_to_event

km = time_to_event.kaplan_meier_estimator(times, events)
cox = time_to_event.cox_ph_model(times, events, covariates)
ci = time_to_event.competing_risks(times, events, event_types)
recur = time_to_event.recurrent_events(times, events, subject_ids)
```
