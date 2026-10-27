# Epidemiology

Mathematical models for infectious disease dynamics including SIR, SEIR, and SIS compartmental models.

## Contents

| File | Purpose |
|------|---------|
| `models.py` | Compartmental epidemic models and reproduction number calculations |

## Key Functions

| Function | Description |
|----------|-------------|
| `basic_reproduction_number()` | Calculate R0 from transmission and recovery rates |
| `effective_reproduction_number()` | Calculate Re accounting for susceptible fraction |
| `herd_immunity_threshold()` | Herd immunity threshold from R0 |
| `sir_step()` | Single time step of SIR model (Susceptible-Infected-Recovered) |
| `seir_step()` | Single time step of SEIR model (adds Exposed compartment) |
| `sis_step()` | Single time step of SIS model (no permanent immunity) |

## Usage

```python
from metainformant.math.epidemiology.models import (
    basic_reproduction_number,
    sir_step,
    herd_immunity_threshold,
)

r0 = basic_reproduction_number(transmission_rate=0.3, recovery_rate=0.1)
threshold = herd_immunity_threshold(r0)
s, i, r = sir_step(susceptible=990, infected=10, recovered=0, transmission_rate=0.3, recovery_rate=0.1, dt=0.1)
```

`sir_step` also accepts the keyword-only aliases `beta` and `gamma` for `transmission_rate` and `recovery_rate`; the compartment arguments are `susceptible`, `infected` and `recovered`, in that order.
