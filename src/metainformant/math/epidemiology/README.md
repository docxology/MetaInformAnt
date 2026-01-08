# Epidemiology

This submodule implements mathematical models for disease spread and control.

## Components

- **models.py**: Core epidemiological models (SIR, SEIR, SIS) and reproduction numbers.

## Usage

```python
from metainformant.math.epidemiology import sir_step, basic_reproduction_number

# Calculate R0
R0 = basic_reproduction_number(transmission_rate=0.3, recovery_rate=0.1)
```
