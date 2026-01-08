# Perception Moduel

This submodule implements mathematical models of sensory perception, psychophysics, and signal detection theory (SDT).

## Components

- **psychophysics.py**: Fundamental laws connecting physical stimuli to perceived sensation (Weber, Fechner, Stevens).
- **signal_detection.py**: Metrics for analyzing detection performance (d', c, beta).

## Usage

```python
from metainformant.math.perception import d_prime, stevens_power_law

# Calculate sensitivity d' from hit and false alarm rates
dp = d_prime(hit_rate=0.8, false_alarm_rate=0.2)

# Calculate perceived brightness using Stevens' Power Law
brightness = stevens_power_law(intensity=50, exponent=0.33)
```

## Dependencies
- **scipy.stats**: Required for Z-score calculations (norm.ppf) in SDT.
