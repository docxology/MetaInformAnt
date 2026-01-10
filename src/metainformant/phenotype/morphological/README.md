# Morphological Phenotype Module

The `morphological` module handles morphometric measurements, shape analysis, and allometric Scaling integration.

## Overview
This module handles:
- **Measurements**: Linear dimensions (e.g., Weber's length, head width).
- **Indices**: Ratios and calculated indices (e.g., CI, SI).
- **Geometric Morphometrics**: Landmark-based shape data (future).

## Components
- `Measurement`: Typed class for individual measurements with units.
- `MorphometricProfile`: Collection of measurements for a specimen.

## Usage
```python
from metainformant.phenotype.morphological import Measurement, MorphometricProfile

# Define measurements
hw = Measurement(value=1.5, unit="mm", name="Head Width")
hl = Measurement(value=1.6, unit="mm", name="Head Length")

# Create profile
profile = MorphometricProfile(specimen_id="spec_001", measurements=[hw, hl])
print(profile.calculate_index("CI", "Head Width", "Head Length"))
```
