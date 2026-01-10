# Chemical Phenotype Module

The `chemical` module handles chemotypical data, such as cuticular hydrocarbon (CHC) profiles, pheromone compositions, and other metabolomic data.

## Overview
This module handles:
- **Chemical Profiles**: Quantitative data on chemical compounds.
- **Compound Libraries**: Definitions of known compounds (retention times, mass spectra).
- **Similarity Analysis**: Comparing chemical profiles between samples.

## Components
- `Compound`: Data class for individual chemical compounds.
- `ChemicalProfile`: Class representing a sample's chemical composition.

## Usage
```python
from metainformant.phenotype.chemical import Compound, ChemicalProfile

# Define a compound
c25 = Compound(name="n-C25", formula="C25H52", retention_time=25.4)

# Create a profile
profile = ChemicalProfile(sample_id="ant_001", compounds={c25: 0.15, ...})
```
