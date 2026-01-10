# Behavior Module

The `behavior` module provides tools for analyzing behavioral data, including ethograms, behavioral sequences, and time-budget analysis.

## Overview
This module handles:
- **Ethograms**: Definitions of behavioral repertoires.
- **Behavior Sequences**: Time-stamped sequences of behaviors.
- **State Analysis**: Markov chains and transition matrices for behavioral states.

## Components
- `Ethogram`: Class for defining a set of valid behaviors and their descriptions.
- `BehaviorSequence`: Class for handling time-series behavioral data.

## Usage
```python
from metainformant.phenotype.behavior import Ethogram, BehaviorSequence

# Define an ethogram
ethogram = Ethogram({
    "rest": "Stationary with no movement",
    "walk": "Locomotion",
    "forage": "Searching for food"
})

# Analyze a sequence
seq = BehaviorSequence(data=..., ethogram=ethogram)
print(seq.calculate_time_budget())
```
