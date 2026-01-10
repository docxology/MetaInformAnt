# Sonic Phenotype Module

The `sonic` module handles acoustic and vibrational signaling data, including stridulation recordings and substrate-borne vibrations.

## Overview
This module handles:
- **Acoustic Signals**: Pulse rate, duration, carrier frequency.
- **Signal Processing**: Basic filtering and property extraction.
- **Classification**: Signal categorization.

## Components
- `AcousticSignal`: Class representing a raw or processed audio signal.
- `SignalAnalyst`: Utilities for signal processing.

## Usage
```python
from metainformant.phenotype.sonic import AcousticSignal

# Load a signal
signal = AcousticSignal.from_file("stridulation_001.wav")

# Analyze properties
print(f"Duration: {signal.duration}s")
print(f"Dominant Frequency: {signal.dominant_frequency}Hz")
```
