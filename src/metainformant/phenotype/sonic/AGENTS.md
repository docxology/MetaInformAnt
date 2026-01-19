# AI Agents in Sonic Phenotype Module

This document outlines AI assistance in developing the sonic phenotype (bioacoustics) analysis module.

## AI Contributions

The **Code Assistant Agent** developed:
- `AcousticSignal` dataclass for waveform and spectral features.
- Functions for frequency spectrum analysis and call detection.
- Integration with scipy for FFT and spectrogram generation.

## Function Index

| Function/Class | Description |
|----------------|-------------|
| `AcousticSignal` | Dataclass for waveform, sample_rate, duration, species_id |
| `extract_spectrogram()` | Generates spectrogram from waveform |
| `detect_calls()` | Identifies acoustic events above threshold |
| `calculate_peak_frequency()` | Returns dominant frequency |

## Related Documentation

- **README**: [README.md](README.md)
- **SPEC**: [SPEC.md](SPEC.md)
- **Parent**: [src/metainformant/phenotype/AGENTS.md](../AGENTS.md)
