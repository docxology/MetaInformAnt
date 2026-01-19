# AI Agents in Morphological Phenotype Module

This document outlines AI assistance in developing the morphological phenotype analysis module.

## AI Contributions

The **Code Assistant Agent** developed:
- `Measurement` dataclass for individual landmarks with value and unit.
- `MorphometricProfile` class for aggregating body measurements.
- Unit conversion and geometric morphometrics functions.

## Function Index

| Function/Class | Description |
|----------------|-------------|
| `Measurement` | Dataclass for name, value, unit |
| `MorphometricProfile` | Aggregates measurements, provides summary statistics |
| `convert_units()` | Converts between mm, cm, inches |
| `calculate_shape_ratios()` | Computes standard morphometric indices |

## Related Documentation

- **README**: [README.md](README.md)
- **SPEC**: [SPEC.md](SPEC.md)
- **Parent**: [src/metainformant/phenotype/AGENTS.md](../AGENTS.md)
