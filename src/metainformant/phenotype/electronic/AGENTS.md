# AI Agents in Electronic Phenotype Module

This document outlines AI assistance in developing the electronic phenotype analysis module (GPS, RFID, motion sensors).

## AI Contributions

The **Code Assistant Agent** developed:
- `TrackingPoint` dataclass for individual readings with coordinates and timestamp.
- `Trajectory` class for aggregating tracking points into movement paths.
- Functions for speed, distance, and activity budget calculations.

## Function Index

| Function/Class | Description |
|----------------|-------------|
| `TrackingPoint` | Dataclass for lat, lon, time, sensor_id |
| `Trajectory` | Aggregates points, calculates total distance, speed profile |
| `calculate_activity_budget()` | Time spent in rest vs. movement |
| `detect_stationary_periods()` | Identifies periods of low activity |

## Related Documentation

- **README**: [README.md](README.md)
- **SPEC**: [SPEC.md](SPEC.md)
- **Parent**: [src/metainformant/phenotype/AGENTS.md](../AGENTS.md)
