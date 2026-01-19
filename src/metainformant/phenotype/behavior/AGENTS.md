# AI Agents in Behavior Module

This document outlines AI assistance in developing the behavior analysis (ethology) module.

## AI Contributions

The **Code Assistant Agent** developed:
- `Ethogram` class for defining behavioral categories and codes.
- `BehaviorSequence` class for temporal sequences of recorded behaviors.
- Functions for transition matrix calculation and bout duration analysis.

## Function Index

| Function/Class | Description |
|----------------|-------------|
| `Ethogram` | Dictionary of behavior codes and descriptions |
| `BehaviorSequence` | Time-stamped list of behaviors per subject |
| `transition_matrix()` | Calculates probability of behavior transitions |
| `calculate_bout_durations()` | Statistics on behavior duration |
| `time_budget()` | Proportion of time spent in each behavior |

## Related Documentation

- **README**: [README.md](README.md)
- **SPEC**: [SPEC.md](SPEC.md)
- **Parent**: [src/metainformant/phenotype/AGENTS.md](../AGENTS.md)
