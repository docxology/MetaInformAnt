# Phenotype Analysis

Life course trajectory analysis linking event sequences to phenotypic traits, critical period identification, and outcome prediction.

## Contents

| File | Purpose |
|------|---------|
| `life_course.py` | Event-to-phenotype mapping, trajectory analysis, critical periods |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `Event` | Dataclass for a single life event with type, age, and metadata |
| `EventSequence` | Ordered collection of events for an individual |
| `extract_phenotypes_from_events()` | Derive phenotypic traits from event sequences |
| `aggregate_temporal_phenotypes()` | Aggregate traits across time windows |
| `analyze_life_course_trajectories()` | Cluster and compare life course patterns |
| `identify_critical_periods()` | Find age ranges with high phenotypic impact |
| `predict_life_course_outcomes()` | Predict future outcomes from event history |
| `create_life_course_report()` | Generate summary report of trajectory analysis |

## Usage

```python
from metainformant.phenotype.analysis.life_course import (
    EventSequence,
    analyze_life_course_trajectories,
    identify_critical_periods,
)

trajectories = analyze_life_course_trajectories(sequences)
critical = identify_critical_periods(sequences, age_ranges=[(0, 5), (5, 15)])
```
