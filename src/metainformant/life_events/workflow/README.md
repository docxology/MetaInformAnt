# Life Events Workflow

End-to-end life course analysis workflows: sequence loading, statistical analysis, population comparison, intervention modeling, and importance ranking.

## Contents

| File | Purpose |
|------|---------|
| `workflow.py` | Complete analysis pipelines for life event sequence data |

## Key Functions

| Function | Description |
|----------|-------------|
| `analyze_life_course()` | Full analysis of event sequences: statistics, patterns, visualizations |
| `compare_populations()` | Statistical comparison of event patterns between two population groups |
| `intervention_analysis()` | Model the effect of interventions on event sequence outcomes |
| `event_importance()` | Rank event types by their impact on target outcomes |

## Workflow Steps

1. Load event sequences from JSON or generate synthetic cohorts
2. Run `analyze_life_course()` for complete statistical profiling
3. Optionally compare populations via `compare_populations()`
4. Model interventions with `intervention_analysis()`
5. Identify key drivers with `event_importance()`

## Usage

```python
from metainformant.life_events.workflow.workflow import (
    analyze_life_course, compare_populations, intervention_analysis
)

results = analyze_life_course(sequences, output_dir="output/life_events/")
comparison = compare_populations(group_a, group_b)
intervention = intervention_analysis(sequences, intervention_type="early_treatment")
```
