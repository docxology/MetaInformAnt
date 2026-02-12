# Life Events: Workflow Pipelines

High-level orchestration functions for end-to-end life event analysis.

## Concept Overview

The workflow module provides four top-level functions that combine data processing,
model training, statistical analysis, and visualization into complete pipelines:

- **analyze_life_course** -- Full analysis with embeddings, prediction, and plots
- **compare_populations** -- Statistical comparison between two population groups
- **intervention_analysis** -- Before/after analysis around an intervention time point
- **event_importance** -- Rank events by frequency, predictive power, or temporal position

Each function accepts event sequences, runs a multi-step pipeline, writes outputs
to a specified directory, and returns a results dictionary.

## Function Reference

### analyze_life_course

Complete life course analysis pipeline.

```python
from metainformant.life_events.workflow.workflow import analyze_life_course

result = analyze_life_course(
    sequences=sequences,
    outcomes=outcomes,
    output_dir="output/life_events/analysis",
    config_obj=config    # Optional LifeEventsWorkflowConfig
)
```

**Pipeline:** Validate sequences, learn embeddings, train predictor, generate
visualizations, compute summary statistics.

**Returns:** `{"embeddings", "model", "predictions", "statistics", "plots"}`

### compare_populations

Statistical comparison between two groups.

```python
from metainformant.life_events.workflow.workflow import compare_populations

result = compare_populations(
    sequences1=treatment_group,
    sequences2=control_group,
    group_names=["Treatment", "Control"],
    output_dir="output/life_events/comparison"
)
```

**Returns:** `{"group1_stats", "group2_stats", "comparison", "plots"}`

### intervention_analysis

Before/after analysis around an intervention time point.

```python
from metainformant.life_events.workflow.workflow import intervention_analysis

result = intervention_analysis(
    sequences=all_sequences,
    intervention_time=datetime(2023, 6, 1),
    output_dir="output/life_events/intervention",
    outcomes=outcome_values    # Optional
)
```

**Returns:** `{"pre_stats", "post_stats", "comparison", "effect_size", "plots"}`

### event_importance

Rank events by importance using different methods.

```python
from metainformant.life_events.workflow.workflow import event_importance

importance = event_importance(
    sequences=sequences,
    method="frequency",    # "frequency", "predictive", or "temporal"
    normalize=True
)
# Returns Dict[str, float] mapping event types to importance scores
```

**Methods:**
- `"frequency"` -- Rank by occurrence count
- `"predictive"` -- Rank by contribution to prediction accuracy
- `"temporal"` -- Rank by temporal position (early events weighted higher)

## End-to-End Example

```python
from metainformant.life_events.core.events import EventDatabase
from metainformant.life_events.core.config import load_life_events_config
from metainformant.life_events.workflow.workflow import (
    analyze_life_course, compare_populations, event_importance,
)

# Load data
db = EventDatabase.load_json("data/cohort_events.json")
sequences = db.sequences
outcomes = [1 if "graduation" in s.get_event_types() else 0 for s in sequences]

# Full analysis
config = load_life_events_config("config/life_events.yaml")
analysis = analyze_life_course(sequences, outcomes, "output/life_events/full", config)

# Compare subgroups
urban = db.filter_by_domain("urban").sequences
rural = db.filter_by_domain("rural").sequences
compare_populations(urban, rural, ["Urban", "Rural"], "output/life_events/compare")

# Event ranking
rankings = event_importance(sequences, method="frequency", normalize=True)
for event_type, score in sorted(rankings.items(), key=lambda x: -x[1])[:10]:
    print(f"  {event_type}: {score:.3f}")
```

## Configuration

Workflow configuration via `LifeEventsWorkflowConfig`:

```yaml
work_dir: output/life_events
embedding_dim: 100
model_type: embedding
learning_rate: 0.001
batch_size: 32
epochs: 100
random_seed: 42
```

Environment overrides use the `LE_` prefix (e.g., `LE_EMBEDDING_DIM=64`).

```python
from metainformant.life_events.core.config import load_life_events_config
config = load_life_events_config("config/life_events.yaml", prefix="LE")
```

## See Also

- [Events](events.md) -- Input data structures
- [Models](models.md) -- Underlying prediction models
- [Survival](survival.md) -- Time-to-event analysis
- [Visualization](visualization.md) -- All available plot types
