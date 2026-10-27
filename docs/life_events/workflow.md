# Life Events: Workflow Pipelines

High-level orchestration functions for end-to-end life event analysis.

## Concept Overview

The workflow module provides four top-level functions that combine data processing,
model training, statistical analysis, and visualization into complete pipelines:

- **analyze_life_course** -- Full analysis with embeddings, optional outcome prediction, and summary plots
- **compare_populations** -- Statistical comparison between two population groups
- **intervention_analysis** -- Before/after analysis around an intervention time point
- **event_importance** -- Rank event types by how early they occur or by their transition connectivity

Every function accepts event sequences and returns a results dictionary;
`analyze_life_course`, `compare_populations`, and `intervention_analysis` also write their outputs to a directory.

## Function Reference

### analyze_life_course

Complete life course analysis pipeline.

```python
from metainformant.life_events.workflow.workflow import analyze_life_course

result = analyze_life_course(
    sequences=sequences,
    outcomes=outcomes,
    output_dir="output/life_events/analysis",
    config_obj=config    # Optional LifeEventsWorkflowConfig or dictionary
)
```

**Pipeline:** convert the sequences to `"domain:event_type"` tokens, learn word2vec embeddings
(`learn_event_embeddings`), compute summary statistics, and, when `outcomes` are given, train an
`EventSequencePredictor` on those embeddings and predict the input sequences. Then write the summary plots and
`analysis_results.json`. `output_dir` defaults to `output/life_events/` (to the configuration's `work_dir` when a
`LifeEventsWorkflowConfig` is passed without `output_dir`). The function raises `ValueError` for an empty `sequences`
or an `outcomes` list of a different length.

**Keyword arguments** override the configuration: `embedding_dim` (100), `window_size` (5), `epochs` (5),
`random_state` (42), `model_type` (`"embedding"` or `"simple"`; default `"embedding"`), `task_type`
(`"classification"` or `"regression"`; default `"classification"`), `save_model` (`True`), and the embedding options
`method`, `learning_rate`, `min_count`, and `negative`, which are passed to `learn_event_embeddings` only when given.
`config_path` loads a configuration file with `load_life_events_config`.

**Returns:**

| Key | Content |
|-----|---------|
| `n_sequences`, `n_events` | Number of sequences and events analysed |
| `output_dir` | Directory the files were written to |
| `config_used` | The `embedding`, `model`, and `workflow` sections that were read (empty without a configuration) |
| `sequence_stats` | `get_event_statistics` output |
| `embeddings` | Path of `event_embeddings.json` |
| `embedding_dim`, `model_type`, `task_type` | Settings used |
| `predictions` | Predictions for the input sequences (empty list without `outcomes`) |
| `visualizations` | One `{"type", "file"}` entry per plot that was written |
| `visualization_errors` | One `{"type", "error_type", "error"}` entry per plot that failed; an empty list when every plot succeeded |
| `model`, `model_path` | Path of `trained_model.json`; present only with `outcomes` and `save_model` true |

**Plots:** `sequence_lengths.png`, `domain_distribution.png`, and, with `outcomes`, `outcome_distribution.png`.
Every plot is made independently: one that fails is recorded in `visualization_errors` and logged, and the other plots and the
analysis continue. `error_type` is the exception class name, or `"NoFigure"` when the plot function returned no figure
(matplotlib is missing or there is nothing to draw).

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

**Returns:** `{"group1_name", "group2_name", "group1_stats", "group2_stats", "group1", "group2", "comparison", "visualizations", "output_dir"}`.
`comparison` includes `length_comparison` (mean and median sequence length per group), `event_types` (unique counts,
overlap, and Jaccard similarity of the event-type sets), `common_event_types`, `unique_to_group1`, and
`unique_to_group2`. The function writes `comparison_results.json` and, when the comparison plot can be drawn,
`population_comparison.png`, and raises `ValueError` when either group is empty.

### intervention_analysis

Before/after analysis around an intervention time point.

```python
from datetime import datetime
from metainformant.life_events.workflow.workflow import intervention_analysis

result = intervention_analysis(
    sequences=all_sequences,
    intervention_time=datetime(2023, 6, 1),
    output_dir="output/life_events/intervention",
    pre_intervention_outcomes=pre_outcomes,     # Optional: one value per sequence with events before the intervention
    post_intervention_outcomes=post_outcomes,   # Optional: one value per sequence with events from the intervention on
)
```

Each sequence is split at `intervention_time` (a `datetime` or POSIX seconds; events at the intervention time count as
after it) and the events before and after it are analysed with `analyze_life_course` (files under `output_dir/pre`
and `output_dir/post`). A list of outcomes must have one value per sequence on its side (`n_pre_sequences` and
`n_post_sequences` of the result), otherwise `analyze_life_course` raises `ValueError`. The single `outcomes` argument is
split by position instead: its first `n_pre_sequences` values belong to the pre-intervention side and the following values to the
post-intervention side.

**Returns:** `{"pre_intervention", "post_intervention", "differences", "outcome_change", "intervention_time", "n_pre_sequences", "n_post_sequences", "output_dir"}`.
`pre_intervention` and `post_intervention` are `analyze_life_course` results (`{}` when no sequence has events on that
side), `differences` holds the post-minus-pre change of every numeric field they share (for example `n_events_change`), and
`outcome_change` holds `pre_mean` and `post_mean`, the mean outcome on each side, and `mean`, their difference (post minus pre); it is empty unless outcomes are given for both sides.

### event_importance

Rank event types with a descriptive score.

```python
from metainformant.life_events.workflow.workflow import event_importance

importance = event_importance(
    sequences=sequences,
    method="frequency",    # "frequency" or "transition"
    normalize=True
)
# Returns Dict[str, float] mapping event types to importance scores
```

**Methods:**
- `"frequency"` -- `1 - mean normalised position`: the position of every occurrence is scaled to [0, 1] within its sequence, so event types that tend to occur early score higher; occurrence counts do not enter the score
- `"transition"` -- number of transitions into and out of the event type between consecutive events

Any other `method` returns an empty dictionary. With `normalize=True` the scores are rescaled to [0, 1] (all scores equal gives 0.5). Keys are event types without the domain prefix. For importance based on a fitted predictor, use `event_importance(predictor, sequences, method="permutation")` and `shapley_values` in `metainformant.life_events.analysis.interpretability`.

## End-to-End Example

```python
from metainformant.life_events.core.utils import generate_synthetic_life_events
from metainformant.life_events.workflow.workflow import (
    analyze_life_course, compare_populations, event_importance,
)

# Synthetic cohort with a binary outcome
sequences, outcomes = generate_synthetic_life_events(
    n_sequences=120, generate_outcomes=True, outcome_relationship="event_based", random_state=0
)

# Full analysis
analysis = analyze_life_course(sequences, outcomes, "output/life_events/full", embedding_dim=32, epochs=5)
print(analysis["n_sequences"], analysis["visualization_errors"])

# Compare subgroups
short = [s for s in sequences if len(s.events) <= 10]
long = [s for s in sequences if len(s.events) > 10]
compare_populations(short, long, ["Short", "Long"], "output/life_events/compare")

# Event ranking
rankings = event_importance(sequences, method="transition", normalize=True)
for event_type, score in sorted(rankings.items(), key=lambda x: -x[1])[:10]:
    print(f"  {event_type}: {score:.3f}")
```

## Configuration

`analyze_life_course` reads the nested `embedding`, `model`, and `workflow` sections of a configuration. The configuration is a
`LifeEventsWorkflowConfig` (`config_obj`), a plain dictionary with the same sections, or a file (`config_path`). A keyword argument
of the same name takes precedence:

```yaml
work_dir: output/life_events/work     # default output_dir when a LifeEventsWorkflowConfig is passed without one
embedding:
  embedding_dim: 100                  # default 100
  window_size: 5                      # default 5
  epochs: 10                          # default 5
  method: skipgram                    # "skipgram" or "cbow"
  learning_rate: 0.01                 # initial word2vec learning rate (learn_event_embeddings default 0.025)
  min_count: 1                        # optional, passed to learn_event_embeddings
  negative: 5                         # optional, passed to learn_event_embeddings
model:
  model_type: embedding               # "embedding" or "simple"
  task_type: classification           # or "regression"
  random_state: 42                    # seeds the embedding training and the predictor
workflow:
  save_model: true                    # write trained_model.json when outcomes are given
```

The complete template is `config/life_events/life_events_template.yaml`.

- With a `LifeEventsWorkflowConfig`, only these nested sections are read: top-level scalar keys such as `embedding_dim` or `model_type` set attributes of the configuration object but do not change what `analyze_life_course` does. A plain dictionary may also carry `embedding_dim`, `window_size`, `epochs`, `learning_rate`, `model_type`, `task_type`, and `random_state` at the top level.
- `load_life_events_config(path, prefix="LE")` requires `work_dir`. The environment variables `LE_THREADS`, `LE_WORK_DIR`, and `LE_LOG_DIR` override the loaded `threads`, `work_dir`, and `log_dir`; no other key has an environment override.

```python
from metainformant.life_events.core.config import load_life_events_config
from metainformant.life_events.workflow.workflow import analyze_life_course

config = load_life_events_config("config/life_events/life_events_template.yaml")
result = analyze_life_course(sequences, outcomes, config_obj=config)
```

## See Also

- [Events](events.md) -- Input data structures
- [Models](models.md) -- Underlying prediction models
- [Survival](survival.md) -- Time-to-event analysis
- [Visualization](visualization.md) -- All available plot types
