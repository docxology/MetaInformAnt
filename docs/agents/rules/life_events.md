# Life Events Module Rules (`metainformant.life_events`)

## Purpose
Life course event analysis and prediction: event extraction, embedding generation, sequence modeling, and interpretability.

## Dependencies
- **Required**: `core`
- **Optional**: PyTorch (the recurrent models in `sequence_models`)

## Source Structure
```
src/metainformant/life_events/
 analysis/
 interpretability.py
 core/
 config.py, events.py, utils.py
 models/
 embeddings.py, predictor.py, sequence_models.py, statistical_models.py
 survival/
 time_to_event.py
 visualization/
 network.py, statistical.py, timeline.py
 workflow/
 workflow.py
```

## Package Management
- **ALWAYS use `uv`** for all Python package management and environment operations
- Use `uv run` to execute commands: `uv run pytest`, `uv run metainformant life-events --help`
- **NEVER use `pip` directly**

## Key Submodules

### Events (`events`)
- Event extraction and analysis
- Event sequence processing
- Temporal patterns

### Embeddings (`embeddings`)
- Event embedding generation (`learn_event_embeddings`): word2vec with negative sampling, trained in NumPy, so events that occur close together in a sequence get similar vectors
- `method` is `"skipgram"` or `"cbow"` (defaults from `sg`); `window_size` is the context window; `embedding_dim`, `window_size`, and `epochs` must be positive
- Options: `min_count` (rarer events are dropped and not returned), `negative` (negative samples per positive example, default 5), `learning_rate` (initial rate, default 0.025, decaying linearly), `random_state` (identical seeds give identical embeddings)
- `combine="sum"` (default) returns each event's input plus output vector; `combine="input"` returns the input vectors only
- Biologically-informed variants: `biological_embedding`, `domain_specific_embeddings`
- **Actual path**: `metainformant.life_events.models.embeddings`

**Patterns**:
```python
# Import from actual path
from metainformant.life_events.models import embeddings

# event_sequences: a list of EventSequence objects (or lists of event-type tokens)
embeddings_dict = embeddings.learn_event_embeddings(
    sequences=event_sequences,
    embedding_dim=100,
    window_size=5,
    method="skipgram",
)
# Returns: {event_token: numpy vector of length embedding_dim}
```

### Models (`models`)
- `EventSequencePredictor` (`predictor`): outcome prediction, classification or regression, with `model_type` `"embedding"` or `"simple"`; any other `model_type`, including `"lstm"`, makes `fit` raise `ValueError`
- `MultiTaskPredictor` and `SurvivalPredictor` (`statistical_models`): NumPy models (ridge and logistic task heads; Cox, exponential, or Kaplan-Meier survival)
- `LSTMSequenceModel` and `GRUSequenceModel` (`sequence_models`): next-event models that require PyTorch and do not predict outcomes

### Survival (`survival`)
- `kaplan_meier_estimator`, `cox_ph_model`, `competing_risks`, `recurrent_events`, `time_varying_covariates`
- `cox_ph_model` is a NumPy Newton-Raphson fit of the Breslow partial likelihood; its list-valued results have `inf` standard errors for unidentified covariates, and invalid input raises `ValueError`

### Interpretability (`interpretability`)
- Model interpretation utilities
- Feature importance
- Event contribution analysis: `shapley_values` is exact up to `max_exact_events` events (default 10) and seeded permutation sampling above; the values sum to `v(all events) - v(empty sequence)`

### Visualization (`visualization`)
- Event timeline visualization
- Sequence plots
- Temporal patterns

### Config (`config`)
- Configuration management
- LifeEventsWorkflowConfig class
- Config validation

### Utils (`utils`)
- Utility functions
- Helper functions for life events
- Common operations

### Workflow (`workflow`)
- End-to-end life events analysis workflows
- Workflow orchestration
- Configuration-based execution

## Patterns

### I/O Operations
- **CRITICAL**: Always use `metainformant.core.io` for JSON/CSV/TSV operations
- Never use direct `import json` or `import csv`

**Pattern**:
```python
from metainformant.core import io

# Save event sequences and embeddings
io.dump_json(sequences, "output/life_events/workflow/sequences.json")
io.dump_json(embeddings, "output/life_events/workflow/embeddings.json")
io.write_csv(df, "output/life_events/workflow/predictions.csv")
```

### Path Handling
- Always use `metainformant.core.io.paths` utilities for path operations
- Always resolve and validate paths using `paths.expand_and_resolve()` and `paths.is_within()`

**Pattern**:
```python
from metainformant.core.io import paths

# Expand and resolve paths
work_dir = paths.expand_and_resolve(config.work_dir)
```

### Type Hints
- Comprehensive type hints throughout
- Use `from __future__ import annotations` for forward references
- Union types: `str | Path | None`

### Configuration
- Use `LifeEventsWorkflowConfig` for configuration
- Config files: `config/life_events/life_events_template.yaml`
- Environment prefix: `LE_`; only `LE_THREADS`, `LE_WORK_DIR`, and `LE_LOG_DIR` are applied (`threads`, `work_dir`, `log_dir`)
- With a `LifeEventsWorkflowConfig`, `analyze_life_course` reads only the nested `embedding`, `model`, and `workflow` sections; top-level scalar keys such as `embedding_dim` do not affect it
- Use `metainformant.core.utils.config.load_mapping_from_file()` and `apply_env_overrides()`

**Config Structure**:
```yaml
work_dir: output/life_events/workflow/work
log_dir: output/life_events/workflow/logs
threads: 4
embedding:
  embedding_dim: 100
  window_size: 5
  epochs: 10
  method: "skipgram"
  learning_rate: 0.01
model:
  model_type: "embedding"        # "embedding" or "simple"
  task_type: "classification"
  random_state: 42
```

### Embeddings
- Embeddings should be configurable (dimension, window size, method)
- Support the two embedding methods (skipgram, CBOW)
- Return embedding dictionaries

### Output Paths
- Base: `output/life_events/<workflow>/`
- Embeddings: `output/life_events/<workflow>/embeddings/`
- Models: `output/life_events/<workflow>/models/`
- Visualizations: `output/life_events/<workflow>/plots/`
- Workflows: `output/life_events/<workflow>/work/`

## Integration

- **Used by**: `phenotype` (life course phenotype extraction)
- **Integrates with**: `ml` (sequence models), `visualization` (event timeline plots)

## Reference

See Global Project Rules (root `.cursorrules`) for:
- Common directory structure and path handling
- Configuration patterns with env overrides
- Testing policy (real-implementation policy)
- Import patterns and code style
- Documentation guidelines

## Testing

- **STRICTLY REAL IMPLEMENTATION**: Test real implementations only (see Global Project Rules (root `.cursorrules`) real-implementation policy)
- Test with real event sequences from `data/`
- Write test outputs to `output/life_events/test/` using `tmp_path` fixture
- Test embedding generation with known results
- Real implementations only - use real computations and file operations


