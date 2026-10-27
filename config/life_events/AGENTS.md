# Agent Directives: config/life_events

## Role

Life events temporal analysis configuration for life course modeling.

## Contents

| File | Description |
| :--- | :--- |
| `life_events_template.yaml` | Template with all options documented |

## Configuration Structure

```yaml
# Life events configuration (see life_events_template.yaml)
work_dir: output/life_events/work
embedding:
  embedding_dim: 100
  window_size: 5
  method: skipgram
model:
  model_type: embedding      # "embedding" or "simple"
  task_type: classification
  random_state: 42
```

## Rules

- Validate with schema before committing new configs
- Follow REAL IMPLEMENTATION policy — tests use real config files
- Use `uv` for dependency management
- Environment overrides use the life-events namespace documented by the module
