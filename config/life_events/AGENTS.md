# Agent Directives: config/life_events

## Role

Life events temporal analysis configuration for life course modeling.

## Contents

| File | Description |
| :--- | :--- |
| `life_events_template.yaml` | Template with all options documented |

## Configuration Structure

```yaml
# Life events configuration
embedding_dim: 100
event_domains: [education, occupation, health]
prediction_target: mortality
```

## Rules

- Validate with schema before committing new configs
- Follow REAL IMPLEMENTATION policy — tests use real config files
- Use `uv` for dependency management
- Environment overrides use `AK_` prefix
