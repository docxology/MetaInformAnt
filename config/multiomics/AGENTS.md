# Agent Directives: config/multiomics

## Role

Multi-omics integration pipeline configuration for cross-platform data harmonization.

## Contents

| File | Description |
| :--- | :--- |
| `multiomics_template.yaml` | Template with integration method options |

## Configuration Structure

```yaml
# Multi-omics configuration
layers: [genomics, transcriptomics, proteomics]
integration_method: joint_pca
n_components: 50
standardize: true
```

## Rules

- Validate with schema before committing new configs
- Follow REAL IMPLEMENTATION policy — tests use real config files
- Use `uv` for dependency management
- Environment overrides use `AK_` prefix
