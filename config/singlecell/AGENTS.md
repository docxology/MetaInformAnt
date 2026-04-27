# Agent Directives: config/singlecell

## Role

Single-cell RNA-seq analysis pipeline configuration for preprocessing, clustering, and trajectory inference.

## Contents

| File | Description |
| :--- | :--- |
| `singlecell_template.yaml` | Template with preprocessing and clustering options |

## Configuration Structure

```yaml
# Single-cell configuration
min_genes: 200
min_cells: 3
normalization: total
clustering_resolution: 0.5
n_pcs: 50
```

## Rules

- Validate with schema before committing new configs
- Follow NO MOCKING policy — tests use real config files
- Use `uv` for dependency management
- Environment overrides use `AK_` prefix
