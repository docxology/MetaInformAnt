# Agent Directives: config/networks

## Role

Biological network analysis pipeline configuration for graph construction and community detection.

## Contents

| File | Description |
| :--- | :--- |
| `networks_template.yaml` | Template with graph algorithm options |

## Configuration Structure

```yaml
# Network analysis configuration
graph_type: undirected
community_method: louvain
centrality_measures: [degree, betweenness]
min_edge_weight: 0.5
```

## Rules

- Validate with schema before committing new configs
- Follow REAL IMPLEMENTATION policy — tests use real config files
- Use `uv` for dependency management
- Environment overrides use `AK_` prefix
