# Specification: visualization

## 🎯 Scope
Ontology enrichment and network visualizations.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 3 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `plots.py`
- `visualization.py`

### Behavior Contracts

#### `visualization.py`
- `plot_go_dag(go_graph, terms=None, *, ..., layout_seed=None)`, `plot_ontology_network(ontology_graph, ..., layout_seed=None)`
  and `create_interactive_go_network(go_graph, term_annotations=None, *, ..., layout_seed=None)` pass `layout_seed` to
  `networkx.spring_layout` (and, for `plot_go_dag`, to the `random_layout` fallback). The same integer gives the same
  node positions; `None` leaves the layout unseeded.

#### `plots.py`
- `pathway_network_plot(network, output_path=None, *, ..., layout_seed=42, ...)` draws the starting positions of its
  spring layout from `random.Random(layout_seed)`: the same seed gives the same layout and the global `random` state is
  neither seeded nor consumed.
