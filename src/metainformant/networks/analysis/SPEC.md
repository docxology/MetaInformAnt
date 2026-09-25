# Specification: analysis

## 🎯 Scope
Network analysis subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 6 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `community.py`
- `graph.py`
- `graph_algorithms.py`
- `graph_core.py`
- `pathway.py`

### Behavior Contracts

#### `graph_algorithms.py`
- `shortest_paths(graph, source=None, target=None)` — returns
  `{source: {target: length}}`. Unreachable pairs are omitted by every branch
  (NetworkX and BiologicalNetwork alike; the former inf-filling wrapper branch
  was removed); a single-pair query with no path returns `{}`. Missing
  `source`/`target` nodes raise `ValueError` naming the node; other errors
  propagate.
