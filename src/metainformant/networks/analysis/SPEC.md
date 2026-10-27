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

#### `community.py`
- `louvain_communities(graph, resolution=1.0, randomize=True, random_state=None, **kwargs)` — with a seed (and
  `randomize=True`) it calls python-louvain's `best_partition(..., random_state=seed)` and does not pass `randomize`;
  `randomize=False` is passed on and ignores the seed. The global `random` and `numpy.random` states are never
  seeded.
- `detect_communities(graph, method, **kwargs)` — when the chosen method raises `ImportError`, Leiden falls back to
  `louvain_communities` (if python-louvain is installed) with the keyword arguments restricted to `resolution`,
  `randomize`, `random_state`, `weight` and `partition`; every other method, and Leiden without python-louvain, falls
  back to `greedy_modularity_communities`, which keeps `weight`, `resolution`, `cutoff` and `best_n`. Any other
  exception falls back to greedy modularity called without options; `method="greedy"` re-raises.

#### `pathway.py`
- `pathway_enrichment(gene_list, pathway_network, ..., method, correction)` — raises `ValueError` for a `method` other
  than `'fisher'` or `'hypergeometric'` (checked first). The p-value is the one-sided `P(X >= overlap)`: from
  `scipy.stats` when scipy is importable, otherwise the exact upper tail summed in integer arithmetic. `correction`
  `'fdr'` (Benjamini-Hochberg) is pure Python and does not depend on scipy.
- `pathway_similarity_analysis(pathways, method)` — `similarities` is keyed `"{name1}_{name2}"`, which is ambiguous when
  a pathway ID contains an underscore; `create_pathway_network` therefore builds its edges from the scored pairs and
  never parses the keys. A pathway with no genes has similarity 0.0 under `'jaccard'`, `'overlap'` and `'cosine'`.
