# Specification: interaction

## 🎯 Scope
Network interaction subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 6 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `_ppi_impl.py`
- `ppi.py`
- `regulatory.py`
- `regulatory_analysis.py`
- `regulatory_core.py`

### Behavior Contracts

#### `_ppi_impl.py` (re-exported by `ppi.py`)
- `load_ppi_network(path, format="tsv")` and `ProteinNetwork.from_file(path, format="tsv")` accept `'tsv'`, `'csv'`,
  `'json'`, `'bioplex'` and `'intact'`; any other format raises `ValueError`. `'json'` reads the node-link layout that
  `save_ppi_network(..., format="json")` writes (edge list under `edges` or, for older networkx releases, `links`),
  keeps node and edge attributes, gives an edge without a `weight` the weight 1.0, and returns an undirected graph; a
  file that is not a node-link object, a node without an `id`, or an edge without `source` or `target` raises
  `ValueError`.
- `ppi_network_analysis(graph)` — the `centrality` section holds `degree_centrality` and `betweenness_centrality`
  limited to the ten most central nodes, most central first, ties ordered by `str(node)`.

#### `regulatory_analysis.py`, `regulatory_core.py` (re-exported by `regulatory.py`)
- `infer_grn(expression_data, gene_names, method="correlation", threshold=0.5, tf_genes=None, random_state=None)` —
  `method` is `'correlation'` or `'mutual_info'`; anything else (including `'granger'`) raises `ValueError`.
  `infer_regulatory_network_from_expression(expression_data, regulators, method, threshold, random_state=None)`
  supports the same two methods, raises `NotImplementedError` for `'grnboost2'` and `ValueError` for any other name.
- `random_state` (an integer seed or a `numpy.random.RandomState`) feeds one generator that every
  `mutual_info_regression` call of a run shares, consumed pair by pair in the order the regulators and targets are
  visited. The global NumPy random state is neither seeded nor consumed; `'correlation'` is deterministic.
- A feed-forward loop is TF1 -> TF2, TF2 -> target and the direct TF1 -> target edge, in `regulatory_motifs()` and in
  `GeneRegulatoryNetwork.regulatory_motifs()`; a chain without the direct edge is not reported.
- `analyze_regulatory_dynamics(graph, initial_states, time_steps=10, random_state=None)` draws its step noise from
  `numpy.random.default_rng(random_state)` (an integer seed or a `numpy.random.Generator`).
