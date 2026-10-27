# Specification: regulatory

## 🎯 Scope
Regulatory network inference and motif analysis subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 3 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `grn_inference.py`
- `motif_analysis.py`

### Behavior Contracts

#### `grn_inference.py`
- `compute_network_motifs(edges, motif_size=3, random_state=None)` — counts motifs of a directed edge list and
  compares them with 100 randomised networks (shuffled targets). `random_state` is an integer seed or a
  `random.Random`; the same seed gives the same z-scores and the global `random` state is neither seeded nor consumed.
  For `motif_size=3` a `feed_forward_loop` is A -> B, B -> C and the direct A -> C; A -> B -> C without A -> C is a
  `cascade`. `motif_size` other than 2 or 3 raises `ValueError`.
- `infer_grn()` and `infer_regulatory_network_from_expression()` are in `networks.interaction.regulatory_analysis`
  and support only the methods `'correlation'` and `'mutual_info'` (`'grnboost2'` raises `NotImplementedError`).
