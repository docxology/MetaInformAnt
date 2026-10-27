# Specification: workflow

## 🎯 Scope
Life events workflow sub-package.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 2 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `workflow.py`

## Behavioral Contracts

- `event_importance(sequences, method="frequency", normalize=True)`: `"frequency"` ranks event types by their number of occurrences, `"temporal"` by earliness (`1 - mean normalised position`, 0 for the first event and 1 for the last) and `"transition"` by the number of transitions into and out of the type. An unknown `method` raises `ValueError`.
- `intervention_analysis(...)`: every sequence is split at `intervention_time` into a pre-intervention part (events before it) and a post-intervention part (events at or after it). `outcomes` must hold exactly one value per pre-intervention sequence followed by one per post-intervention sequence, otherwise `ValueError`; `pre_intervention_outcomes` and `post_intervention_outcomes` give each side its own list and take precedence over `outcomes`.
- `compare_populations(...)`: the result lists the plots that were written in `visualizations` and the plots that failed in `visualization_errors` (`{"type", "error_type", "error"}`, an empty list when the plot was written). A failed plot does not stop the comparison.
