# AGENTS.md — `MetaInformAnt/tests/life_events/`

tests for the life events module.
Files (verified 2026-10-27): __init__.py, test_life_events.py, test_life_events_cli.py, test_life_events_compare_populations_plot_errors.py, test_life_events_config.py, test_life_events_cox_ph_model.py, test_life_events_domain_embeddings.py, test_life_events_embedding_training.py, test_life_events_embeddings.py, test_life_events_events.py, test_life_events_generation_loading.py, test_life_events_integration.py, test_life_events_interpretability.py, test_life_events_intervention_outcomes.py, test_life_events_models.py, test_life_events_multitask_predictor.py, test_life_events_predictor_simple.py, test_life_events_script_options.py, test_life_events_sequence_models_torch.py, test_life_events_shapley.py, test_life_events_shapley_background.py, test_life_events_simulation.py, test_life_events_simulation_advanced.py, test_life_events_survival.py, test_life_events_survival_predictor.py, test_life_events_temporal_embeddings.py, test_life_events_utils.py, test_life_events_visualization.py, test_life_events_visualization_extended.py, test_life_events_workflow.py, test_life_events_workflow_event_importance.py, test_life_events_workflow_outcome_plot.py, test_life_events_workflow_plot_failures.py (32 test files).

## Contract tests

- `test_life_events_workflow_event_importance.py`: `event_importance(method="frequency")` ranks by occurrence count, `"temporal"` is the earliness score, `"transition"` counts transitions, and an unknown method raises `ValueError`.
- `test_life_events_intervention_outcomes.py`: `intervention_analysis(outcomes=...)` requires one value per pre-intervention sequence followed by one per post-intervention sequence and raises `ValueError` naming both counts otherwise.
- `test_life_events_compare_populations_plot_errors.py`: a failed comparison plot is recorded in `visualization_errors` (provoked by a directory squatting on the plot's file name) and the comparison is still returned.
- `test_life_events_temporal_embeddings.py`: `biological_embedding(..., "temporal")` works on real `EventSequence` objects (intervals in days), matches hand-computed statistics, and raises `ValueError` for an `embedding_dim` below 6.
- `test_life_events_shapley_background.py`: a non-empty `background_sequences` of `shapley_values` emits a `DeprecationWarning` and is ignored.

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/life_events -p no:cacheprovider`
Repo-wide policy: see the repository-root `AGENTS.md`.
