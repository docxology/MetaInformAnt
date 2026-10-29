# AGENTS.md — `MetaInformAnt/tests/ml/`

tests for the ML module.
Files (verified 2026-10-28): `__init__.py`, `test_ml_automl.py`, `test_ml_automl_cross_validation.py`,
`test_ml_comprehensive.py`, `test_ml_deep_learning.py`, `test_ml_dimensionality.py`,
`test_ml_ensemble_parallel_jobs.py`, `test_ml_evaluation.py`, `test_ml_evaluation_arguments.py`,
`test_ml_feature_selection_depth.py`, `test_ml_feature_selection_exact_tests.py`, `test_ml_features.py`,
`test_ml_interpretability.py`, `test_ml_models.py`, `test_ml_seeded_randomness.py`,
`test_ml_sequence_cnn_short_inputs.py`, `test_ml_stability_selection_failures.py`,
`test_ml_validation_strategies_and_ensemble_cv.py` (17 test files plus `__init__.py`).

## Children (documented on disk)
- llm/ (test_ollama_chains.py, test_ollama_client.py, test_ollama_config_retries.py, test_ollama_local_server.py, test_ollama_prompts.py)

## Coverage notes
- `test_ml_seeded_randomness.py`: the routines that draw random numbers take a seed and neither seed nor consume the global
  random state: `evaluation.validation` (`permutation_importance_biological`, `validate_model_stability`,
  `bootstrap_validate`, `learning_curve`), `models.regression.analyze_prediction_uncertainty`, `automl.optimization`
  (`random_search`, `bayesian_optimization`), `interpretability.explainers` (`compute_permutation_importance`,
  `compute_shap_values_kernel`, `compute_lime_explanation`) and `interpretability.feature_selection` (`boruta_selection`,
  including the seed 0, and `stability_selection`). Oracles are exact where the maths allows (scikit-learn
  cross-validated accuracy, exact Shapley values of a linear model) and structural otherwise.
- `test_ml_validation_strategies_and_ensemble_cv.py`: `compare_validation_strategies` runs the `"stratified_5fold"`
  strategy (oracle: scikit-learn's `cross_val_score` over `StratifiedKFold(5, shuffle=True, random_state)`), and
  `cross_validate_biological(method="ensemble", random_state=<int>)` runs (oracle: a hand-assembled soft-voting ensemble
  with the same members).
- `test_ml_ensemble_parallel_jobs.py`: the ensemble in `cross_validate_biological` follows `parallel_jobs` unless an
  explicit `n_jobs` is given; with `parallel_jobs=1` no worker process is started (checked in a subprocess).
- `test_ml_sequence_cnn_short_inputs.py`: `conv1d_forward` and `predict_sequences` raise a `ValueError` that names the
  sequence length and the kernel size, or the minimum length of the whole stack, for sequences that are too short; the
  convolution oracle is a direct sliding-window computation.
- `llm/test_ollama_config_retries.py`: `OllamaConfig` rejects `max_retries` below 1, the client refuses a value lowered
  after construction before it touches the network, and the retry loop is checked against `file://` URLs (no server and
  no network).

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Ollama integration tests in `llm/test_ollama_client.py` and
  `llm/test_ollama_chains.py` skip unless `METAINFORMANT_RUN_OLLAMA_INTEGRATION=1`.
- Tests that train the voting ensemble with the default `n_jobs=-1` start joblib worker processes and need an
  environment that allows process pools; the regression tests above pass `n_jobs=1` or `parallel_jobs=1`.
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/ml -p no:cacheprovider`.
Repo-wide policy: see the repository-root `AGENTS.md`.
