# ml tests

pytest suite for the `ml` domain of METAINFORMANT. Tests import from `src/metainformant/ml` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_ml_automl.py`
- `test_ml_comprehensive.py`
- `test_ml_evaluation.py`
- `test_ml_features.py`
- `test_ml_interpretability.py`
- `test_ml_models.py`

Subdirectories: `llm`.

Run from the repo root: `pytest tests/ml/ -v`. Tests follow the real-implementation policy (no mocks).
