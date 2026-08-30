# AGENTS.md — `MetaInformAnt/tests/ml/`

tests for the ML module.
Files (verified 2026-08-29): __init__.py, test_ml_automl.py, test_ml_comprehensive.py, test_ml_evaluation.py, test_ml_features.py, test_ml_interpretability.py, test_ml_models.py.


## Children (documented on disk)
- llm/

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/ml` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.