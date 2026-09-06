# AGENTS.md — `MetaInformAnt/tests/ml/`

tests for the ML module.
Files (verified 2026-09-05): __init__.py, test_ml_automl.py, test_ml_comprehensive.py, test_ml_evaluation.py, test_ml_features.py, test_ml_interpretability.py, test_ml_models.py, test_ml_deep_learning.py, test_ml_feature_selection_depth.py.


## Children (documented on disk)
- llm/ (test_ollama_chains.py, test_ollama_client.py, test_ollama_local_server.py, test_ollama_prompts.py)

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Ollama integration tests in `llm/test_ollama_client.py` and
  `llm/test_ollama_chains.py` skip unless `METAINFORMANT_RUN_OLLAMA_INTEGRATION=1`.
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/ml -p no:cacheprovider`.
Repo-wide policy: see the repository-root `AGENTS.md`.
