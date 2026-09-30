# AGENTS.md — `MetaInformAnt/tests/ml`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `ml` domain of METAINFORMANT. Tests import from `src/metainformant/ml` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_ml_automl.py`
- `test_ml_comprehensive.py`
- `test_ml_evaluation.py`
- `test_ml_features.py`
- `test_ml_interpretability.py`
- `test_ml_models.py`
- Subdirectories: `llm`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/ml/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
