# AGENTS.md — `MetaInformAnt/tests/life_events`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `life_events` domain of METAINFORMANT. Tests import from `src/metainformant/life_events` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_life_events.py`
- `test_life_events_cli.py`
- `test_life_events_config.py`
- `test_life_events_embeddings.py`
- `test_life_events_events.py`
- `test_life_events_integration.py`
- `test_life_events_interpretability.py`
- `test_life_events_models.py`
- `test_life_events_simulation.py`
- `test_life_events_simulation_advanced.py`
- `test_life_events_utils.py`
- `test_life_events_visualization.py`
- `test_life_events_visualization_extended.py`
- … (+1 more test modules)

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/life_events/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
