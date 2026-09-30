# AGENTS.md — `MetaInformAnt/tests/information`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `information` domain of METAINFORMANT. Tests import from `src/metainformant/information` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_information_comprehensive.py`
- `test_information_geometry_decomposition.py`
- `test_information_integration.py`
- `test_information_new_modules.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/information/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
