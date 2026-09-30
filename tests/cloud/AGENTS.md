# AGENTS.md — `MetaInformAnt/tests/cloud`

Verified against disk 2026-08-30 (doc-realization fleet pass). Cloud-related tests (`test_cloud.py`) for METAINFORMANT's cloud helper surface.

## Layout

- `__init__.py`
- `test_cloud.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/cloud/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
