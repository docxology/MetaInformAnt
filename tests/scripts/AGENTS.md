# AGENTS.md — `MetaInformAnt/tests/scripts`

Verified against disk 2026-08-30 (doc-realization fleet pass). Tests for utility scripts (`test_verify_utils.py`).

## Layout

- `test_verify_utils.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/scripts/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
