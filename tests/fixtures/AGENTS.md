# AGENTS.md — `MetaInformAnt/tests/fixtures`

Verified against disk 2026-08-30 (doc-realization fleet pass). Shared pytest fixture package for the test suite (init + domain fixture subpackages).

## Layout

- `__init__.py`
- Subdirectories: `gwas`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/fixtures/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
- Fixtures construct real in-memory/local data; keep them deterministic.
