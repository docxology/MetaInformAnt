# AGENTS.md — `MetaInformAnt/tests/ecology`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `ecology` domain of METAINFORMANT. Tests import from `src/metainformant/ecology` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_ecology_basic.py`
- `test_ecology_comprehensive.py`
- `test_ecology_functional.py`
- `test_ecology_macroecology.py`
- `test_ecology_ordination.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/ecology/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
