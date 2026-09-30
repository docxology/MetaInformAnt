# AGENTS.md — `MetaInformAnt/tests/simulation`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `simulation` domain of METAINFORMANT. Tests import from `src/metainformant/simulation` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_simulation.py`
- `test_simulation_agents.py`
- `test_simulation_popgen.py`
- `test_simulation_rna_advanced.py`
- `test_simulation_workflow.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/simulation/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
