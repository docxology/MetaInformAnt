# AGENTS.md — `MetaInformAnt/tests/singlecell`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `singlecell` domain of METAINFORMANT. Tests import from `src/metainformant/singlecell` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_singlecell_basic.py`
- `test_singlecell_celltyping.py`
- `test_singlecell_differential.py`
- `test_singlecell_dimensionality.py`
- `test_singlecell_preprocessing.py`
- `test_singlecell_velocity.py`
- `test_singlecell_visualization.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/singlecell/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
