# AGENTS.md — `MetaInformAnt/tests/multiomics`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `multiomics` domain of METAINFORMANT. Tests import from `src/metainformant/multiomics` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_multiomics_comprehensive.py`
- `test_multiomics_integration.py`
- `test_multiomics_methods.py`
- `test_multiomics_pathways.py`
- `test_multiomics_sample_mapping.py`
- `test_multiomics_survival.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/multiomics/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
