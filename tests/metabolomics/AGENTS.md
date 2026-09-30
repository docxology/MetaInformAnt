# AGENTS.md — `MetaInformAnt/tests/metabolomics`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `metabolomics` domain of METAINFORMANT. Tests import from `src/metainformant/metabolomics` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_metabolomics_identification.py`
- `test_metabolomics_io.py`
- `test_metabolomics_pathways.py`
- `test_metabolomics_visualization.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/metabolomics/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
