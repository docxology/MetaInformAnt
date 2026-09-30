# AGENTS.md — `MetaInformAnt/tests/spatial`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `spatial` domain of METAINFORMANT. Tests import from `src/metainformant/spatial` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_spatial_analysis.py`
- `test_spatial_autocorrelation.py`
- `test_spatial_communication.py`
- `test_spatial_deconvolution_advanced.py`
- `test_spatial_integration.py`
- `test_spatial_io.py`
- `test_spatial_neighborhood.py`
- `test_spatial_visualization.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/spatial/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
