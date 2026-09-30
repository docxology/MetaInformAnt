# AGENTS.md — `MetaInformAnt/tests/pharmacogenomics`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `pharmacogenomics` domain of METAINFORMANT. Tests import from `src/metainformant/pharmacogenomics` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_pharmacogenomics_alleles.py`
- `test_pharmacogenomics_annotations.py`
- `test_pharmacogenomics_clinical.py`
- `test_pharmacogenomics_metabolism.py`
- `test_pharmacogenomics_visualization.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/pharmacogenomics/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
