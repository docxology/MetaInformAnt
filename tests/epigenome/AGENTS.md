# AGENTS.md — `MetaInformAnt/tests/epigenome`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `epigenome` domain of METAINFORMANT. Tests import from `src/metainformant/epigenome` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_epigenome.py`
- `test_epigenome_analysis.py`
- `test_epigenome_assays.py`
- `test_epigenome_chromatin.py`
- `test_epigenome_peak_calling.py`
- `test_epigenome_visualization.py`
- `test_epigenome_workflow.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/epigenome/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
