# AGENTS.md — `MetaInformAnt/tests/longread`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `longread` domain of METAINFORMANT. Tests import from `src/metainformant/longread` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_longread.py`
- `test_longread_analysis.py`
- `test_longread_assembly.py`
- `test_longread_io.py`
- `test_longread_methylation.py`
- `test_longread_quality.py`
- `test_longread_visualization.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/longread/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
