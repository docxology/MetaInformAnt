# AGENTS.md — `MetaInformAnt/tests/structural_variants`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `structural_variants` domain of METAINFORMANT. Tests import from `src/metainformant/structural_variants` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_structural_variants.py`
- `test_structural_variants_detection_advanced.py`
- `test_structural_variants_population.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/structural_variants/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
