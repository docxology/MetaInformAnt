# AGENTS.md — `MetaInformAnt/tests/phenotype`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `phenotype` domain of METAINFORMANT. Tests import from `src/metainformant/phenotype` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_multivariate.py`
- `test_phenotype_basic.py`
- `test_phenotype_behavior.py`
- `test_phenotype_comprehensive.py`
- `test_phenotype_integration.py`
- `test_phenotype_life_course.py`
- `test_phenotype_modules.py`
- `test_phenotype_morphological.py`
- `test_phenotype_scraper.py`
- `test_phenotype_workflow.py`
- `test_plots.py`
- `test_statistical.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/phenotype/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
