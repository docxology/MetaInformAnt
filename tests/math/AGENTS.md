# AGENTS.md — `MetaInformAnt/tests/math`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `math` domain of METAINFORMANT. Tests import from `src/metainformant/math` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_math.py`
- `test_math_bayesian.py`
- `test_math_coalescent.py`
- `test_math_coalescent_expectations.py`
- `test_math_coalescent_extras.py`
- `test_math_comprehensive.py`
- `test_math_decision.py`
- `test_math_demography.py`
- `test_math_drift_migration.py`
- `test_math_dynamics.py`
- `test_math_effective_size_extras.py`
- `test_math_egt_epi_fst_ne.py`
- `test_math_epidemiology.py`
- … (+12 more test modules)

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/math/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
