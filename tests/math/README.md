# math tests

pytest suite for the `math` domain of METAINFORMANT. Tests import from `src/metainformant/math` and follow the repo's real-implementation policy.

## Files

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

Run from the repo root: `pytest tests/math/ -v`. Tests follow the real-implementation policy (no mocks).
