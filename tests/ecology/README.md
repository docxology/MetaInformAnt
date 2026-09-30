# ecology tests

pytest suite for the `ecology` domain of METAINFORMANT. Tests import from `src/metainformant/ecology` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_ecology_basic.py`
- `test_ecology_comprehensive.py`
- `test_ecology_functional.py`
- `test_ecology_macroecology.py`
- `test_ecology_ordination.py`

Run from the repo root: `pytest tests/ecology/ -v`. Tests follow the real-implementation policy (no mocks).
