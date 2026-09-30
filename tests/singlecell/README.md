# singlecell tests

pytest suite for the `singlecell` domain of METAINFORMANT. Tests import from `src/metainformant/singlecell` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_singlecell_basic.py`
- `test_singlecell_celltyping.py`
- `test_singlecell_differential.py`
- `test_singlecell_dimensionality.py`
- `test_singlecell_preprocessing.py`
- `test_singlecell_velocity.py`
- `test_singlecell_visualization.py`

Run from the repo root: `pytest tests/singlecell/ -v`. Tests follow the real-implementation policy (no mocks).
