# longread tests

pytest suite for the `longread` domain of METAINFORMANT. Tests import from `src/metainformant/longread` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_longread.py`
- `test_longread_analysis.py`
- `test_longread_assembly.py`
- `test_longread_io.py`
- `test_longread_methylation.py`
- `test_longread_quality.py`
- `test_longread_visualization.py`

Run from the repo root: `pytest tests/longread/ -v`. Tests follow the real-implementation policy (no mocks).
