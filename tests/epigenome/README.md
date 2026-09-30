# epigenome tests

pytest suite for the `epigenome` domain of METAINFORMANT. Tests import from `src/metainformant/epigenome` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_epigenome.py`
- `test_epigenome_analysis.py`
- `test_epigenome_assays.py`
- `test_epigenome_chromatin.py`
- `test_epigenome_peak_calling.py`
- `test_epigenome_visualization.py`
- `test_epigenome_workflow.py`

Run from the repo root: `pytest tests/epigenome/ -v`. Tests follow the real-implementation policy (no mocks).
