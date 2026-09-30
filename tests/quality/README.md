# quality tests

pytest suite for the `quality` domain of METAINFORMANT. Tests import from `src/metainformant/quality` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_documentation_verifier.py`
- `test_quality_contamination.py`
- `test_quality_fastq.py`
- `test_quality_metrics.py`
- `test_quality_reporting.py`
- `test_real_implementation_policy.py`

Run from the repo root: `pytest tests/quality/ -v`. Tests follow the real-implementation policy (no mocks).
