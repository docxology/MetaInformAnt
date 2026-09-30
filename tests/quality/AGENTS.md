# AGENTS.md — `MetaInformAnt/tests/quality`

Verified against disk 2026-08-30 (doc-realization fleet pass). pytest suite for the `quality` domain of METAINFORMANT. Tests import from `src/metainformant/quality` and follow the repo's real-implementation policy.

## Layout

- `__init__.py`
- `test_documentation_verifier.py`
- `test_quality_contamination.py`
- `test_quality_fastq.py`
- `test_quality_metrics.py`
- `test_quality_reporting.py`
- `test_real_implementation_policy.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/quality/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
