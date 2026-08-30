# AGENTS.md — `MetaInformAnt/tests/quality/`

tests for the quality module.
Files (verified 2026-08-29): __init__.py, test_documentation_verifier.py, test_quality_contamination.py, test_quality_fastq.py, test_quality_metrics.py, test_quality_reporting.py, test_real_implementation_policy.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/quality` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.