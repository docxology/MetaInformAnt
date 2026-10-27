# AGENTS.md — `MetaInformAnt/tests/quality/`

tests for the quality module.
Files (verified 2026-10-27): __init__.py, test_contamination_depth.py, test_doc_calls_and_links.py, test_doc_imports.py, test_documentation_verifier.py, test_fastq_depth.py, test_quality_batch_analysis.py, test_quality_batch_combat.py, test_quality_batch_detection.py, test_quality_contamination.py, test_quality_fastq.py, test_quality_metrics.py, test_quality_reporting.py, test_real_implementation_policy.py (13 test files).

## Contract tests

- `test_quality_batch_analysis.py`: `batch_quality_analysis` accepts only `data_type="fastq"`; `"vcf"` and `"bam"` raise `NotImplementedError` and any other value raises `ValueError`, before a file is touched and even for an empty list.
- `test_quality_batch_combat.py`: `correct_batch_combat` is a parametric empirical-Bayes ComBat that does not depend on batch order or labels.
- `test_doc_calls_and_links.py`: calls quoted in documents are bound against the current signatures, and relative links resolve against the real directory trees.

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/quality`
Repo-wide policy: see the repository-root `AGENTS.md`.
