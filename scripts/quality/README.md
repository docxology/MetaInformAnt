# QUALITY

## Overview
Command-line helpers for Quality workflows. Scripts should remain thin wrappers around `src/metainformant/` implementations and be run from the repository root with `uv`.

## Contents
- [audit_docs.py](audit_docs.py)
- [audit_docstrings.py](audit_docstrings.py)
- [audit_repo.py](audit_repo.py)
- [check_exports.py](check_exports.py)
- [fix_mermaid_violations.py](fix_mermaid_violations.py)
- [generate_docs_coverage.py](generate_docs_coverage.py)
- [run_quality_control.py](run_quality_control.py)

## Usage
```bash
uv run python scripts/quality/audit_docs.py --help
```
