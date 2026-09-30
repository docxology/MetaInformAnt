# AGENTS.md — `MetaInformAnt/config/quality`

Verified against disk 2026-08-30 (doc-realization fleet pass). Quality-gate configuration: the mypy error budget for the codebase.

## Layout

- `mypy_error_budget.txt` — current allowed mypy error count/baseline used by quality checks.

## Gotchas

- Reduce the budget as errors are fixed; do not raise it casually.
