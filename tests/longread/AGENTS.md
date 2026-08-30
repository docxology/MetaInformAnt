# AGENTS.md — `MetaInformAnt/tests/longread/`

tests for the longread module.
Files (verified 2026-08-29): __init__.py, test_longread.py, test_longread_analysis.py, test_longread_assembly.py, test_longread_io.py, test_longread_methylation.py, test_longread_quality.py, test_longread_visualization.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/longread` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.