# AGENTS.md — `MetaInformAnt/tests/epigenome/`

tests for the epigenome module.
Files (verified 2026-08-29): __init__.py, test_epigenome.py, test_epigenome_analysis.py, test_epigenome_assays.py, test_epigenome_chromatin.py, test_epigenome_peak_calling.py, test_epigenome_visualization.py, test_epigenome_workflow.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/epigenome` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.