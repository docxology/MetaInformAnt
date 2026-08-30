# AGENTS.md — `MetaInformAnt/tests/metagenomics/`

tests for the metagenomics module.
Files (verified 2026-08-29): __init__.py, test_metagenomics_amplicon.py, test_metagenomics_comparative.py, test_metagenomics_diversity.py, test_metagenomics_functional.py, test_metagenomics_shotgun.py, test_metagenomics_visualization.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/metagenomics` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.