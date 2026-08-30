# AGENTS.md — `MetaInformAnt/tests/simulation/`

tests for the simulation module.
Files (verified 2026-08-29): __init__.py, test_simulation.py, test_simulation_agents.py, test_simulation_popgen.py, test_simulation_rna_advanced.py, test_simulation_workflow.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/simulation` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.