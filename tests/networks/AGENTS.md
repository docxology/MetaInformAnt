# AGENTS.md — `MetaInformAnt/tests/networks/`

tests for the networks module.
Files (verified 2026-08-29): __init__.py, test_networks_community.py, test_networks_comprehensive.py, test_networks_graph.py, test_networks_pathway.py, test_networks_ppi.py, test_networks_regulatory.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/networks` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.