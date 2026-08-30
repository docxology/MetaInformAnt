# AGENTS.md — `MetaInformAnt/tests/mcp/`

tests for the MCP server/tools.
Files (verified 2026-08-29): test_mcp_monitor.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/mcp` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.