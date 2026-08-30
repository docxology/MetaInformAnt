# AGENTS.md — `MetaInformAnt/tests/menu/`

tests for the interactive menu.
Files (verified 2026-08-29): __init__.py, test_menu_discovery.py, test_menu_display.py, test_menu_executor.py, test_menu_navigation.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/menu` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.