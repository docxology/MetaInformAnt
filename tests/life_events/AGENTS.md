# AGENTS.md — `MetaInformAnt/tests/life_events/`

tests for the life-events module.
Files (verified 2026-08-29): __init__.py, test_life_events.py, test_life_events_cli.py, test_life_events_config.py, test_life_events_embeddings.py, test_life_events_events.py, test_life_events_integration.py, test_life_events_interpretability.py … (15 files).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/life_events` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.