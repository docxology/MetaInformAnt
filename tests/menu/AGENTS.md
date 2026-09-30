# AGENTS.md — `MetaInformAnt/tests/menu`

Verified against disk 2026-08-30 (doc-realization fleet pass). Tests for the interactive menu system: discovery, display, executor, navigation.

## Layout

- `__init__.py`
- `test_menu_discovery.py`
- `test_menu_display.py`
- `test_menu_executor.py`
- `test_menu_navigation.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/menu/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
