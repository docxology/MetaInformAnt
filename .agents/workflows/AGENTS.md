# AGENTS.md — `MetaInformAnt/.agents/workflows`

Verified against disk 2026-08-30 (doc-realization fleet pass). Named agent workflows for this repo (plain markdown recipes): machine setup, running the test suite, pipeline monitoring, pushing the repo.

## Layout

- `run_tests.md` — how to run the real-implementation test suite (`uv run pytest tests/ -v`, per-file, coverage).
- `monitor_pipeline.md` (purpose per file; unverified beyond name).
- `push_repo.md`, `setup_machine.md` (purpose per file; unverified beyond name).

## Gotchas

- Keep each recipe runnable as written; update when commands change.
