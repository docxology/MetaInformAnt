# AGENTS.md — `MetaInformAnt/.agents/workflows/`

Four markdown playbooks (verified 2026-08-29): `setup_machine.md`,
`run_tests.md`, `monitor_pipeline.md`, `push_repo.md`. Machine setup, test
running, pipeline monitoring, and repo push procedures for agents working in
this repo.

## Gotchas
- `push_repo.md` performs git writes — only follow it with explicit owner
  authorization; doc passes must never push.
Repo-wide policy: see the repository-root `AGENTS.md`.