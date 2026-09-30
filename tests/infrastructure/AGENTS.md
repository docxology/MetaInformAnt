# AGENTS.md — `MetaInformAnt/tests/infrastructure`

Verified against disk 2026-08-30 (doc-realization fleet pass). Meta-tests of the repo itself: build, CLI, dependency verification, repo structure, and test-infrastructure checks.

## Layout

- `__init__.py`
- `test_build.py`
- `test_cli.py`
- `test_dependency_verifier.py`
- `test_repo_structure.py`
- `test_test_infrastructure.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/infrastructure/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
- These test the testing/build scaffolding, not scientific code.
