# infrastructure tests

Meta-tests of the repo itself: build, CLI, dependency verification, repo structure, and test-infrastructure checks.

## Files

- `__init__.py`
- `test_build.py`
- `test_cli.py`
- `test_dependency_verifier.py`
- `test_repo_structure.py`
- `test_test_infrastructure.py`

Run from the repo root: `pytest tests/infrastructure/ -v`. Tests follow the real-implementation policy (no mocks).
