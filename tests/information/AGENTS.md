# AGENTS.md — `MetaInformAnt/tests/information/`

tests for the information-theory module.
Files (verified 2026-09-05): `__init__.py`, `test_information_chao_shen.py`,
`test_information_comprehensive.py`, `test_information_estimation_depth.py`,
`test_information_geometry_decomposition.py`, `test_information_integration.py`,
`test_information_new_modules.py`.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/information -p no:cacheprovider`
  (verified 2026-09-05: 376 passed).
Repo-wide policy: see the repository-root `AGENTS.md`.