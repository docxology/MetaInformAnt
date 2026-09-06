# AGENTS.md — `MetaInformAnt/tests/cloud/`

tests for the cloud/GCP worker scripts (mirror of `src/metainformant/cloud/`).
Files (verified 2026-09-05): __init__.py, test_cloud.py, test_download_results_depth.py, test_gcp_deployer_depth.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/cloud -p no:cacheprovider` (verified 2026-09-05: 32 passed).
Repo-wide policy: see the repository-root `AGENTS.md`.