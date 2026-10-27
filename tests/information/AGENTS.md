# AGENTS.md — `MetaInformAnt/tests/information/`

tests for the information-theory module.
Files (verified 2026-10-26): `__init__.py`, `test_continuous_entropy.py`, `test_information_chao_shen.py`,
`test_information_comprehensive.py`, `test_information_estimation_depth.py`,
`test_information_geometry_decomposition.py`, `test_information_integration.py`,
`test_information_new_modules.py`, `test_information_sequence_scans.py`
(8 test files plus `__init__.py`).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/information -p no:cacheprovider`
  (verified 2026-10-26: 482 passed).
- `test_continuous_entropy.py` pins the continuous (differential) measures in nats: histogram and
  Kozachenko-Leonenko entropy, and mutual information, conditional entropy and transfer entropy in the joint space.
- `test_information_sequence_scans.py` checks the repeat and palindrome scans behind `analyze_sequence_information`
  against hand-worked records and brute-force references, and bounds the run time on 2 kb inputs (5 s) and a 600 bp
  tandem repeat (10 s); its `pytest.mark.timeout` guards need `pytest-timeout`.
Repo-wide policy: see the repository-root `AGENTS.md`.