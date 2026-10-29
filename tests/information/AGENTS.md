# AGENTS.md — `MetaInformAnt/tests/information/`

tests for the information-theory module.
Files (verified 2026-10-28): `__init__.py`, `test_continuous_entropy.py`, `test_information_chao_shen.py`,
`test_information_comprehensive.py`, `test_information_estimation_depth.py`,
`test_information_geometry_decomposition.py`, `test_information_integration.py`, `test_information_new_modules.py`,
`test_information_seeded_randomness.py`, `test_information_sequence_scans.py`,
`test_information_units_and_options.py`, `test_integration_failed_platform.py`,
`test_integration_network_entropy_von_neumann.py`, `test_integration_network_flow_and_motifs.py`,
`test_network_entropy_von_neumann.py` (14 test files plus `__init__.py`).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/information -p no:cacheprovider`
  (verified 2026-10-27: 513 passed; the run used `-p no:rerunfailures` to keep the plugin from binding a local socket).
- `test_continuous_entropy.py` pins the continuous (differential) measures in nats: histogram and
  Kozachenko-Leonenko entropy, and mutual information, conditional entropy and transfer entropy in the joint space.
- `test_information_sequence_scans.py` checks the repeat and palindrome scans behind `analyze_sequence_information`
  against hand-worked records and brute-force references, and bounds the run time on 2 kb inputs (5 s) and a 600 bp
  tandem repeat (10 s); its `pytest.mark.timeout` guards need `pytest-timeout`.
- `test_information_seeded_randomness.py`: `transfer_entropy`, `information_flow_network`, `mutual_information_network`
  and `entropy_bootstrap_confidence` take a seed (an integer, or a `random.Random` for the permutation tests); the same
  seed reproduces the output, different seeds differ, and the global `random` / `numpy.random` states are identical
  before and after. The non-random quantities (the TE value, the MI matrix) are checked against
  `sklearn.metrics.mutual_info_score` on separately binned data.
- `test_network_entropy_von_neumann.py`: `information_flow.network_entropy` is the Von Neumann entropy of the Laplacian
  density matrix `rho = L / trace(L)` (closed forms for complete graphs and stars, and NumPy's eigensolver on
  `networkx`'s Laplacian), depends on the structure of the network, ignores self-loops, uses absolute edge weights, and
  raises `ValueError` for a matrix that is not square or not symmetric; the pure-Python Jacobi solver is checked too.
Repo-wide policy: see the repository-root `AGENTS.md`.