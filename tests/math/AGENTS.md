# AGENTS.md — `MetaInformAnt/tests/math/`

tests for the math module (population genetics, Bayesian, decision theory…).
Files (verified 2026-08-29): __init__.py, test_math.py, test_math_bayesian.py, test_math_coalescent.py, test_math_coalescent_expectations.py, test_math_coalescent_extras.py, test_math_comprehensive.py, test_math_decision.py … (26 files).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/math` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.