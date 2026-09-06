# AGENTS.md — `MetaInformAnt/tests/math/`

tests for the math module (population genetics, Bayesian, decision theory…).
Files (verified 2026-09-05): 28 test modules + `__init__.py`, covering population genetics (coalescent, demography, effective size, Fst, LD, selection, statistics, popgen shim), core utilities, Bayesian inference, decision theory, epidemiology, evolutionary dynamics, perception, and quantitative genetics.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/math -p no:cacheprovider` (verified 2026-09-05: 274 passed pre-review).
Repo-wide policy: see the repository-root `AGENTS.md`.