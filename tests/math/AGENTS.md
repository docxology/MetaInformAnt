# AGENTS.md — `MetaInformAnt/tests/math/`

tests for the math module (population genetics, Bayesian, decision theory…).
Files (verified 2026-10-27): 30 test modules + `__init__.py`, covering population genetics (coalescent, demography, effective size, Fst, LD, selection, statistics, popgen shim), core utilities, Bayesian inference, decision theory, epidemiology, evolutionary dynamics, perception, and quantitative genetics.

`test_math_r3_fixes.py` checks these behaviours against independent oracles (closed forms, SciPy and simulation): the argument types and values `effective_size_from_family_size_variance` accepts, the `effective_size_sex_ratio` formula, `fisher_exact_test` raising `ImportError` without SciPy, the `ld_decay_r2` fit, the DDM closed forms, and the `random_state` seeding of `metropolis_hastings` and `simulate_coalescent` plus the global-state isolation of `calculate_confidence_intervals`.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/math -p no:cacheprovider` (verified 2026-10-27: 504 passed).
Repo-wide policy: see the repository-root `AGENTS.md`.