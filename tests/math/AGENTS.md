# AGENTS.md — `MetaInformAnt/tests/math/`

tests for the math module (population genetics, Bayesian, decision theory…).
Files: 32 test modules + `__init__.py`, covering population genetics (coalescent, demography, effective size, Fst, LD, selection, statistics, popgen shim), core utilities, Bayesian inference, decision theory, epidemiology, evolutionary dynamics, perception, and quantitative genetics.

`test_math_r3_fixes.py` checks these behaviours against independent oracles (closed forms, SciPy and simulation): the argument types and values `effective_size_from_family_size_variance` accepts, the `effective_size_sex_ratio` formula, `fisher_exact_test` raising `ImportError` without SciPy, the `ld_decay_r2` fit, the DDM closed forms, and the `random_state` seeding of `metropolis_hastings` and `simulate_coalescent` plus the global-state isolation of `calculate_confidence_intervals`.

`test_math_ddm_wiener_density.py` checks the exact Wiener first-passage-time density and `fit_ddm_parameters` against independent oracles: the scale-function boundary probabilities, the Laplace transform of the absorption time from the backward equation, the closed-form mean decision time, the inverse Gaussian and leading-eigenfunction limits, an Euler-Maruyama simulation (chi-square test), the series truncation bound, and recovery of planted parameters within four observed-information standard errors (plus a 24-dataset replicate study of bias).

`test_math_r4_fixes.py` covers the round-4 defects: the standard `(4N - 2) / (Vk + 2)` in list mode of `effective_size_from_family_size_variance` (hand-computed values), no invented p-values, variances or trial counts (`permutation_test`, `compare_population_statistic`, `conjugate_normal`, `sdt_metrics`), exact Beta credible intervals, `compute_dic`, `metropolis_hastings` guards, input types that raised a bare `AssertionError`, `hardy_weinberg_genotype_freqs`, the replicator derivative at zero mean fitness, and the coalescent expectations (E[pi] = theta against a Kingman simulation, exact harmonic numbers, effective-size validation).

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/math -p no:cacheprovider` (passes: 2026-10-28).
Repo-wide policy: see the repository-root `AGENTS.md`.