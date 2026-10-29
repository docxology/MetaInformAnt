# AGENTS.md — `MetaInformAnt/tests/multiomics/`

tests for the multiomics module.
Files (verified 2026-10-28): __init__.py, test_concordance_median.py, test_multiomics_cca_degenerate.py, test_multiomics_comprehensive.py, test_multiomics_converters_summary.py, test_multiomics_descriptive_summary.py, test_multiomics_integration.py, test_multiomics_methods.py, test_multiomics_methods_hygiene.py, test_multiomics_pathways.py, test_multiomics_pathways_determinism.py, test_multiomics_pathways_fallbacks.py, test_multiomics_pathways_omic_contributions.py, test_multiomics_pathways_random_state.py, test_multiomics_sample_mapping.py, test_multiomics_survival.py, test_multiomics_survival_oracles.py, test_multiomics_visualization_data.py.

- `test_concordance_median.py` — `cross_omic_pathway_concordance` splits concordant and discordant pathways at the
  true median CV (the mean of the two middle values for an even number of pathways).
- `test_multiomics_pathways.py` — multi-omic enrichment, active modules, topology analysis and concordance scoring.
- `test_multiomics_pathways_determinism.py` — `pathway_topology_analysis`, `active_module_detection` and
  `cross_omic_pathway_concordance` return identical results in child processes started with different
  `PYTHONHASHSEED` values, and tied pathways come back in identifier order.
- `test_multiomics_pathways_omic_contributions.py` — `active_module_detection(..., omic_scores=...)` reports each
  module's `omic_contributions` as the module score recomputed from one omic's p-values (`sum(log(alpha / p))`),
  leaves modules and scores unchanged, and gives `{}` when `omic_scores` is omitted.
- `test_multiomics_pathways_random_state.py` — `active_module_detection` and `pathway_topology_analysis` take `random_state`
  (default 42 reproduces the former fixed-seed results exactly) and refuse `n_permutations < 1`, a bad `alpha` and
  out-of-range p-values.
- `test_multiomics_pathways_fallbacks.py` — the scipy-free chi-squared survival function is exact for even degrees of
  freedom, and the `p = 1.0` entry of an unmeasured omic is pinned.
- `test_multiomics_cca_degenerate.py` — unregularized CCA raises `ValueError` where every correlation would be 1
  (a layer of rank `n - 1`) or leading correlations are forced to 1
  (combined ranks greater than `n - 1`), in both `canonical_correlation` functions, which run one shared core and agree; per-feature
  imputation of missing values.
- `test_multiomics_methods_hygiene.py` — no global `random.seed`, no NumPy-version-specific calls (AST guards), a
  private RNG in the pure-Python NMF, the consensus CDF area, and the no-op SNF `alpha` warns.
- `test_multiomics_converters_summary.py` — `from_rna_expression` / `from_metabolomics` choose their log transform from
  the observed scale (zeros and gaps do not switch it off); the descriptive summary keeps real feature names.
- `test_multiomics_survival_oracles.py` — `cox_regression` and `log_rank_test` match statsmodels' `PHReg` (Breslow) and
  `survdiff` (2-5 groups, tied times) and refuse collinear or separated fits.
- `test_multiomics_visualization_data.py` — heatmap row labels, blank (not p = 1) enrichment cells, the dashboard's real
  network layout, seeded and hash-seed-independent plots.
- `test_multiomics_integration.py`, `test_multiomics_sample_mapping.py`, `test_multiomics_descriptive_summary.py` —
  data integration, sample ID mapping and descriptive summaries of the omic layers.
- `test_multiomics_methods.py` — factorization (joint NMF, MOFA, tensor, SNF, CCA) and multi-omic clustering.
- `test_multiomics_survival.py` — Cox models, Kaplan-Meier, log-rank tests, risk stratification, concordance index.
- `test_multiomics_comprehensive.py` — integration workflows across the submodules.

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/multiomics` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.
