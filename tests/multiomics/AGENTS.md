# AGENTS.md — `MetaInformAnt/tests/multiomics/`

tests for the multiomics module.
Files (verified 2026-10-27): __init__.py, test_concordance_median.py, test_multiomics_comprehensive.py, test_multiomics_descriptive_summary.py, test_multiomics_integration.py, test_multiomics_methods.py, test_multiomics_pathways.py, test_multiomics_pathways_determinism.py, test_multiomics_pathways_omic_contributions.py, test_multiomics_sample_mapping.py, test_multiomics_survival.py.

- `test_concordance_median.py` — `cross_omic_pathway_concordance` splits concordant and discordant pathways at the
  true median CV (the mean of the two middle values for an even number of pathways).
- `test_multiomics_pathways.py` — multi-omic enrichment, active modules, topology analysis and concordance scoring.
- `test_multiomics_pathways_determinism.py` — `pathway_topology_analysis`, `active_module_detection` and
  `cross_omic_pathway_concordance` return identical results in child processes started with different
  `PYTHONHASHSEED` values, and tied pathways come back in identifier order.
- `test_multiomics_pathways_omic_contributions.py` — `active_module_detection(..., omic_scores=...)` reports each
  module's `omic_contributions` as the module score recomputed from one omic's p-values (`sum(log(alpha / p))`),
  leaves modules and scores unchanged, and gives `{}` when `omic_scores` is omitted.
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
