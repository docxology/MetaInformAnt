# AGENTS.md — `src/metainformant/dna/population`

Sub-package `population` of the `dna` domain module
(verified 2026-08-29: 8 files, 1 subdirs).
Part of the METAINFORMANT package; keep public API in sync with docs and
README. `__pycache__/` is generated.
Repo-wide policy: see the repository-root `AGENTS.md`.

## Rules
- `analysis.py` never imports `metainformant.simulation`. Simulated p-values
  (`p_value_method="simulation"`) go through `dna/integration/neutrality_null.py`, which
  `analysis.py` imports inside the functions that need it, never at module scope, so
  `import metainformant.dna` and the default `"normal"` method load no simulation
  module (pinned by `tests/dna/test_dna_population_simulated_null.py`).
- Keep the neutrality statistic computation shared: the observed statistic and the
  simulated null statistics both go through `core._neutrality_statistics_from_spectra`.