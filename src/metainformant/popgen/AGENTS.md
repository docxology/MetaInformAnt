# Agent Directives: src/metainformant/popgen

Module `popgen` of the METAINFORMANT package. Extracted from
`scripts/popgen/` (thin-orchestrator migration, 2026-09-01); also repairs
the defunct `metainformant.dna.population_analysis` and
`metainformant.simulation.popgen` import paths the scripts relied on.

## Role
Reusable population-genetics analysis: scenario statistics suites,
two-population comparison, genotype structure (PCA/kinship/HWE), LD
summaries, demographic model comparisons, and the dataset pipeline
(`analyze_dataset`).

## Contents
- `analysis.py` - all analysis methods + re-exported simulation generators

## Rules
- Scripts under `scripts/popgen/` must only orchestrate; add methods here.
- Core statistics live in `metainformant.dna.population`; this module
  composes them, it does not reimplement them.
- All outputs are descriptive. Label inferential statistics (if ever
  added) as post-evidence-freeze gated.
- Repo-wide policy: see the repository-root `AGENTS.md`.
