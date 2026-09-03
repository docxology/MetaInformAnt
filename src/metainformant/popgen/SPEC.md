# Specification: popgen

## Scope

Population genetics analysis for METAINFORMANT. Sequence summary
statistics, neutrality tests, population comparison, genotype structure
analysis, linkage disequilibrium summaries, and demographic model
comparisons.

## Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## Data Structures

- **Submodules**: workflow.analysis
- **Key Concepts**: nucleotide diversity, Watterson's theta, Tajima's D,
  Fu & Li D*/F*, Fay & Wu's H, Fst, PCA/kinship/VanRaden, Hardy-Weinberg
  equilibrium, r-squared linkage disequilibrium, bottleneck/expansion
  effective-size models

## API Definition

### Exports — `workflow/analysis.py`

- `summarize_scenario` — Full descriptive suite for one sequence set
- `sequence_scenario_suite` — FASTA-loading variant of summarize_scenario
- `compare_two_population_sequences` — Fst-based two-population comparison
- `genotype_structure_analysis` — PCA + kinship + HWE on a dosage matrix
- `ld_summary` — Adjacent-site LD (r-squared) summary
- `demographic_model_comparisons` — Bottleneck/expansion Ne vs observed diversity
- `analyze_dataset` — Full scenario-dataset analysis pipeline (writes JSON)
- Re-exported generators: `generate_population_sequences`,
  `generate_two_populations`, `simulate_bottleneck_population`,
  `simulate_population_expansion`, `generate_genotype_matrix`,
  `generate_linkage_disequilibrium_data`, `generate_site_frequency_spectrum`
