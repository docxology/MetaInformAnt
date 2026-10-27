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

### Scenario models (the re-exported generators, defined in `simulation.models.popgen`)

- `generate_population_sequences`, and `simulate_bottleneck_population` and
  `simulate_population_expansion` with `n_sequences`: samples of the neutral Kingman
  coalescent with infinite-sites mutation and no recombination, one genealogy per call
  (constant size; a sample drawn at the end of a bottleneck; or exponential growth
  ending at the sample). The statistics of one call are one draw; they scatter around
  the model's expectations.
- `generate_two_populations`: Balding-Nichols allele frequencies at a target Fst.
- `generate_genotype_matrix`: independent genotypes at per-site allele frequencies.
- `generate_linkage_disequilibrium_data`: a first-order Markov chain along the sites;
  adjacent-site r-squared equals `r_squared_target`.
- `generate_site_frequency_spectrum`: independent draws from the neutral 1/i
  expectation (`"constant"`) or heuristic shapes (`"expansion"`, `"bottleneck"`).
- `summarize_scenario` reports the default (normal-approximation) neutrality p-values.
