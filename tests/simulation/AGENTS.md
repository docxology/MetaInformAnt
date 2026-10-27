# AGENTS.md — `MetaInformAnt/tests/simulation/`

tests for the simulation module.
Files (verified 2026-10-27): __init__.py, test_simulation.py, test_simulation_agents.py, test_simulation_benchmark.py, test_simulation_coalescent.py, test_simulation_methylation.py, test_simulation_popgen.py, test_simulation_r3_fixes.py, test_simulation_rng.py, test_simulation_rna_advanced.py, test_simulation_workflow.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Statistical checks use seeded samples and tolerances measured in standard errors, so
  they are deterministic.
- `test_simulation_coalescent.py` validates the coalescent engine of
  `simulation/models/popgen.py` (`PopulationEpoch`, `simulate_coalescent_spectra`,
  `simulate_coalescent_haplotypes`, `haplotypes_to_sequences`) against closed forms
  (Watterson 1975, Tajima 1989, Hudson 1990) and numerical integration.
- `test_simulation_popgen.py` pins the population-genetics generators: the sequence
  generators (neutral, bottleneck, expansion) over 150 seeded replicates, the
  mutation-parameter precedence and the infinite-sites capacity checks, plus the Fst,
  genotype-matrix, site-frequency-spectrum and admixture generators and the statistics
  forms of the bottleneck and expansion functions.
- `test_simulation_r3_fixes.py` holds regression tests for signed fold changes in
  `simulate_differential_expression` and for the `generate_linkage_disequilibrium_data`
  contract (adjacent r-squared, exact allele frequencies, the attainable cap).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/simulation` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.
