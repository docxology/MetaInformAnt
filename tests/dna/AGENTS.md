# AGENTS.md — `MetaInformAnt/tests/dna/`

tests for DNA/sequence/alignment/phylogeny modules.
Files (verified 2026-10-27): __init__.py plus 53 test modules (test_dna_accession.py, test_dna_alignment.py, test_dna_annotation.py, test_dna_codon_usage.py, test_dna_compatibility_facades.py, test_dna_comprehensive.py, test_dna_consensus.py, …, test_dna_variation_new.py).


## Children (documented on disk)
- data/

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Statistical checks use seeded samples and tolerances measured in standard errors, so
  they are deterministic.
- `test_dna_population_simulated_null.py` covers the simulated null distributions of the
  neutrality tests (`p_value_method="simulation"` of `dna/population/analysis.py` and the
  adapter `dna/integration/neutrality_null.py`): p-value uniformity on neutral coalescent
  samples (drawn with `metainformant.simulation.models.popgen`), small p-values after strong
  expansions and bottlenecks, the Monte Carlo arithmetic, an exactly integrated conditional
  null for 4 sequences, option validation, and a subprocess check that the default
  `"normal"` method loads no `metainformant.simulation` module.
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/dna` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.
