# Specification: integration

## 🎯 Scope
Cross-domain integration sub-package: DNA-RNA integration (`rna.py`) and the
DNA-simulation adapter that supplies simulated null distributions for the neutrality
tests (`neutrality_null.py`).

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code
- **Boundary**: cross-domain imports are allowed only from adapter surfaces
  (`integration` and `workflow` modules, CLI entry points and the MCP tool layer;
  `scripts/quality/check_module_boundaries.py`). `neutrality_null.py` is the only module
  of `dna` that imports `simulation`.

## 💾 Data Structures
- **Modules**: 3 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py` (imports `rna` only)
- `rna.py`
- `neutrality_null.py` (`FOLDED_TESTS`, `POLARIZED_TESTS`, `monte_carlo_two_sided_p`,
  `simulated_neutrality_tests`)

### Function contracts (behavior pinned by tests)

Pinned by `tests/dna/test_dna_population_simulated_null.py`.

- `neutrality_null.simulated_neutrality_tests(sequences, tests, *, n_simulations,
  random_state=None, outgroup=None) -> {test: (statistic, p_value)}`: the observed
  statistic of each requested test and its two-sided Monte Carlo p-value against the
  standard neutral coalescent (constant size, infinite sites, no recombination) for the
  same number of sequences. Every replicate has exactly as many mutations as the
  observed spectrum has, placed on an unweighted genealogy in proportion to branch
  length (`simulation.models.popgen.simulate_coalescent_spectra(...,
  segregating_sites=...)`, the `ms -s` approach). `FOLDED_TESTS` (`tajima_d`,
  `fu_li_d_star`, `fu_li_f_star`) score the folded spectrum and need no outgroup;
  `POLARIZED_TESTS` (`fay_wu_h`) need `outgroup` and the polarized spectrum. Each group
  is scored against its own simulation. The observed and the replicate statistics go
  through one computation (`core._neutrality_statistics_from_spectra`). A sample without
  mutations returns `(0.0, 1.0)`. An unknown test name, or a polarized test without an
  outgroup, raises `ValueError`.
- `neutrality_null.monte_carlo_two_sided_p(observed, simulated) -> float`: each
  one-sided tail counts the simulated values at least as extreme as `observed` (ties
  included, up to rounding) with the +1 correction, `(1 + count) / (R + 1)`; the
  p-value is twice the smaller tail, at most 1 and at least `2 / (R + 1)`.
- Import isolation: `dna.population.analysis` imports `neutrality_null` inside the
  functions that need it, when `p_value_method="simulation"`; `import
  metainformant.dna` and the default `p_value_method="normal"` do not load any
  `metainformant.simulation` module.
