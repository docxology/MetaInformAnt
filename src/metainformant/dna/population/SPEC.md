# Specification: population

## 🎯 Scope
Population genetics analysis sub-package.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 6 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `analysis.py`
- `core.py`
- `visualization.py`
- `visualization_core.py`
- `visualization_stats.py`

### Function contracts (behavior pinned by tests)

- `analysis.mcdonald_kreitman_test(sequences) -> (alpha, omega)`: real
  McDonald-Kreitman 2x2 contingency over aligned coding sequences (last
  sequence = outgroup; preceding = ingroup population sample). Sites are
  counted per codon column; ambiguous codons are skipped. `alpha = 1 - NI`
  with `NI = (Pn/Ps) / (Dn/Ds)`; `omega = Dn/Ds`. Fewer than 3 sequences
  returns `(0.0, 0.0)`. See `analysis.mcdonald_kreitman_contingency` for
  the full table (`Ps`, `Pn`, `Ds`, `Dn`, `neutral_ratio`, `alpha`,
  `omega`, `fisher_p` with a scipy Fisher exact two-sided p-value, falling
  back to an exact 2x2 hypergeometric computation).
- `core.linkage_disequilibrium(seqs, pos1, pos2)`: `D = f_AB - p_A * p_B`
  where A and B are the most common alleles at each site with ties broken
  alphabetically -- deterministic across processes (never derived from set
  iteration order). Repulsion haplotypes yield negative D.
- `core.hudson_fst(pop1, pop2)`: Hudson's (1992) moment estimator
  (Bhatia et al. 2013 formulation): `F_ST = sum(num)/sum(den)` with
  per-site sampling corrections; sites with < 2 valid copies per
  population are skipped. Distinct from the heterozygosity-based per-site
  estimator in `analysis.calculate_fst`.
- `analysis.calculate_tajima_d(sequences)`, `calculate_fu_li_d(sequences)` (D*),
  `calculate_fu_li_f(sequences)` (F*), `calculate_fay_wu_h(sequences, outgroup=None)`
  and `neutrality_test_suite(sequences, outgroup=None)` take the keyword-only options
  `p_value_method` (`"normal"`, the default, or `"simulation"`), `n_simulations`
  (a positive integer, default `DEFAULT_NULL_SIMULATIONS` = 2000) and `random_state`
  (an integer, `random.Random`, `numpy.random.Generator` or `None`). An unknown method
  or an invalid `n_simulations` raises `ValueError` before any test runs, in the suite
  as well. The first four return `(statistic, p_value)` and `(0.0, 1.0)` for fewer than
  4 sequences. `"normal"` gives `erfc(abs(z) / sqrt(2))`; `"simulation"` gives the
  two-sided Monte Carlo p-value of `dna.integration.neutrality_null` (constant-size
  neutral coalescent conditioned on the observed number of segregating sites, twice
  the smaller tail with the +1 correction, so never below `2 / (n_simulations + 1)`).
  `calculate_fay_wu_h` returns the raw per-site H with the p-value of Zeng et al.'s
  normalized H, `None` for a polymorphic sample without an outgroup and `1.0` for a
  monomorphic one. The suite adds the keys `p_value_method` and `n_simulations` only
  when the method is `"simulation"`; `calculate_summary_statistics` and
  `detect_selection` always use the normal approximation. `analysis` imports
  `dna.integration.neutrality_null` inside the functions that need it, so
  `import metainformant.dna` and the `"normal"` method load no `metainformant.simulation`
  module. Pinned by `tests/dna/test_dna_population_simulated_null.py`.
