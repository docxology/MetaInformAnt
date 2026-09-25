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
