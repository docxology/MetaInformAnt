# Specification: pathway_enrichment

## 🎯 Scope
Pathway enrichment analysis subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 2 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `enrichment.py`

### Behavior Contracts

#### `enrichment.py`
- `over_representation_analysis(gene_list, gene_sets, background=None, correction="fdr_bh")` — one-sided hypergeometric
  tail `P(X >= overlap)` summed in log space; the log-factorials are exact (`math.lgamma(n + 1)`). `correction` is
  `"fdr_bh"`, `"bonferroni"` or `"none"`; results are sorted by `adjusted_p`.
- `gsea(ranked_genes, gene_sets, n_permutations=1000, min_size=15, max_size=500, random_state=None)` — preranked GSEA
  (Subramanian et al. 2005). The enrichment score is the maximum deviation of the weighted running sum
  (`|rank metric|` weights for members, `1 / (number of non-members)` for non-members). The null is the score of
  `n_permutations` random gene sets of the same size drawn from the ranked list with the ranking and weights fixed. The
  p-value is `(1 + #(null ES at least as extreme)) / (1 + #(null ES of the observed sign))`, one-sided in the direction
  of the observed ES. NES is the observed ES divided by the mean magnitude of the null ES of the same sign, and is NaN
  when no permutation produced that sign. FDR is Benjamini-Hochberg over the p-values; results are sorted by FDR.
  `random_state` is an integer seed or a `random.Random`; the same seed gives the same results and the global `random`
  state is neither seeded nor consumed.
