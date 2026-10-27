# Pathway Enrichment

Over-representation analysis (ORA), Gene Set Enrichment Analysis (GSEA), pathway similarity networks, and cross-condition enrichment comparison using pure Python statistical implementations.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `enrichment` module |
| `enrichment.py` | ORA, GSEA, pathway networks, enrichment comparison |

## Key Functions

| Function | Description |
|----------|-------------|
| `over_representation_analysis()` | Hypergeometric test-based over-representation analysis (log-factorials are exact, `lgamma(n + 1)`) |
| `compute_enrichment_score()` | Weighted running enrichment score of one gene set (Subramanian et al. 2005) |
| `gsea()` | Preranked Gene Set Enrichment Analysis with a gene-set permutation null; `random_state` seeds the random gene sets |
| `pathway_network()` | Build pathway similarity network from shared gene membership |
| `compare_enrichments()` | Compare enrichment results across conditions or experiments |

## Usage

```python
from metainformant.ontology.pathway_enrichment import enrichment

# gene_sets: term ID -> list or set of genes; ranked_genes: (gene, rank metric) tuples, sorted by the metric
ora = enrichment.over_representation_analysis(
    gene_list=significant_genes,
    gene_sets=gene_sets,
    background=all_genes,
)
gsea_result = enrichment.gsea(ranked_genes, gene_sets, n_permutations=1000, random_state=0)
net = enrichment.pathway_network(ora, gene_sets, similarity_threshold=0.3)
comparison = enrichment.compare_enrichments(results_a, results_b)
```

## GSEA

`gsea()` implements the preranked gene-set analysis of Subramanian et al. (2005):

- The enrichment score (ES) of a gene set is the maximum deviation from zero of a weighted running sum down the ranked
  list: each member adds `|rank metric|` divided by the sum of the members' `|rank metric|`, and each non-member
  subtracts `1 / (number of non-members)`.
- The null is a gene-set permutation: `n_permutations` random gene sets of the same size are drawn from the ranked
  list, with the ranking and the weights fixed, and scored the same way. `random_state` (an integer seed or a
  `random.Random`) draws them: the same seed gives the same p-values, NES and FDR, and the global `random` state is
  neither seeded nor consumed. `None` draws fresh entropy.
- The nominal p-value is one-sided in the direction of the observed ES, taken against the part of the null with the same
  sign and corrected by +1: `(1 + #(null ES at least as extreme)) / (1 + #(null ES of that sign))`.
- NES is the observed ES divided by the mean magnitude of the null ES of the same sign, and is NaN when no permutation
  produced an ES of that sign.
- FDR is the Benjamini-Hochberg adjustment of the p-values. The result list is sorted by FDR.

The stage-10 workflow passes `ontology.random_state` of its configuration to `gsea()` (see
[`../workflow/`](../workflow/)).
