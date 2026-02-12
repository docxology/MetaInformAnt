# Pathway Enrichment Analysis

Over-representation analysis (ORA), Gene Set Enrichment Analysis (GSEA),
pathway similarity networks, and cross-condition enrichment comparison. All
statistical computations use pure Python with no external dependencies.

## Key Concepts

**Over-representation analysis** tests whether a gene list is enriched in any
gene set using the hypergeometric test (Fisher's exact test). Multiple testing
is corrected with Benjamini-Hochberg FDR or Bonferroni.

**GSEA** computes a running enrichment score against a ranked gene list.
Significance is estimated via permutation of rank weights, and normalised
enrichment scores (NES) enable comparison across gene sets.

**Pathway similarity networks** connect gene sets that share genes above a
Jaccard similarity threshold. Connected components identify clusters of
related pathways.

## Function Reference

### over_representation_analysis

```python
def over_representation_analysis(
    gene_list: list[str],
    gene_sets: dict,
    background: list[str] | None = None,
    correction: str = "fdr_bh",
) -> list[dict]
```

Tests each gene set for over-representation of `gene_list` genes. If
`background` is None, the union of all gene set members plus the query genes
is used. Returns a list of result dicts sorted by adjusted p-value, each with:
`term_id`, `term_name`, `p_value`, `adjusted_p`, `odds_ratio`, `n_genes`,
`n_overlap`, `overlap_genes`.

Correction options: `"fdr_bh"` (default), `"bonferroni"`, `"none"`.

### compute_enrichment_score

```python
def compute_enrichment_score(
    ranked_list: list[str],
    gene_set: set,
    weighted: bool = True,
    weights: list[float] | None = None,
) -> dict
```

Computes the running enrichment score for a single gene set. Returns `es`
(maximum deviation), `running_es` (list), `hit_indices`, and
`leading_edge_index`.

### gsea

```python
def gsea(
    ranked_genes: list[tuple],
    gene_sets: dict,
    n_permutations: int = 1000,
    min_size: int = 15,
    max_size: int = 500,
) -> list[dict]
```

Full GSEA pipeline. `ranked_genes` is a list of `(gene_name, rank_metric)`
tuples sorted by metric. Gene sets outside the `[min_size, max_size]` range
are skipped. Returns per-gene-set results with `es`, `nes` (normalised),
`p_value`, `fdr`, `leading_edge`, and `leading_edge_genes`.

### pathway_network

```python
def pathway_network(
    enrichment_results: list[dict],
    gene_sets: dict,
    similarity_threshold: float = 0.3,
) -> dict
```

Builds a Jaccard similarity network from enrichment results. Returns `nodes`
(with p-values and gene counts), `edges` (with source, target, jaccard,
n_shared), and `clusters` (connected components via union-find).

### compare_enrichments

```python
def compare_enrichments(
    results_a: list[dict],
    results_b: list[dict],
) -> dict
```

Compares two sets of enrichment results (e.g. two conditions or timepoints).
Identifies `shared_terms`, `unique_a`, `unique_b`, computes direction
`concordance` (fraction with matching effect direction), and produces a full
`comparison_table`.

Significance threshold: adjusted p < 0.05.

## Internal Helpers

| Function | Purpose |
|----------|---------|
| `_hypergeometric_sf` | Survival function for the hypergeometric distribution |
| `_fdr_correction` | Benjamini-Hochberg FDR correction |
| `_bonferroni_correction` | Bonferroni correction |
| `_log_factorial` / `_log_comb` | Log-space combinatorics (Stirling for large n) |

## Usage Example

```python
from metainformant.ontology.pathway_enrichment.enrichment import (
    over_representation_analysis,
    gsea,
    pathway_network,
    compare_enrichments,
)

# ORA
gene_sets = {"apoptosis": ["BAX", "BCL2", "CASP3"], "proliferation": ["MYC", "CDK2"]}
results = over_representation_analysis(["BAX", "CASP3", "MYC"], gene_sets)
for r in results:
    print(f"{r['term_id']}: p={r['adjusted_p']:.4f}, overlap={r['n_overlap']}")

# GSEA
ranked = [("BAX", 3.2), ("MYC", 2.1), ("CDK2", 0.5), ("BCL2", -1.0)]
gsea_results = gsea(ranked, gene_sets, n_permutations=500, min_size=1)

# Pathway network
net = pathway_network(results, gene_sets, similarity_threshold=0.2)
print(f"Clusters: {len(net['clusters'])}")

# Compare conditions
comparison = compare_enrichments(results_a, results_b)
print(f"Concordance: {comparison['concordance']:.2f}")
```

## Related Modules

- `metainformant.ontology.query.query` -- ontology graph traversal
- `metainformant.ontology.core.go` -- Gene Ontology utilities
- `metainformant.multiomics.pathways` -- multi-omic pathway enrichment
