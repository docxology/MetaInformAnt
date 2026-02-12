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
| `over_representation_analysis()` | Hypergeometric test-based over-representation analysis |
| `compute_enrichment_score()` | Running enrichment score computation for GSEA |
| `gsea()` | Full Gene Set Enrichment Analysis with permutation testing |
| `pathway_network()` | Build pathway similarity network from shared gene membership |
| `compare_enrichments()` | Compare enrichment results across conditions or experiments |

## Usage

```python
from metainformant.ontology.pathway_enrichment import enrichment

ora = enrichment.over_representation_analysis(
    gene_list=significant_genes,
    pathways=pathway_db,
    background=all_genes,
)
gsea_result = enrichment.gsea(ranked_genes, pathway_db, n_permutations=1000)
net = enrichment.pathway_network(enrichment_results, overlap_threshold=0.3)
comparison = enrichment.compare_enrichments(results_a, results_b)
```
