# Multi-Omic Pathway Analysis

Methods for combining pathway-level signals across multiple omic layers, including enrichment analysis, active module detection, topology analysis, and cross-omic concordance assessment.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `enrichment` module |
| `enrichment.py` | Multi-omic enrichment, active modules, topology analysis, concordance scoring |

## Key Functions

| Function | Description |
|----------|-------------|
| `multi_omic_enrichment()` | Combine pathway p-values across omics via Fisher/Stouffer/min-p |
| `active_module_detection()` | Greedy seed-and-grow active subnetwork module detection; optional `omic_scores={omic: {node: p}}` gives each module an `omic_contributions` entry per omic (`{}` without it); `random_state` (default 42) seeds every module's permutation test, `n_permutations` must be at least 1 and `alpha` in (0, 1] |
| `pathway_topology_analysis()` | Topology-weighted pathway impact analysis |
| `cross_omic_pathway_concordance()` | Assess cross-omic concordance at the pathway level |

## Usage

```python
from metainformant.multiomics.pathways import enrichment

combined = enrichment.multi_omic_enrichment(
    gene_sets={"CELL_CYCLE": ["CDK1", "CDK2", "TP53"]},
    omic_results={"rna": rna_gene_pvalues, "protein": protein_gene_pvalues},
    method="fisher_combined",
)
modules = enrichment.active_module_detection(network, node_scores, omic_scores=omic_node_pvalues)
topology = enrichment.pathway_topology_analysis(pathway_graph, gene_scores)
concordance = enrichment.cross_omic_pathway_concordance(
    {"rna": rna_pathway_pvalues, "protein": protein_pathway_pvalues}
)
```

## Behaviour

- **`omic_contributions`**: for each module, `{omic: score}` is the module score recomputed from that omic's
  p-values alone with the same `-log(p) - (-log(alpha))` transform (a node absent from an omic counts as
  `p = 1`). A layer that supports the module scores above zero and one that does not scores below zero. The
  argument does not change which modules are found or how they are scored.
- **Determinism**: node and pathway identifiers are visited in sorted order, ties are broken by identifier, and
  the permutation tests start from a `random.Random(random_state)` (default 42, also for `pathway_topology_analysis`,
  which takes `n_permutations` as well), so results do not depend on `PYTHONHASHSEED`.

See [docs/multiomics/pathways.md](../../../../docs/multiomics/pathways.md) for the full reference.
