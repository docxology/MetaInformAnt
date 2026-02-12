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
| `active_module_detection()` | Greedy seed-and-grow active subnetwork module detection |
| `pathway_topology_analysis()` | Topology-weighted pathway impact analysis |
| `cross_omic_pathway_concordance()` | Assess cross-omic concordance at the pathway level |

## Usage

```python
from metainformant.multiomics.pathways import enrichment

combined = enrichment.multi_omic_enrichment(
    pathway_pvalues={"rna": rna_pvals, "protein": prot_pvals},
    method="fisher"
)
modules = enrichment.active_module_detection(network, node_scores)
topology = enrichment.pathway_topology_analysis(pathway_graph, gene_scores)
concordance = enrichment.cross_omic_pathway_concordance(rna_results, prot_results)
```
