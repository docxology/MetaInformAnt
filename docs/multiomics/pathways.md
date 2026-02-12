# Multi-Omic Pathway Analysis

Methods for combining pathway-level signals across omic layers, detecting
active subnetwork modules, analysing pathway topology effects, and assessing
cross-omic concordance.

## Key Concepts

**P-value combination** merges per-omic significance for each pathway using
one of three strategies:

| Method | Statistic | Distribution |
|--------|-----------|-------------|
| Fisher's combined | -2 sum(log p_i) | Chi-squared (2k df) |
| Stouffer's Z | Weighted sum of Phi^{-1}(1 - p_i) | Standard normal |
| Minimum-p | min(p_i) * n_omics | Bonferroni-corrected |

**Active module detection** uses a greedy seed-and-grow algorithm on a gene
interaction network. Node scores are derived from -log(p) minus a significance
threshold (-log(alpha)). Seeds are nodes with positive transformed scores;
modules expand by adding neighbours that improve the aggregate score.

**Topology-weighted analysis** incorporates graph structure by weighting gene
perturbation scores by degree centrality, producing an impact factor that
reflects both statistical significance and network position.

**Cross-omic concordance** measures whether pathway signals agree across omic
layers using the coefficient of variation (CV) of -log10(p) values. Low CV
indicates concordant signals.

## Function Reference

### multi_omic_enrichment

```python
def multi_omic_enrichment(
    gene_sets: dict[str, list[str]],
    omic_results: dict[str, dict[str, float]],
    method: str = "fisher_combined",
) -> list[dict[str, Any]]
```

Combines per-omic gene-level p-values into pathway-level significance. For
each pathway, Fisher's method is first applied within each omic to obtain
per-omic p-values, which are then combined across omics using the specified
method.

Returns a list sorted by `combined_p`, each with `pathway_id`,
`pathway_name`, `combined_p`, `per_omic_p`, `n_genes`, and `leading_edge`
(top 10 most significant genes).

Methods: `"fisher_combined"`, `"stouffer"`, `"min_p"`.

### active_module_detection

```python
def active_module_detection(
    network: dict[str, list[str]],
    scores: dict[str, float],
    alpha: float = 0.05,
    n_permutations: int = 1000,
) -> list[dict[str, Any]]
```

Finds active subnetwork modules via seed-and-grow. Returns modules sorted by
permutation p-value, each with `module_genes`, `module_score`, `p_value`, and
`omic_contributions`.

### pathway_topology_analysis

```python
def pathway_topology_analysis(
    pathway_graph: dict[str, list[str]],
    gene_scores: dict[str, float],
) -> dict[str, Any]
```

Computes a topology-weighted impact factor:
`IF = sum(score_i * centrality_i) / sum(centrality_i)` where centrality is
normalised degree centrality and score is -log10(p). Significance is assessed
via 1000 permutations.

Returns `impact_factor`, `p_value`, `perturbed_genes` (p < 0.05), and
`pathway_perturbation` (per-gene score, centrality, weighted score).

### cross_omic_pathway_concordance

```python
def cross_omic_pathway_concordance(
    pathway_results: dict[str, dict[str, float]],
) -> dict[str, Any]
```

Assesses signal concordance across >= 2 omic layers. Classifies pathways as
concordant (CV <= median) or discordant (CV > median).

Returns `concordant_pathways`, `discordant_pathways`, global
`concordance_score` (1 - mean CV), and `heatmap_data` for visualisation
(sorted by concordance).

## Usage Example

```python
from metainformant.multiomics.pathways.enrichment import (
    multi_omic_enrichment,
    active_module_detection,
    pathway_topology_analysis,
    cross_omic_pathway_concordance,
)

# Multi-omic enrichment
gene_sets = {"apoptosis": ["BAX", "BCL2", "CASP3"], "cycle": ["CDK2", "CDK4"]}
omic_results = {
    "rna": {"BAX": 0.001, "BCL2": 0.05, "CDK2": 0.8},
    "protein": {"BAX": 0.01, "CASP3": 0.002, "CDK4": 0.1},
}
enriched = multi_omic_enrichment(gene_sets, omic_results, method="fisher_combined")

# Active modules
network = {"BAX": ["BCL2", "CASP3"], "BCL2": ["BAX"], "CASP3": ["BAX"]}
scores = {"BAX": 0.001, "BCL2": 0.05, "CASP3": 0.002}
modules = active_module_detection(network, scores, alpha=0.05)

# Topology analysis
topo = pathway_topology_analysis(network, scores)
print(f"Impact factor: {topo['impact_factor']:.4f}, p={topo['p_value']:.4e}")

# Concordance
pathway_pvals = {"rna": {"apoptosis": 0.01}, "protein": {"apoptosis": 0.02}}
conc = cross_omic_pathway_concordance(pathway_pvals)
print(f"Concordance score: {conc['concordance_score']:.3f}")
```

## Related Modules

- `metainformant.ontology.pathway_enrichment` -- single-omic pathway enrichment
- `metainformant.multiomics.methods` -- factorization and clustering
- `metainformant.multiomics.survival` -- survival analysis
- `metainformant.networks` -- biological network analysis
