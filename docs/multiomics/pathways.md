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
modules expand by adding neighbours that improve the aggregate score. The
combined node p-values carry no per-omic information, so the optional
`omic_scores` argument supplies per-omic p-values from which each module's
`omic_contributions` are computed.

**Topology-weighted analysis** incorporates graph structure by weighting gene
perturbation scores by degree centrality, producing an impact factor that
reflects both statistical significance and network position.

**Cross-omic concordance** measures whether pathway signals agree across omic
layers using the coefficient of variation (CV) of -log10(p) values. Low CV
indicates concordant signals.

**Determinism.** The pathway analyses return the same result for the same
input in every process. Sets of node and pathway identifiers are visited in
sorted order (Python randomises string hashes per process, so set order would
otherwise vary), ties are broken by identifier, and the two permutation tests
(`active_module_detection`, `pathway_topology_analysis`) draw from
`random.Random(42)`. Results do not depend on `PYTHONHASHSEED`.

## Function Reference

### multi_omic_enrichment

```python-snippet
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

```python-snippet
def active_module_detection(
    network: dict[str, list[str]],
    scores: dict[str, float],
    alpha: float = 0.05,
    n_permutations: int = 1000,
    omic_scores: dict[str, dict[str, float]] | None = None,
) -> list[dict[str, Any]]
```

Finds active subnetwork modules via seed-and-grow. `scores` holds one combined
p-value per node. The module score is the sum of the node scores
`-log(p) - threshold` with `threshold = -log(alpha)`; a node missing from
`scores` counts as `p = 1`. Ties between equally scored neighbours are broken
by node identifier, and the permutation test uses a fixed seed, so the same
input always returns the same modules.

Returns modules sorted by permutation p-value, each with:

- `module_genes`: the sorted node identifiers of the module
- `module_score`: the aggregate score of the module
- `p_value`: the permutation p-value
- `omic_contributions`: `{omic: score}`, the module score recomputed from that
  omic's p-values alone with the same `-log(p) - threshold` transform (a node
  absent from an omic counts as `p = 1`, as in `scores`). It shows which omic
  layers support the module: a layer that does not support it scores below
  zero. It is an empty dict `{}` when `omic_scores` is not given.

`omic_scores` is the per-omic mapping `{omic: {node: p_value}}`, the same shape
that `multi_omic_enrichment` takes as `omic_results`. It does not change which
modules are found or how they are scored; it only supplies the data for
`omic_contributions`.

### pathway_topology_analysis

```python-snippet
def pathway_topology_analysis(
    pathway_graph: dict[str, list[str]],
    gene_scores: dict[str, float],
) -> dict[str, Any]
```

Computes a topology-weighted impact factor:
`IF = sum(score_i * centrality_i) / sum(centrality_i)` where centrality is
normalised degree centrality and score is -log10(p). Significance is assessed
via 1000 permutations of the gene scores over the nodes in sorted order, drawn
from a fixed seed, so the p-value is reproducible.

Returns `impact_factor`, `p_value`, `perturbed_genes` (p < 0.05, sorted), and
`pathway_perturbation` (per-gene score, centrality, weighted score).

### cross_omic_pathway_concordance

```python-snippet
def cross_omic_pathway_concordance(
    pathway_results: dict[str, dict[str, float]],
) -> dict[str, Any]
```

Assesses signal concordance across >= 2 omic layers. Classifies pathways as
concordant (CV <= median CV, the true median: the mean of the two middle values
when the number of pathways is even) or discordant (CV > median).

Returns `concordant_pathways`, `discordant_pathways` (both sorted),
global `concordance_score` (1 - mean CV), and `heatmap_data` for
visualisation, ordered by concordance (low CV first) with ties in pathway
identifier order.

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

# Active modules; omic_scores adds the per-omic contribution of each module
network = {"BAX": ["BCL2", "CASP3"], "BCL2": ["BAX"], "CASP3": ["BAX"]}
scores = {"BAX": 0.001, "BCL2": 0.05, "CASP3": 0.002}
modules = active_module_detection(network, scores, alpha=0.05, omic_scores=omic_results)
print(modules[0]["module_genes"])          # ['BAX', 'CASP3']
print(modules[0]["omic_contributions"])    # {'rna': 0.916..., 'protein': 4.828...}
# protein supports the module more than rna, which measured only BAX (CASP3 counts as p = 1 there)

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
