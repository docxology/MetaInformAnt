# Networks: Pathway Analysis

The pathway analysis module provides tools for biological pathway enrichment analysis and pathway network construction.

## Pathway Network Construction

### PathwayNetwork Class

Container for biological pathways with gene sets and optional per-pathway metadata:

```python
from metainformant.networks.analysis.pathway import PathwayNetwork

# Create pathway network
pathway_network = PathwayNetwork(name="demo_pathways")

# Add pathways from database
pathways = {
    "KEGG:00010": {
        "name": "Glycolysis / Gluconeogenesis",
        "genes": ["HK1", "PFKL", "FBA", "PGK1", "PGAM1"],
        "description": "Central carbohydrate metabolism pathway"
    },
    "KEGG:00020": {
        "name": "Citrate cycle (TCA cycle)",
        "genes": ["CS", "ACO1", "IDH1", "OGDH", "SUCLG1"],
        "description": "Tricarboxylic acid cycle"
    }
}

for pathway_id, info in pathways.items():
    pathway_network.add_pathway(
        pathway_id,
        info["genes"],
        metadata={"name": info["name"], "description": info["description"]},
    )

print(f"Network contains {len(pathway_network)} pathways")
print(pathway_network.get_statistics())  # num_pathways, num_genes, avg/max/min_pathway_size
```

`network.pathways` maps each pathway ID to its gene list, and `network.pathway_metadata` holds the metadata passed to `add_pathway`. Other helpers include `get_pathway_genes(pathway_id)`, `get_gene_pathways(gene)`, `find_pathways_containing_gene(gene)`, `pathway_similarity(pathway1, pathway2)` (Jaccard), `find_overlapping_pathways(pathway_id)`, `filter_pathways_by_size(min_size=5, max_size=None)` and `pathway_overlap_matrix()`.

### load_pathway_database()

Build a `PathwayNetwork` from structured pathway data that you have already parsed from a database export (KEGG, Reactome, GO, WikiPathways, ...). The function does not query any database itself:

```python
from metainformant.networks.analysis.pathway import load_pathway_database

# Structured data: {"pathways": {pathway_id: {"name": ..., "genes": [...]}}, ...}
pathway_data = {
    "pathways": pathways,    # the dictionary defined above
    "source": "KEGG",        # stored in network.metadata["source"]
}
kegg_network = load_pathway_database(pathway_data, name="kegg_subset")

# A pathway entry may also be a plain gene list
simple_network = load_pathway_database({"pathway_a": ["G1", "G2", "G3"], "pathway_b": ["G3", "G4"]})

# Keep pathways with 5 to 200 genes
kegg_network = kegg_network.filter_pathways_by_size(min_size=5, max_size=200)
print(len(kegg_network), kegg_network.metadata)
```

Per-pathway keys other than `genes` (such as `name`) end up in `network.pathway_metadata`. `PathwayNetwork.load_from_database(path, format="json")` reads a file instead: JSON with a top-level `"pathways"` mapping of pathway ID to gene list, or `format="gmt"` for GMT gene-set files.

## Enrichment Analysis

### pathway_enrichment()

Perform pathway enrichment analysis:

```python
from metainformant.networks.analysis.pathway import pathway_enrichment

# Gene list of interest (e.g., differentially expressed genes)
gene_list = ["HK1", "PFKL", "FBA", "PGK1", "CS", "ACO1"]

# Perform enrichment analysis
enriched_pathways = pathway_enrichment(
    gene_list,
    pathway_network,
    background_genes=None,    # All genes in the network, or specify background
    method="hypergeometric",  # Statistical method
    correction="bonferroni",  # Multiple testing correction
    min_overlap=2             # Minimum genes in common
)

# Display results: a dict keyed by pathway ID, best p-values first
ranked = sorted(enriched_pathways.items(), key=lambda item: item[1]["p_value"])
for pathway_id, result in ranked[:10]:  # Top 10 enriched pathways
    name = pathway_network.pathway_metadata.get(pathway_id, {}).get("name", pathway_id)
    print(f"{pathway_id}: {name}")
    print(f"  p-value: {result['p_value']:.2e}")
    print(f"  Genes: {result['overlap_genes']}")
    print(f"  Enrichment: {result['enrichment_ratio']:.2f}")
```

Each result holds `overlap_size`, `pathway_size`, `query_size`, `background_size`, `p_value`, `enrichment_ratio`, `expected_overlap` and the overlapping genes (`overlap_genes`, also available as `overlapping_genes`). When a `correction` is requested, every result also carries the adjusted value as `corrected_p_value`.

**Available methods** (any other name raises `ValueError`):
- `'fisher'`: One-sided Fisher's exact test (default)
- `'hypergeometric'`: Hypergeometric test

Both methods report the one-sided enrichment p-value `P(X >= overlap)`, so they test the same upper tail. With scipy installed the p-value comes from `scipy.stats`; without scipy, `pathway_enrichment` computes the same upper tail exactly in integer arithmetic, so the result does not depend on whether scipy is available. `correction` is `'bonferroni'` (default), `'fdr'` (Benjamini-Hochberg) or `None`. Benjamini-Hochberg is pure Python, so `'fdr'` is applied with or without scipy. An unrecognised `correction` value copies the raw p-values into `corrected_p_value` unchanged.

### network_enrichment_analysis()

Network-based enrichment analysis. The function takes the same arguments as `pathway_enrichment` and simply forwards to it, so the network neighbourhood is whatever gene list you pass in. To enrich a 1-hop neighbourhood, expand the list from an interaction network first:

```python
from metainformant.networks.analysis.pathway import network_enrichment_analysis
from metainformant.networks.interaction.ppi import ProteinNetwork

# Protein interaction network
ppi_network = ProteinNetwork(name="demo_ppi")
ppi_network.add_interaction("HK1", "PFKL", 0.9)
ppi_network.add_interaction("PFKL", "FBA", 0.8)
ppi_network.add_interaction("TP53", "MDM2", 0.95)

# Expand the genes of interest by their direct partners (1-hop neighborhood)
gene_list = ["HK1", "TP53"]
neighborhood = set(gene_list)
for gene in gene_list:
    neighborhood.update(ppi_network.get_protein_partners(gene))

# Network-based enrichment
network_enriched = network_enrichment_analysis(
    sorted(neighborhood),
    pathway_network,
    method="fisher",
    correction="fdr",
)

ranked = sorted(network_enriched.items(), key=lambda item: item[1]["p_value"])
for pathway_id, result in ranked[:5]:
    print(f"{pathway_id}: p={result['p_value']:.2e}")
    print(f"  Overlapping genes: {result['overlap_genes']}")
```

## Pathway Similarity

### pathway_similarity_analysis() and create_pathway_network()

`pathway_similarity_analysis(pathways, method="jaccard")` scores every unordered pair of pathways, and `create_pathway_network(pathways, similarity_threshold=0.1, method="jaccard")` turns those scores into a `networkx` graph. Both take a mapping of pathway ID to a gene list (or to a dict with a `genes` entry, as in `network.pathways`); `method` is `'jaccard'`, `'overlap'` or `'cosine'`, and any other value raises `ValueError`:

```python
from metainformant.networks.analysis.pathway import create_pathway_network, pathway_similarity_analysis

# Pathway ID -> gene list; the IDs may contain underscores
pathway_sets = {
    "KEGG_GLYCOLYSIS": ["HK1", "PFKL", "FBA", "PGK1", "PGAM1"],
    "KEGG_GLUCONEOGENESIS": ["FBA", "PGK1", "PGAM1", "PCK1"],
    "KEGG_TCA_CYCLE": ["CS", "ACO1", "IDH1", "OGDH", "SUCLG1"],
    "KEGG_UNANNOTATED": [],  # no genes: similarity 0.0 with every other pathway
}

similarity = pathway_similarity_analysis(pathway_sets, method="jaccard")
print(similarity["n_pathways"], similarity["pathway_names"])

# One node per pathway; an edge joins two pathways whose similarity reaches the threshold
graph = create_pathway_network(pathway_sets, similarity_threshold=0.2, method="jaccard")
print(sorted(graph.edges(data="similarity")))  # [('KEGG_GLYCOLYSIS', 'KEGG_GLUCONEOGENESIS', 0.5)]
```

`pathway_similarity_analysis` returns `similarities`, `method`, `n_pathways` and `pathway_names`. The keys of `similarities` have the form `"{name1}_{name2}"`. A key cannot be split back into the two pathway names when a pathway ID contains an underscore, and two different pairs can then share a key, so do not parse the keys. `create_pathway_network` does not: it scores the pairs directly, so pathway IDs with underscores (`KEGG_GLYCOLYSIS`, `HALLMARK_APOPTOSIS`, ...) produce edges between the right nodes. Each edge stores the pair's score as `similarity`, and a dict-valued pathway entry becomes the attributes of its node.

A pathway with no genes (an empty list, a dict with an empty `genes` entry, or an entry that is neither a list nor a dict with `genes`) has similarity 0.0 with every other pathway under all three methods. The cosine of an empty gene set is undefined and is defined as 0.0, like the Jaccard and overlap scores.

## Pathway Integration

### pathway_activity_score()

Summarize pathway activity from expression data (aggregated expression of the pathway's genes):

```python-snippet
import pandas as pd
from metainformant.networks.analysis.pathway import pathway_activity_score

# Expression of one sample as a gene -> value mapping
expression_data = pd.read_csv("expression_matrix.csv", index_col=0)
sample_expression = expression_data["sample1"].to_dict()

# Activity of every pathway in the network
pathway_activities = {
    pathway_id: pathway_activity_score(
        pathway_network,
        pathway_id,
        sample_expression,
        method="mean"           # 'mean', 'max', 'sum' or 'median'
    )
    for pathway_id in pathway_network.pathways
}

# Pathway activities as new features
activity_df = pd.DataFrame.from_dict(pathway_activities, orient="index", columns=["activity"])
activity_df.to_csv("pathway_activities.csv")
```

## Pathway Visualization

### Pathway Heatmap

Visualize pathway enrichment results:

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Create enrichment heatmap
top_pathways = sorted(enriched_pathways.items(), key=lambda item: item[1]["p_value"])[:20]  # Top 20 pathways
all_genes = sorted({gene for genes in pathway_network.pathways.values() for gene in genes})[:100]  # First 100 genes

pathway_matrix = []
for pathway_id, result in top_pathways:
    # Create binary vector for pathway genes
    pathway_vector = [1 if gene in result['overlap_genes'] else 0
                     for gene in all_genes]
    pathway_matrix.append(pathway_vector)

# Plot heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(
    pathway_matrix,
    xticklabels=all_genes,
    yticklabels=[pathway_id[:30] for pathway_id, _ in top_pathways],
    cmap='RdBu_r',
    center=0
)

plt.title('Pathway Enrichment Heatmap')
plt.xlabel('Genes')
plt.ylabel('Pathways')
plt.show()
```

### Pathway Network Visualization

```python
import matplotlib.pyplot as plt
import networkx as nx

# Visualize pathway relationships
fig, ax = plt.subplots(figsize=(14, 10))

# Pathway similarity: Jaccard similarity between the gene sets of every pathway pair
pathway_similarities = pathway_network.pathway_overlap_matrix()

# Draw pathway network
G = nx.Graph()
G.add_nodes_from(pathway_network.pathways)
for pathway1, similarities in pathway_similarities.items():
    for pathway2, similarity in similarities.items():
        if pathway1 < pathway2 and similarity > 0.2:
            G.add_edge(pathway1, pathway2, weight=similarity)

pos = nx.spring_layout(G, k=3, iterations=100)
nx.draw(
    G,
    pos=pos,
    ax=ax,
    node_size=1000,
    node_color='lightblue',
    edge_color='gray',
    with_labels=True,
    font_size=8,
    alpha=0.7
)

plt.title('Pathway Similarity Network')
plt.show()
```

## Integration Examples

### With Expression Data

```python-snippet
import pandas as pd
from metainformant.rna.analysis import expression
from metainformant.networks.analysis.pathway import pathway_enrichment

# counts: raw count matrix with genes as rows and samples as columns
counts = pd.read_csv("counts.csv", index_col=0)
conditions = ["control", "control", "control", "treated", "treated", "treated"]

# Find differentially expressed genes
de_results = expression.differential_expression(counts, conditions)
diff_genes = de_results.loc[de_results["adjusted_p_value"] < 0.05, "gene"].tolist()

# Pathway enrichment
enriched = pathway_enrichment(
    diff_genes,
    pathway_network,
    background_genes=counts.index.tolist(),
    correction="fdr"
)

# Export results
enrichment_df = pd.DataFrame([
    {
        'pathway_id': pathway_id,
        'pathway_name': pathway_network.pathway_metadata.get(pathway_id, {}).get('name', pathway_id),
        'p_value': result['p_value'],
        'corrected_p_value': result.get('corrected_p_value', result['p_value']),
        'n_overlap': result['overlap_size'],
        'overlap_genes': ','.join(result['overlap_genes'])
    }
    for pathway_id, result in enriched.items()
])

enrichment_df.to_csv("pathway_enrichment_results.csv", index=False)
```

### With Protein Data

```python-snippet
import pandas as pd
from metainformant.networks.analysis.pathway import pathway_enrichment

# Protein abundance table (indexed by gene symbol, with a p_value column)
protein_abundances = pd.read_csv("proteomics_data.csv", index_col=0)

# Find differentially abundant proteins
diff_proteins = protein_abundances[
    protein_abundances['p_value'] < 0.01
].index.tolist()

# Pathway enrichment for proteins
protein_enriched = pathway_enrichment(
    diff_proteins,
    pathway_network,
    method="fisher",
    correction="bonferroni"
)

# Compare with transcriptomic enrichment (diff_genes comes from the expression example)
transcriptomic_enriched = pathway_enrichment(diff_genes, pathway_network)
shared_pathways = set(protein_enriched) & set(transcriptomic_enriched)
```

## Testing

Pathway analysis functionality is tested in `tests/networks/test_networks_pathway.py`:

```bash
# Run pathway analysis tests
uv run pytest tests/networks/test_networks_pathway.py -v

# Test specific classes and functions
uv run pytest tests/networks/test_networks_pathway.py::TestPathwayEnrichment::test_pathway_enrichment_basic -v
uv run pytest tests/networks/test_networks_pathway.py::TestLoadPathwayDatabase::test_load_from_dict -v
```

## Related Documentation

- [Networks Overview](./index.md): Main networks module documentation
- [Community Detection](./community.md): Community detection in pathway networks
- [Graph Analysis](./graph.md): Basic graph operations and metrics
- [PPI Networks](./ppi.md): Protein-protein interaction analysis
