# Networks: Pathway Analysis

The pathway analysis module provides tools for biological pathway enrichment analysis and pathway network construction.

## Pathway Network Construction

### PathwayNetwork Class

Container for biological pathways with gene sets and interactions:

```python
from metainformant.networks import PathwayNetwork

# Create pathway network
pathway_network = PathwayNetwork()

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

pathway_network.add_pathways(pathways)
print(f"Network contains {pathway_network.n_pathways} pathways")
```

### load_pathway_database()

Load pathways from standard biological databases:

```python
from metainformant.networks import load_pathway_database

# Load KEGG pathways
kegg_pathways = load_pathway_database(
    database="kegg",
    organism="hsa",           # Human
    pathway_ids=None,         # All pathways, or specify list
    min_genes=5,             # Minimum genes per pathway
    max_genes=200            # Maximum genes per pathway
)

# Load Reactome pathways
reactome_pathways = load_pathway_database(
    database="reactome",
    organism="Homo sapiens",
    pathway_types=["metabolic", "signaling"]
)

# Load GO biological processes
go_pathways = load_pathway_database(
    database="go",
    go_domain="biological_process",
    min_genes=10
)
```

**Available databases:**
- `'kegg'`: KEGG pathway database
- `'reactome'`: Reactome pathway database
- `'go'`: Gene Ontology (biological_process, molecular_function, cellular_component)
- `'wikipathways'`: WikiPathways database

## Enrichment Analysis

### pathway_enrichment()

Perform pathway enrichment analysis:

```python
from metainformant.networks import pathway_enrichment

# Gene list of interest (e.g., differentially expressed genes)
gene_list = ["HK1", "PFKL", "FBA", "PGK1", "CS", "ACO1"]

# Perform enrichment analysis
enriched_pathways = pathway_enrichment(
    gene_list,
    pathway_network,
    background_genes=None,    # All genes in network, or specify background
    method="hypergeometric",  # Statistical method
    correction="bonferroni",  # Multiple testing correction
    min_overlap=3            # Minimum genes in common
)

# Display results
for pathway in enriched_pathways[:10]:  # Top 10 enriched pathways
    print(f"{pathway['id']}: {pathway['name']}")
    print(f"  p-value: {pathway['p_value']:.2e}")
    print(f"  Genes: {pathway['overlap_genes']}")
    print(f"  Enrichment: {pathway['enrichment_ratio']:.2f}")
```

**Available methods:**
- `'hypergeometric'`: Hypergeometric test (default)
- `'fisher'`: Fisher's exact test
- `'binomial'`: Binomial test
- `'chi_square'`: Chi-square test

### network_enrichment_analysis()

Network-based enrichment analysis:

```python
from metainformant.networks import network_enrichment_analysis

# Analyze enrichment in network neighborhoods
gene_list = ["TP53", "MDM2", "CDKN1A"]
network = load_ppi_network()  # Protein interaction network

# Network-based enrichment
network_enriched = network_enrichment_analysis(
    gene_list,
    network,
    pathway_network,
    neighborhood_size=1,     # 1-hop neighborhood
    method="subgraph"        # Enrichment method
)

for pathway in network_enriched[:5]:
    print(f"{pathway['id']}: p={pathway['p_value']:.2e}")
    print(f"  Network genes: {pathway['network_genes']}")
```

## Pathway Integration

### multi_omics_pathway_analysis()

Integrate pathway analysis across multiple omics layers:

```python
from metainformant.multiomics import integration
from metainformant.networks import multi_omics_pathway_analysis

# Load multi-omics data
multiomics_data = integration.load_multiomics_data([
    "genomics_data.csv",
    "transcriptomics_data.csv",
    "proteomics_data.csv"
])

# Pathway analysis across omics
pathway_results = multi_omics_pathway_analysis(
    multiomics_data,
    pathway_network,
    analysis_type="integrated",  # Analysis type
    min_consensus=2             # Minimum omics layers showing enrichment
)

# Results show pathways enriched in multiple omics layers
for result in pathway_results:
    print(f"Pathway: {result['pathway_name']}")
    print(f"  Enriched in: {result['omics_layers']}")
    print(f"  Combined p-value: {result['combined_p_value']:.2e}")
```

### pathway_activity_inference()

Infer pathway activity from expression data:

```python
from metainformant.networks import pathway_activity_inference

# Infer pathway activities
expression_data = pd.read_csv("expression_matrix.csv", index_col=0)

pathway_activities = pathway_activity_inference(
    expression_data,
    pathway_network,
    method="gsva",           # Gene Set Variation Analysis
    min_genes=5,            # Minimum genes per pathway
    max_genes=500           # Maximum genes per pathway
)

# Pathway activities as new features
activity_df = pd.DataFrame(pathway_activities)
activity_df.to_csv("pathway_activities.csv")
```

## Pathway Visualization

### Pathway Heatmap

Visualize pathway enrichment results:

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Create enrichment heatmap
pathway_matrix = []
for pathway in enriched_pathways[:20]:  # Top 20 pathways
    # Create binary vector for pathway genes
    pathway_vector = [1 if gene in pathway['overlap_genes'] else 0
                     for gene in all_genes[:100]]  # First 100 genes
    pathway_matrix.append(pathway_vector)

# Plot heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(
    pathway_matrix,
    xticklabels=[f"Gene_{i}" for i in range(100)],
    yticklabels=[p['name'][:30] for p in enriched_pathways[:20]],
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
# Visualize pathway relationships
fig, ax = plt.subplots(figsize=(14, 10))

# Create pathway similarity network
pathway_similarities = calculate_pathway_similarities(enriched_pathways)

# Draw pathway network
G = nx.Graph()
for i, pathway1 in enumerate(enriched_pathways[:15]):
    for j, pathway2 in enumerate(enriched_pathways[:15]):
        if i != j and pathway_similarities[i,j] > 0.3:
            G.add_edge(pathway1['name'][:20], pathway2['name'][:20],
                      weight=pathway_similarities[i,j])

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

```python
from metainformant.rna import workflow
from metainformant.networks import pathway_enrichment

# Load expression data
expression_data = workflow.extract_expression_patterns(rna_data)

# Find differentially expressed genes
diff_genes = expression_data[
    expression_data['p_value'] < 0.05
].index.tolist()

# Pathway enrichment
enriched = pathway_enrichment(
    diff_genes,
    pathway_network,
    background_genes=expression_data.index.tolist(),
    correction="fdr"
)

# Export results
enrichment_df = pd.DataFrame([
    {
        'pathway_id': p['id'],
        'pathway_name': p['name'],
        'p_value': p['p_value'],
        'q_value': p['q_value'],
        'n_overlap': len(p['overlap_genes']),
        'overlap_genes': ','.join(p['overlap_genes'])
    }
    for p in enriched
])

enrichment_df.to_csv("pathway_enrichment_results.csv", index=False)
```

### With Protein Data

```python
from metainformant.protein import proteomes
from metainformant.networks import pathway_enrichment

# Load protein data
protein_abundances = proteomes.load_protein_data("proteomics_data.csv")

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

# Compare with transcriptomic enrichment
transcriptomic_enriched = pathway_enrichment(diff_genes, pathway_network)
```

## Testing

Pathway analysis functionality is tested in `tests/test_networks_pathway.py`:

```bash
# Run pathway analysis tests
uv run pytest tests/test_networks_pathway.py -v

# Test specific functions
uv run pytest tests/test_networks_pathway.py::test_pathway_enrichment -v
uv run pytest tests/test_networks_pathway.py::test_load_pathway_database -v
```

## Related Documentation

- [Networks Overview](./index.md): Main networks module documentation
- [Community Detection](./community.md): Community detection in pathway networks
- [Graph Analysis](./graph.md): Basic graph operations and metrics
- [PPI Networks](./ppi.md): Protein-protein interaction analysis
