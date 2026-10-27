# Networks: Protein-Protein Interaction Analysis

The protein-protein interaction (PPI) module provides tools for analyzing protein interaction networks, including data loading, network construction, and functional analysis.

The PPI classes and functions below are importable from `metainformant.networks.interaction.ppi`, the public facade over the implementation module.

## Network Construction

### ProteinNetwork Class

Container for protein interaction networks with metadata. Interactions are added one at a time with `add_interaction(protein1, protein2, confidence=1.0, evidence_types=None, **attributes)`:

```python
from metainformant.networks.interaction.ppi import ProteinNetwork

# Create PPI network from interactions
interactions = [
    ("P12345", "P67890", 0.9, ["experimental"]),
    ("P12345", "P11111", 0.7, ["database"]),
    ("P67890", "P22222", 0.8, ["experimental", "textmining"]),
]

ppi_network = ProteinNetwork(name="demo_ppi")
for protein1, protein2, confidence, evidence in interactions:
    ppi_network.add_interaction(protein1, protein2, confidence, evidence)

# Per-protein metadata is stored on the graph nodes
ppi_network.add_protein_metadata("P12345", gene_name="GENE1", function="kinase")

print(f"Network has {len(ppi_network)} proteins and {len(ppi_network.interactions)} interactions")
```

`len(network)` is the number of proteins, `network.proteins` is the set of protein IDs, `network.interactions` holds `(protein1, protein2, attributes)` records, and `network.graph` is the underlying `networkx.Graph`. Two class methods build a network in one call: `ProteinNetwork.from_interactions([(protein1, protein2), ...])` (no confidence or evidence) and `ProteinNetwork.from_file(path, format="tsv")`, where `format` is one of `'tsv'`, `'csv'`, `'bioplex'` or `'intact'`.

### load_string_interactions()

Build a network from a STRING-style table. The function takes a pandas DataFrame with the columns `protein1`, `protein2` and `combined_score` (0-1000); the optional `experimental`, `database` and `textmining` columns become the evidence labels of each interaction. It does not read files or download anything.

```python
import pandas as pd
from metainformant.networks.interaction.ppi import load_string_interactions

# A STRING links table (for example pd.read_csv("string_links_detailed.txt", sep=" "))
interactions_df = pd.DataFrame({
    "protein1": ["P1", "P1", "P2", "P3"],
    "protein2": ["P2", "P3", "P3", "P4"],
    "experimental": [900, 0, 500, 100],
    "database": [0, 800, 0, 0],
    "textmining": [0, 0, 600, 300],
    "combined_score": [950, 800, 700, 350],
})

# Keep interactions with combined_score >= 700 (the default threshold is 400)
string_network = load_string_interactions(interactions_df, confidence_threshold=700)

for protein1, protein2, attrs in string_network.interactions:
    print(protein1, protein2, attrs["confidence"], attrs["evidence_types"])  # confidence = score / 1000
```

An optional `proteins_df` argument (columns `protein_id`, `gene_name`, `protein_name`) attaches protein metadata to the network.

### predict_interactions()

Predict potential interactions for a list of target proteins. The result maps each target protein to a list of `{"partner", "confidence", "method"}` dictionaries:

```python
from metainformant.networks.interaction.ppi import predict_interactions

# known_network is the ProteinNetwork built above
predicted = predict_interactions(
    ["P12345"],                       # Proteins to predict partners for
    known_network=ppi_network,        # Reference network (needed by topology-based methods)
    method="similarity",              # Prediction method
    threshold=0.1,                    # Minimum prediction confidence
    max_predictions_per_protein=10,   # Limit predictions per target
)

for target, candidates in predicted.items():
    print(f"{target}: {len(candidates)} predicted partners")

# Add predictions to the network, flagged as predicted
for target, candidates in predicted.items():
    for candidate in candidates:
        ppi_network.add_interaction(
            target, candidate["partner"], candidate["confidence"], ["predicted"], interaction_type="predicted"
        )
```

**Available prediction methods** (any other value raises `ValueError`):
- `'similarity'`: Jaccard similarity of network neighbours (needs `known_network`; a target that is not in the network gets an empty list). `'domain'` is accepted as an alias
- `'correlation'`: absolute Pearson correlation between per-protein feature vectors (needs `features={protein: [values, ...]}`). `'coexpression'` is accepted as an alias
- `'guilt-by-association'`: second-degree neighbours scored by shared neighbours (needs `known_network`)
- `'ml'`: placeholder that logs a warning and falls back to the `'correlation'` method (needs `features` and `known_network`)

## Network Analysis

### get_protein_partners()

Find interaction partners for a protein (a `ProteinNetwork` method):

```python
# Get direct interaction partners
partners = ppi_network.get_protein_partners("P12345")
print(f"Protein P12345 has {len(partners)} direct partners")

# Get only partners with confidence >= 0.8
confident_partners = ppi_network.get_protein_partners("P12345", min_confidence=0.8)

for partner in confident_partners:
    confidence = ppi_network.graph["P12345"][partner]["confidence"]
    print(f"  {partner}: {confidence:.2f}")
```

### filter_by_confidence()

Filter interactions by confidence or evidence (`ProteinNetwork` methods that return a new network):

```python
# Filter high-confidence interactions
high_conf_network = ppi_network.filter_by_confidence(0.8)

# Keep only interactions supported by experimental evidence
experimental_only = ppi_network.filter_by_evidence("experimental")
```

## Functional Analysis

### Functional Enrichment

Analyze functional enrichment of protein modules. `detect_communities` returns a `{protein: module_id}` mapping, so group it into modules first, then test each module against the `function` annotation stored in the protein metadata (`functional_enrichment_ppi`) or against GO term annotations (`go.enrich_genes`):

```python
from metainformant.networks.analysis.community import detect_communities
from metainformant.networks.interaction.ppi import ProteinNetwork, functional_enrichment_ppi
from metainformant.ontology.core import go

# A small annotated network: two kinase-rich and one transporter-rich cluster
ppi_network = ProteinNetwork(name="annotated_ppi")
for protein1, protein2, confidence in [
    ("A", "B", 0.9), ("B", "C", 0.8), ("A", "C", 0.85),
    ("D", "E", 0.9), ("E", "F", 0.7), ("D", "F", 0.8),
    ("C", "D", 0.4), ("G", "H", 0.9),
]:
    ppi_network.add_interaction(protein1, protein2, confidence, ["experimental"])
for protein, function in zip("ABCDEFGH", ["kinase"] * 3 + ["transport"] * 3 + ["kinase", "transport"]):
    ppi_network.add_protein_metadata(protein, function=function)

# Protein modules (community detection returns protein -> module id)
membership = detect_communities(ppi_network, method="louvain")
modules = {}
for protein, module_id in membership.items():
    modules.setdefault(module_id, []).append(protein)

# GO term annotations: GO term -> set of proteins
go_annotations = {"GO:0016301": {"A", "B", "C", "G"}, "GO:0005215": {"D", "E", "F", "H"}}

for module_id, proteins in modules.items():
    print(f"\nModule {module_id} ({len(proteins)} proteins):")

    # Enrichment of the 'function' metadata annotation
    enriched_functions = functional_enrichment_ppi(proteins, ppi_network, function_key="function", min_overlap=2)
    for function, stats in enriched_functions.items():
        print(f"  {function}: p={stats['p_value']:.2e}, enrichment={stats['enrichment_ratio']:.2f}")

    # GO enrichment (a DataFrame sorted by p-value)
    go_enrichment = go.enrich_genes(proteins, ppi_network.get_proteins(), go_annotations)
    for _, term in go_enrichment.head(5).iterrows():  # Top 5 terms
        print(f"  {term['go_term']}: p={term['p_value']:.2e}")
```

`detect_communities` uses Louvain when the optional `python-louvain` package is installed and falls back to greedy modularity otherwise (see [Community Detection](./community.md)). `functional_enrichment_ppi` returns `{function: {"count", "test_count", "background_count", "enrichment_ratio", "p_value", "proteins"}}`; it keeps functions seen at least `min_overlap` times in the module. `go.enrich_genes(genes, background, annotations, alpha=0.05, method="fisher")` returns a DataFrame with the columns `go_term`, `p_value`, `enrichment`, `observed`, `expected`, `term_size`, `query_size` and `bonferroni_p`.

### Network Statistics

Get network statistics:

```python
# Basic network statistics
summary = ppi_network.network_summary()
print(f"Proteins: {summary['n_proteins']}")
print(f"Interactions: {summary['n_interactions']}")
print(f"Average degree: {summary['average_degree']:.2f}")
print(f"Connected components: {summary['n_components']}")

# Confidence summary
stats = ppi_network.get_network_statistics()
print(f"Mean confidence: {stats['avg_confidence']:.2f}")

# Evidence type distribution
from collections import Counter

evidence_counts = Counter(
    evidence for _, _, attrs in ppi_network.interactions for evidence in attrs["evidence_types"]
)
print("\nEvidence types:")
for evidence_type, count in evidence_counts.items():
    print(f"  {evidence_type}: {count}")
```

For degree, connectivity, clustering and centrality summaries of the underlying graph, use `ppi_network_analysis(ppi_network.graph)`; `find_ppi_hubs(ppi_network.graph, percentile=95.0)` returns `(protein, degree)` tuples for hub proteins.

## Integration Examples

### With Expression Data

```python-snippet
import pandas as pd
from metainformant.rna.analysis import expression

# counts: raw count matrix with genes as rows and samples as columns
counts = pd.read_csv("counts.csv", index_col=0)
conditions = ["control", "control", "control", "treated", "treated", "treated"]

# Find differentially expressed proteins (genes) between the two conditions
de_results = expression.differential_expression(counts, conditions)
diff_expr_proteins = de_results.loc[de_results["adjusted_p_value"] < 0.05, "gene"].tolist()

# Get their interaction partners
partners = set()
for protein in diff_expr_proteins:
    if protein in ppi_network:
        partners.update(ppi_network.get_protein_partners(protein))

print(f"Found {len(partners)} interaction partners of differentially expressed proteins")
```

### With Functional Annotation

```python-snippet
import pandas as pd
from metainformant.networks.analysis.community import detect_communities
from metainformant.ontology.core import go

# Protein annotations table: one row per (protein, GO term) pair
annotation_table = pd.read_csv("protein_annotations.tsv", sep="\t")  # columns: protein, go_term

# GO term -> set of annotated proteins
go_annotations = {
    term: set(group["protein"]) for term, group in annotation_table.groupby("go_term")
}

# Analyze functional coherence of network modules
membership = detect_communities(ppi_network)
modules = {}
for protein, module_id in membership.items():
    modules.setdefault(module_id, []).append(protein)

for module_id, proteins in modules.items():
    enrichment = go.enrich_genes(proteins, ppi_network.get_proteins(), go_annotations)
    n_enriched = int((enrichment["bonferroni_p"] <= 0.05).sum()) if len(enrichment) else 0
    print(f"Module {module_id}: {n_enriched} enriched GO terms")
```

## Data Sources

### STRING Database

```python-snippet
import pandas as pd
from metainformant.networks.interaction.ppi import load_string_interactions

# STRING exports are downloaded separately from string-db.org; the loader takes the parsed table
string_data = pd.read_csv("string_links_detailed.txt", sep=" ")  # protein1 protein2 ... combined_score

# Create comprehensive PPI network (combined_score >= 700)
comprehensive_ppi = load_string_interactions(string_data, confidence_threshold=700)
```

### Custom Interaction Data

```python
import pandas as pd
from metainformant.networks.interaction.ppi import ProteinNetwork

# Load custom interaction data (for example pd.read_csv("custom_ppi_data.csv"))
custom_interactions = pd.DataFrame({
    "protein_a": ["P1", "P2"],
    "protein_b": ["P2", "P3"],
    "confidence": [0.9, 0.6],
    "evidence": ["experimental", "database"],
    "method": ["yeast_two_hybrid", "curated"],
})

# Add to network; extra columns can ride along as edge attributes
ppi_network = ProteinNetwork(name="custom_ppi")
for _, row in custom_interactions.iterrows():
    ppi_network.add_interaction(
        row["protein_a"],
        row["protein_b"],
        row["confidence"],
        [row["evidence"]],
        method=row["method"],
    )
```

## Visualization

### Network Visualization

```python
import matplotlib.pyplot as plt
import networkx as nx
from metainformant.networks.interaction.ppi import ProteinNetwork

# A small annotated network
ppi_network = ProteinNetwork(name="viz_ppi")
for protein1, protein2, confidence in [("A", "B", 0.9), ("B", "C", 0.8), ("C", "D", 0.7), ("D", "A", 0.6)]:
    ppi_network.add_interaction(protein1, protein2, confidence, ["experimental"])
for protein, category in {"A": "kinase", "B": "kinase", "C": "transport"}.items():
    ppi_network.add_protein_metadata(protein, functional_category=category)
category_colors = {"kinase": "tomato", "transport": "steelblue"}

# Visualize PPI network
fig, ax = plt.subplots(figsize=(12, 10))

# Color by functional annotation (the networkx graph is ppi_network.graph)
node_colors = []
for protein in ppi_network.graph.nodes:
    category = ppi_network.get_protein_metadata(protein).get("functional_category")
    node_colors.append(category_colors.get(category, "lightblue"))

# Draw network
nx.draw(
    ppi_network.graph,
    ax=ax,
    node_color=node_colors,
    node_size=100,
    edge_color='gray',
    alpha=0.7,
    with_labels=False
)

plt.title('Protein-Protein Interaction Network')
plt.show()
```

## Testing

PPI functionality is tested in `tests/networks/test_networks_ppi.py`:

```bash
# Run PPI tests
uv run pytest tests/networks/test_networks_ppi.py -v

# Test specific classes and functions
uv run pytest tests/networks/test_networks_ppi.py::TestLoadStringInteractions -v
uv run pytest tests/networks/test_networks_ppi.py::TestPredictInteractions::test_predict_interactions_basic -v
```

## Related Documentation

- [Networks Overview](./index.md): Main networks module documentation
- [Community Detection](./community.md): Community detection in PPI networks
- [Graph Analysis](./graph.md): Basic graph operations and metrics
- [Regulatory Networks](./regulatory.md): Gene regulatory network analysis
