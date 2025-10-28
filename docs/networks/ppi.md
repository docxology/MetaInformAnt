# Networks: Protein-Protein Interaction Analysis

The protein-protein interaction (PPI) module provides tools for analyzing protein interaction networks, including data loading, network construction, and functional analysis.

## Network Construction

### ProteinNetwork Class

Container for protein interaction networks with metadata:

```python
from metainformant.networks import ProteinNetwork

# Create PPI network from interactions
interactions = [
    {"protein1": "P12345", "protein2": "P67890", "confidence": 0.9, "method": "yeast_two_hybrid"},
    {"protein1": "P12345", "protein2": "P11111", "confidence": 0.7, "method": "co_ip"},
    {"protein1": "P67890", "protein2": "P22222", "confidence": 0.8, "method": "mass_spec"}
]

ppi_network = ProteinNetwork()
ppi_network.add_interactions(interactions)

print(f"Network has {ppi_network.n_proteins} proteins and {ppi_network.n_interactions} interactions")
```

### load_string_interactions()

Load interactions from STRING database format:

```python
from metainformant.networks import load_string_interactions

# Load STRING database interactions
interactions = load_string_interactions(
    "string_interactions.tsv",
    min_confidence=0.7,     # Minimum confidence score
    evidence_filter=["experimental", "database"],  # Evidence types to include
    organism="9606"         # NCBI taxonomy ID (Homo sapiens)
)

# Create network from loaded interactions
ppi_network = ProteinNetwork()
ppi_network.add_interactions(interactions)
```

### predict_interactions()

Predict potential protein interactions using various methods:

```python
from metainformant.networks import predict_interactions

# Predict interactions using different methods
predicted = predict_interactions(
    known_network,
    prediction_method="domain_fusion",  # Method to use
    confidence_threshold=0.6,           # Minimum prediction confidence
    max_predictions=1000               # Limit number of predictions
)

print(f"Predicted {len(predicted)} potential interactions")

# Add predictions to network
ppi_network.add_interactions(predicted, interaction_type="predicted")
```

**Available prediction methods:**
- `'domain_fusion'`: Domain architecture similarity
- `'gene_fusion'`: Rosetta stone method
- `'neighborhood'`: Gene neighborhood conservation
- `'co_expression'`: Gene co-expression patterns
- `'co_occurrence'`: Phylogenetic co-occurrence

## Network Analysis

### get_protein_partners()

Find interaction partners for a protein:

```python
from metainformant.networks import get_protein_partners

# Get direct interaction partners
partners = get_protein_partners(ppi_network, "P12345")
print(f"Protein P12345 has {len(partners)} direct partners")

# Get partners with confidence scores
partners_with_confidence = get_protein_partners(
    ppi_network,
    "P12345",
    include_confidence=True
)

for partner, confidence in partners_with_confidence.items():
    print(f"  {partner}: {confidence:.2f}")
```

### filter_by_confidence()

Filter interactions by confidence or evidence:

```python
from metainformant.networks import filter_by_confidence

# Filter high-confidence interactions
high_conf_network = filter_by_confidence(
    ppi_network,
    min_confidence=0.8,
    evidence_types=["experimental"]  # Only experimental evidence
)

# Filter by interaction method
yeast_two_hybrid_only = filter_by_confidence(
    ppi_network,
    method_filter=["yeast_two_hybrid"]
)
```

## Functional Analysis

### Functional Enrichment

Analyze functional enrichment of protein modules:

```python
from metainformant.ontology import go

# Get protein modules (from community detection)
modules = community.detect_communities(ppi_network, method='leiden')

# Functional enrichment for each module
for module_id, proteins in modules.items():
    print(f"\nModule {module_id} ({len(proteins)} proteins):")

    # GO enrichment analysis
    enriched_terms = go.enrich_proteins(proteins, organism="human")

    for term in enriched_terms[:5]:  # Top 5 terms
        print(f"  {term['name']}: p={term['p_value']:.2e}")

    # Pathway enrichment
    pathway_enrichment = pathway.enrich_proteins(proteins)
    print(f"  Enriched pathways: {len(pathway_enrichment)}")
```

### Network Statistics

Get comprehensive network statistics:

```python
# Basic network statistics
stats = ppi_network.get_statistics()
print(f"Proteins: {stats['n_proteins']}")
print(f"Interactions: {stats['n_interactions']}")
print(f"Average degree: {stats['avg_degree']:.2f}")
print(f"Connected components: {stats['n_components']}")

# Interaction type distribution
print("\nInteraction types:")
for interaction_type, count in stats['interaction_types'].items():
    print(f"  {interaction_type}: {count}")

# Evidence type distribution
print("\nEvidence types:")
for evidence_type, count in stats['evidence_types'].items():
    print(f"  {evidence_type}: {count}")
```

## Integration Examples

### With Expression Data

```python
from metainformant.rna import workflow
from metainformant.networks import ppi

# Load expression data
expression_data = workflow.extract_expression_patterns(rna_data)

# Find differentially expressed proteins
diff_expr_proteins = expression_data[
    expression_data['p_value'] < 0.05
].index.tolist()

# Get their interaction partners
partners = set()
for protein in diff_expr_proteins:
    if protein in ppi_network.nodes:
        partners.update(get_protein_partners(ppi_network, protein))

print(f"Found {len(partners)} interaction partners of differentially expressed proteins")
```

### With Functional Annotation

```python
from metainformant.protein import proteomes
from metainformant.networks import ppi
from metainformant.ontology import go

# Load protein annotations
annotations = proteomes.load_annotations("protein_annotations.tsv")

# Analyze functional coherence of network modules
modules = community.detect_communities(ppi_network)

for module_id, proteins in modules.items():
    # Get annotations for module proteins
    module_annotations = annotations[annotations.index.isin(proteins)]

    # Analyze functional enrichment
    enriched_functions = go.enrich_annotations(module_annotations)

    print(f"Module {module_id}: {len(enriched_functions)} enriched functions")
```

## Data Sources

### STRING Database

```python
# Download STRING data (requires internet)
string_data = load_string_interactions(
    species="9606",  # Human
    version="11.5",
    min_confidence=0.7
)

# Create comprehensive PPI network
comprehensive_ppi = ProteinNetwork()
comprehensive_ppi.add_interactions(string_data)
```

### Custom Interaction Data

```python
# Load custom interaction data
custom_interactions = pd.read_csv("custom_ppi_data.csv")

# Convert to required format
interactions_list = []
for _, row in custom_interactions.iterrows():
    interactions_list.append({
        "protein1": row["protein_a"],
        "protein2": row["protein_b"],
        "confidence": row["confidence"],
        "method": row["method"],
        "evidence": row["evidence"]
    })

# Add to network
ppi_network.add_interactions(interactions_list)
```

## Visualization

### Network Visualization

```python
import matplotlib.pyplot as plt

# Visualize PPI network
fig, ax = plt.subplots(figsize=(12, 10))

# Color by functional annotation
node_colors = []
for protein in ppi_network.nodes:
    if protein in annotations.index:
        # Color by functional category
        category = annotations.loc[protein, 'functional_category']
        node_colors.append(category_colors.get(category, 'gray'))
    else:
        node_colors.append('lightblue')

# Draw network
nx.draw(
    ppi_network,
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

PPI functionality is tested in `tests/test_networks_ppi.py`:

```bash
# Run PPI tests
uv run pytest tests/test_networks_ppi.py -v

# Test specific functions
uv run pytest tests/test_networks_ppi.py::test_load_string_interactions -v
uv run pytest tests/test_networks_ppi.py::test_predict_interactions -v
```

## Related Documentation

- [Networks Overview](./index.md): Main networks module documentation
- [Community Detection](./community.md): Community detection in PPI networks
- [Graph Analysis](./graph.md): Basic graph operations and metrics
- [Regulatory Networks](./regulatory.md): Gene regulatory network analysis
