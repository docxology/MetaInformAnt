# Information Theory Examples

This document provides real-world examples and workflow tutorials for using the information module in biological data analysis.

## Example 1: DNA Sequence Complexity Analysis

### Problem
Analyze the information content and complexity of DNA sequences from different genomic regions.

### Solution

```python
from metainformant.dna.sequence.core import read_fasta
from metainformant.information.metrics.analysis import (
    analyze_sequence_information,
    information_profile,
)

# Load sequences
dna_seqs = read_fasta("data/genomic_regions.fasta")

# Analyze each sequence
for seq_id, seq in dna_seqs.items():
    analysis = analyze_sequence_information(seq, k_values=[1, 2, 3])

    # Extract entropy for different k-mer sizes
    entropy_1mer = analysis["k_mer_analysis"]["k1"]["entropy"]
    entropy_2mer = analysis["k_mer_analysis"]["k2"]["entropy"]

    print(f"{seq_id}: 1-mer entropy={entropy_1mer:.3f}, 2-mer entropy={entropy_2mer:.3f}")

# Overall profile (sequences must be aligned to the same length)
profile = information_profile(list(dna_seqs.values()), k=2)
print(f"Overall mean entropy: {profile['statistics']['mean_entropy']:.3f} bits")
```

### Performance Tips
- For large sequences, use k=1 or k=2 to reduce computation
- Process sequences in batches for memory efficiency

## Example 2: Gene Expression Entropy Analysis

### Problem
Identify genes with highly variable expression patterns using entropy.

### Solution

```python
from metainformant.information.integration import rna_integration
import numpy as np

# Load expression data (samples x genes)
expression = np.load("data/expression_matrix.npy")

# Calculate discretized entropy for each gene
results = rna_integration(expression, method="entropy")

gene_entropies = results["gene_entropies"]  # list of floats, one per gene
metrics = results["integrated_metrics"]
mean_entropy = metrics["gene_entropy_mean"]
std_entropy = metrics["gene_entropy_std"]

# Find highly variable genes (high entropy)
high_variability_genes = [
    (i, h) for i, h in enumerate(gene_entropies) if h > mean_entropy + std_entropy
]

print(f"Found {len(high_variability_genes)} highly variable genes")
for gene_index, entropy in high_variability_genes[:10]:
    print(f"gene_{gene_index}: entropy={entropy:.3f}")
```

## Example 3: Feature Selection using Mutual Information

### Problem
Select informative features for a classification task using mutual information.

### Solution

```python
from metainformant.information.integration import ml_integration
import numpy as np
from sklearn.model_selection import train_test_split

# Load data
X = np.load("data/features.npy")
y = np.load("data/labels.npy")

# Split data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# Calculate MI for each feature
results = ml_integration(X_train, y_train, method="feature_mi")

# Select top features
top_k = 50
top_features = results["top_features"][:top_k]  # list of {"index": int, "mi": float}
top_indices = [f["index"] for f in top_features]

# Use selected features
X_train_selected = X_train[:, top_indices]
X_test_selected = X_test[:, top_indices]

print(f"Selected {len(top_indices)} features based on mutual information")
```

## Example 4: Multi-Omics Information Integration

### Problem
Compare information content across different omics platforms and analyze cross-platform relationships.

### Solution

```python
from metainformant.information.integration import multiomics_integration
import numpy as np

# Load multi-omics data
genomics = np.load("data/genomics.npy")  # (samples, variants)
transcriptomics = np.load("data/transcriptomics.npy")  # (samples, genes)
proteomics = np.load("data/proteomics.npy")  # (samples, proteins)

# Calculate platform entropy
results = multiomics_integration(
    genomics_data=genomics,
    transcriptomics_data=transcriptomics,
    proteomics_data=proteomics,
    method="platform_entropy"
)

# Compare information content
print("Platform Information Content:")
print(f"Genomics: {results['genomics_entropy']:.3f} bits")
print(f"Transcriptomics: {results['transcriptomics_entropy']:.3f} bits")
print(f"Proteomics: {results['proteomics_entropy']:.3f} bits")

# Cross-platform mutual information with feature selection
# Use first feature (default)
mi_results = multiomics_integration(
    genomics_data=genomics,
    transcriptomics_data=transcriptomics,
    method="cross_platform_mi"
)

# Results always in matrix format
mi_matrix = mi_results["genomics_transcriptomics_mi_matrix"]
print(f"Genomics-Transcriptomics MI Matrix shape: {len(mi_matrix)} x {len(mi_matrix[0])}")
print(f"Mean MI: {mi_results['genomics_transcriptomics_mean_mi']:.3f}")
print(f"Max MI: {mi_results['genomics_transcriptomics_max_mi']:.3f}")
print(f"Min MI: {mi_results['genomics_transcriptomics_min_mi']:.3f}")

# Use specific features for comprehensive analysis
mi_results_multi = multiomics_integration(
    genomics_data=genomics,
    transcriptomics_data=transcriptomics,
    method="cross_platform_mi",
    feature_indices={
        "genomics": [0, 1, 2, 3],  # Compare first 4 variants
        "transcriptomics": [0, 5, 10, 15]  # Compare specific genes
    }
)

# Results include MI matrix and summary statistics
mi_matrix_multi = mi_results_multi["genomics_transcriptomics_mi_matrix"]
print(f"MI Matrix shape: {len(mi_matrix_multi)} x {len(mi_matrix_multi[0])}")
print(f"Mean MI: {mi_results_multi['genomics_transcriptomics_mean_mi']:.3f}")
```

## Example 5: Network Information Analysis

### Problem
Analyze information flow and structure in a network built from time-series observations.

### Solution

```python
from metainformant.information.network_info.information_flow import (
    information_flow_network,
    mutual_information_network,
)
from metainformant.information.integration.networks import network_entropy
import networkx as nx

# Build a network
G = nx.Graph()
G.add_edges_from([("P1", "P2"), ("P2", "P3"), ("P3", "P4"), ("P1", "P3")])

# Calculate network entropy
entropy = network_entropy(G)
print(f"Network entropy: {entropy:.3f} bits")

# Build a directed information flow network from multivariate time series
time_series = {
    "P1": [...],  # equal-length observation series per node
    "P2": [...],
    "P3": [...],
    "P4": [...],
}
flow = information_flow_network(time_series, method="transfer_entropy", threshold=0.05)
print(f"Significant directed edges: {flow['edges']}")
print(f"Hub nodes: {flow['hub_nodes']}")

# Build an undirected MI network from a data matrix (observations x variables)
data_matrix = np.array([...])  # shape (n_observations, 4)
mi_net = mutual_information_network(
    data_matrix, variable_names=["P1", "P2", "P3", "P4"], threshold=0.05
)
print(f"Significant pairs: {mi_net['significant_pairs']}")
```

## Example 6: Single-Cell Data Analysis

### Problem
Analyze information content of cell type distributions and gene expression patterns.

### Solution

```python
from metainformant.information.integration import singlecell_integration
import numpy as np

# Load single-cell data (cells x genes count matrix)
count_matrix = np.load("data/counts.npy")
cell_types = ["T cell", "B cell", "T cell", "monocyte"]  # one label per cell

# Calculate cell type entropy
results = singlecell_integration(
    count_matrix,
    cell_types=cell_types,
    method="cell_type_entropy"
)

print(f"Cell type entropy: {results['cell_type_entropy']:.3f} bits")
print(f"Number of cell types: {results['num_cell_types']}")
print(f"Per-type mean gene entropy: {results['per_type_entropy']}")

# Calculate gene expression entropy
gene_results = singlecell_integration(count_matrix, method="gene_entropy")
gene_entropies = gene_results["gene_entropies"]  # list of floats, one per gene

print(f"Analyzed {len(gene_entropies)} genes")
```

## Example 7: Batch Processing Multiple Datasets

### Problem
Process multiple FASTA files and generate a comprehensive report.

### Solution

```python
from pathlib import Path
from metainformant.dna.sequence.core import read_fasta
from metainformant.information.workflow import information_workflow, information_report

fasta_files = list(Path("data/sequences").glob("*.fasta"))
all_results = {}

for fasta_file in fasta_files:
    # Read sequences
    seqs = read_fasta(str(fasta_file))
    seq_list = list(seqs.values())

    # Analyze
    results = information_workflow(
        seq_list,
        k_values=[1, 2],
        output_dir=f"output/information/{fasta_file.stem}"
    )

    all_results[fasta_file.stem] = results

# Generate a report for the last processed dataset
information_report(
    results,
    output_path="docs/information/comprehensive_report.md",
    format="markdown"
)
```

## Example 8: Rényi Entropy Spectrum Analysis

### Problem
Analyze how entropy varies with different orders using Rényi entropy.

### Solution

```python
from metainformant.information.metrics.core.syntactic import renyi_entropy, shannon_entropy
from collections import Counter

# Get k-mer distribution from sequence
sequence = "ATCGATCGATCG" * 100
kmer_counts = Counter()
k = 2
for i in range(len(sequence) - k + 1):
    kmer_counts[sequence[i:i+k]] += 1

# Convert to probabilities
total = sum(kmer_counts.values())
probs = [count / total for count in kmer_counts.values()]

# Compute the Rényi spectrum (alpha must differ from 1)
for alpha in [0.1, 0.5, 1.5, 2.0, 3.0, 5.0, 10.0]:
    h_alpha = renyi_entropy(probs, alpha=alpha)
    print(f"Rényi (α={alpha}): {h_alpha:.3f}")

# Compare with Shannon entropy (the alpha -> 1 limit)
h_shannon = shannon_entropy(probs)
print(f"Shannon (α=1 limit): {h_shannon:.3f}")
```

Note: `renyi_entropy` raises `ValueError` for `alpha == 1`; use `shannon_entropy`
for that limit. A plotting helper is not provided by this module.

## Example 9: Continuous Data Analysis

### Problem
Analyze continuous gene expression data using differential entropy.

### Solution

```python
from metainformant.information.metrics.core.continuous import (
    differential_entropy,
    mutual_information_continuous,
)
import numpy as np

# Load continuous expression data
gene1_expr = np.load("data/gene1_expression.npy")
gene2_expr = np.load("data/gene2_expression.npy")

# Calculate differential entropy
h1 = differential_entropy(gene1_expr, method="histogram", bins=20)
h2 = differential_entropy(gene2_expr, method="histogram", bins=20)

print(f"Gene 1 entropy: {h1:.3f} nats")
print(f"Gene 2 entropy: {h2:.3f} nats")

# Calculate continuous mutual information
mi = mutual_information_continuous(gene1_expr, gene2_expr, method="histogram", bins=20)
print(f"Mutual information: {mi:.3f} nats")
```

## Example 10: Semantic Similarity for GO Terms

### Problem
Calculate semantic similarity between Gene Ontology terms using information content.

### Solution

```python
from metainformant.information.metrics.advanced.semantic import (
    information_content,
    information_content_from_annotations,
    semantic_similarity,
    semantic_similarity_matrix,
)

# Annotate genes with GO terms
gene_annotations = {
    "gene1": {"GO:0008150", "GO:0003674"},
    "gene2": {"GO:0008150", "GO:0005524"},
    "gene3": {"GO:0003674", "GO:0005524"},
}

# IC of a single term from annotations (-log2 of annotation frequency)
ic_bp = information_content_from_annotations(gene_annotations, "GO:0008150")
print(f"IC(GO:0008150): {ic_bp:.3f} bits")

# Build per-term IC values and a term -> parent-terms hierarchy
term_frequencies = {"GO:0008150": 3, "GO:0003674": 2, "GO:0005524": 2}
term_ic = {t: information_content(term_frequencies, t) for t in term_frequencies}
hierarchy = {"GO:0003674": {"GO:0008150"}, "GO:0005524": {"GO:0008150"}}

# Similarity between specific terms
go_term1 = "GO:0003674"
go_term2 = "GO:0005524"
similarity = semantic_similarity(go_term1, go_term2, term_ic, hierarchy)
print(f"Similarity between {go_term1} and {go_term2}: {similarity:.3f}")

# Pairwise similarity matrix (nested lists, not a numpy array)
all_terms = list(term_ic)
similarity_matrix = semantic_similarity_matrix(all_terms, term_ic, hierarchy)
print(f"Similarity matrix: {len(similarity_matrix)} x {len(similarity_matrix[0])}")
```

## Performance Tips

1. **Batch Processing**: Use `batch_entropy_analysis()` for processing many sequences
2. **Memory Management**: For large datasets, process in chunks
3. **Binning**: Use appropriate bin counts for continuous data (Sturges' rule: 1 + log2(n))
4. **Bias Correction**: Always use bias correction for small samples (n < 100)
5. **Parallel Processing**: Consider using parallel processing for batch operations

## Integration Patterns

### Pattern 1: DNA → Information

```python
from metainformant.dna.sequence.core import read_fasta
from metainformant.information.metrics.analysis import information_profile

seqs = read_fasta("data/sequences.fasta")
profile = information_profile(list(seqs.values()), k=2)
print(f"Mean entropy: {profile['statistics']['mean_entropy']:.3f} bits")
```

### Pattern 2: Expression → Information → ML

```python
from metainformant.information.integration import rna_integration, ml_integration

# Get entropy for feature selection
expression = np.load("data/expression.npy")
entropy_results = rna_integration(expression, method="entropy")

# Use entropy for feature selection
X = expression
y = labels
mi_results = ml_integration(X, y, method="feature_mi")
```

These examples demonstrate practical applications of information theory in biological data analysis.
