# Information Theory Examples

This document provides real-world examples and workflow tutorials for using the information module in biological data analysis.

## Example 1: DNA Sequence Complexity Analysis

### Problem
Analyze the information content and complexity of DNA sequences from different genomic regions.

### Solution

```python
from metainformant.dna import sequences
from metainformant.information import information_profile, analyze_sequence_information

# Load sequences
dna_seqs = sequences.read_fasta("data/genomic_regions.fasta")

# Analyze each sequence
for seq_id, seq in dna_seqs.items():
    analysis = analyze_sequence_information(seq, k_values=[1, 2, 3])
    
    # Extract entropy for different k-mer sizes
    entropy_1mer = analysis["kmer_analyses"][1]["entropy"]
    entropy_2mer = analysis["kmer_analyses"][2]["entropy"]
    
    print(f"{seq_id}: 1-mer entropy={entropy_1mer:.3f}, 2-mer entropy={entropy_2mer:.3f}")

# Overall profile
profile = information_profile(list(dna_seqs.values()), k=2)
print(f"Overall entropy: {profile['entropy']:.3f} bits")
print(f"Sequence complexity: {profile['sequence_complexity']:.3f}")
```

### Performance Tips
- For large sequences, use k=1 or k=2 to reduce computation
- Process sequences in batches for memory efficiency

## Example 2: Gene Expression Entropy Analysis

### Problem
Identify genes with highly variable expression patterns using entropy.

### Solution

```python
from metainformant.information import rna_integration
import numpy as np

# Load expression data (samples x genes)
expression = np.load("data/expression_matrix.npy")
gene_names = np.load("data/gene_names.npy")

# Calculate entropy for each gene
results = rna_integration(expression, gene_names=gene_names.tolist(), method="entropy")

# Find highly variable genes (high entropy)
mean_entropy = results["mean_entropy"]
high_variability_genes = [
    g for g in results["gene_entropies"] 
    if g["entropy"] > mean_entropy + results["std_entropy"]
]

print(f"Found {len(high_variability_genes)} highly variable genes")
for gene in high_variability_genes[:10]:
    print(f"{gene['gene_name']}: entropy={gene['entropy']:.3f}")
```

## Example 3: Feature Selection using Mutual Information

### Problem
Select informative features for a classification task using mutual information.

### Solution

```python
from metainformant.information import ml_integration
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
top_features = results["top_features"][:top_k]

# Use selected features
X_train_selected = X_train[:, top_features]
X_test_selected = X_test[:, top_features]

print(f"Selected {len(top_features)} features based on mutual information")
```

## Example 4: Multi-Omics Information Integration

### Problem
Compare information content across different omics platforms.

### Solution

```python
from metainformant.information import multiomics_integration
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

# Cross-platform mutual information
mi_results = multiomics_integration(
    genomics_data=genomics,
    transcriptomics_data=transcriptomics,
    method="cross_platform_mi"
)

if "genomics_transcriptomics_mi" in mi_results:
    print(f"Genomics-Transcriptomics MI: {mi_results['genomics_transcriptomics_mi']:.3f}")
```

## Example 5: Network Information Analysis

### Problem
Analyze information flow and structure in a protein-protein interaction network.

### Solution

```python
from metainformant.networks import load_string_interactions
from metainformant.information.networks import (
    network_entropy,
    information_flow,
    mutual_information_network,
)
from metainformant.information.visualization import plot_mi_network

# Load PPI network
ppi_network = load_string_interactions("data/string_interactions.tsv", score_threshold=400)
network = ppi_network.create_network()

# Calculate network entropy
entropy = network_entropy(network)
print(f"Network entropy: {entropy:.3f} bits")

# Analyze information flow
flow = information_flow(network)
print(f"Path length entropy: {flow['path_length_entropy']:.3f}")
print(f"Mean path length: {flow['mean_path_length']:.2f}")

# Calculate MI matrix
mi_matrix = mutual_information_network(network)

# Visualize
protein_ids = list(network.nodes())[:50]  # Limit for visualization
plot_mi_network(
    mi_matrix[:50, :50],
    labels=protein_ids,
    threshold=0.1,
    output_path="output/information/ppi_mi_network.png"
)
```

## Example 6: Single-Cell Data Analysis

### Problem
Analyze information content of cell type distributions and gene expression patterns.

### Solution

```python
from metainformant.singlecell import preprocessing
from metainformant.information import singlecell_integration

# Load single-cell data
count_matrix = preprocessing.load_count_matrix("data/counts.h5ad")
cell_types = preprocessing.get_cell_types(count_matrix)

# Calculate cell type entropy
results = singlecell_integration(
    count_matrix,
    cell_types=cell_types,
    method="cell_type_entropy"
)

print(f"Cell type entropy: {results['cell_type_entropy']:.3f} bits")
print(f"Number of cell types: {results['num_cell_types']}")

# Calculate gene expression entropy
gene_results = singlecell_integration(
    count_matrix,
    method="gene_entropy"
)

# Find highly variable genes
high_entropy_genes = [
    g for g in gene_results["gene_entropies"]
    if g["entropy"] > gene_results["mean_entropy"]
]

print(f"Found {len(high_entropy_genes)} highly variable genes")
```

## Example 7: Batch Processing Multiple Datasets

### Problem
Process multiple FASTA files and generate a comprehensive report.

### Solution

```python
from pathlib import Path
from metainformant.information import batch_entropy_analysis, information_workflow, information_report

# Process multiple files
fasta_files = list(Path("data/sequences").glob("*.fasta"))
all_results = {}

for fasta_file in fasta_files:
    # Read sequences
    from metainformant.dna import sequences
    seqs = sequences.read_fasta(str(fasta_file))
    seq_list = list(seqs.values())
    
    # Analyze
    results = information_workflow(
        seq_list,
        k_values=[1, 2],
        output_dir=f"output/information/{fasta_file.stem}"
    )
    
    all_results[fasta_file.stem] = results

# Generate comprehensive report
information_report(
    all_results,
    output_path="output/information/comprehensive_report.md",
    format="markdown"
)
```

## Example 8: Rényi Entropy Spectrum Analysis

### Problem
Analyze how entropy varies with different orders using Rényi entropy.

### Solution

```python
from metainformant.information import renyi_entropy
from metainformant.information.visualization import plot_renyi_spectrum
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

# Plot Rényi spectrum
plot_renyi_spectrum(
    probs,
    alpha_range=[0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0],
    output_path="output/information/renyi_spectrum.png"
)

# Compare with Shannon entropy
from metainformant.information import shannon_entropy
h_shannon = shannon_entropy(probs)
h_renyi_1 = renyi_entropy(probs, alpha=1.0)
print(f"Shannon: {h_shannon:.3f}, Rényi (α=1): {h_renyi_1:.3f}")
```

## Example 9: Continuous Data Analysis

### Problem
Analyze continuous gene expression data using differential entropy.

### Solution

```python
from metainformant.information import differential_entropy, mutual_information_continuous
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
from metainformant.information import (
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

# Calculate information content
term_ic = information_content_from_annotations(gene_annotations)

# Calculate similarity between specific terms
go_term1 = "GO:0008150"
go_term2 = "GO:0003674"
similarity = semantic_similarity(go_term1, go_term2, term_ic)
print(f"Similarity between {go_term1} and {go_term2}: {similarity:.3f}")

# Pairwise similarity matrix
all_terms = list(term_ic.keys())
similarity_matrix = semantic_similarity_matrix(all_terms, term_ic)
print(f"Similarity matrix shape: {similarity_matrix.shape}")
```

## Performance Tips

1. **Batch Processing**: Use `batch_entropy_analysis()` for processing many sequences
2. **Memory Management**: For large datasets, process in chunks
3. **Binning**: Use appropriate bin counts for continuous data (Sturges' rule: 1 + log2(n))
4. **Bias Correction**: Always use bias correction for small samples (n < 100)
5. **Parallel Processing**: Consider using parallel processing for batch operations

## Integration Patterns

### Pattern 1: DNA → Information → Visualization

```python
from metainformant.dna import sequences
from metainformant.information import information_profile
from metainformant.information.visualization import plot_information_profile

seqs = sequences.read_fasta("data/sequences.fasta")
profile = information_profile(list(seqs.values()), k=2)
plot_information_profile(profile, output_path="output/information/profile.png")
```

### Pattern 2: Expression → Information → ML

```python
from metainformant.information import rna_integration, ml_integration

# Get entropy for feature selection
expression = np.load("data/expression.npy")
entropy_results = rna_integration(expression, method="entropy")

# Use entropy for feature selection
X = expression
y = labels
mi_results = ml_integration(X, y, method="feature_mi")
```

These examples demonstrate practical applications of information theory in biological data analysis.

