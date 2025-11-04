# Information Theory Module

The `information` module provides comprehensive information-theoretic methods for analyzing biological data, including both syntactic (Shannon entropy, mutual information) and semantic (information content, semantic similarity) information measures.

## Overview

This module implements fundamental information theory concepts for biological data analysis:
- **Syntactic Information Theory**: Shannon entropy, mutual information, KL divergence, transfer entropy
- **Semantic Information Theory**: Information content, semantic similarity, semantic entropy
- **Network Integration**: Information flow, network entropy, information-based community detection
- **Visualization**: Integration with visualization module for plotting information measures

## Key Components

### Syntactic Information Theory (`syntactic.py`)

Fundamental information-theoretic measures for analyzing data distributions.

#### Shannon Entropy

```python
from metainformant.information import shannon_entropy

# Calculate entropy of a probability distribution
probs = [0.5, 0.3, 0.2]
entropy = shannon_entropy(probs)  # bits

# From counts
from metainformant.information.syntactic import shannon_entropy_from_counts
counts = {"A": 50, "T": 30, "G": 20}
entropy = shannon_entropy_from_counts(counts)
```

#### Mutual Information

```python
from metainformant.information import mutual_information

# Calculate mutual information between two sequences
x = [0, 1, 0, 1, 0, 1]
y = [1, 0, 1, 0, 1, 0]  # Perfect negative correlation
mi = mutual_information(x, y)
```

#### Kullback-Leibler Divergence

```python
from metainformant.information import kl_divergence

# Compare two probability distributions
p = [0.5, 0.3, 0.2]
q = [0.4, 0.4, 0.2]
divergence = kl_divergence(p, q)
```

#### Transfer Entropy

```python
from metainformant.information import transfer_entropy

# Measure information transfer between time series
x = [0, 1, 0, 1, 0, 1]
y = [0, 0, 1, 1, 0, 0]
te = transfer_entropy(x, y, lag=1)
```

### Semantic Information Theory (`semantic.py`)

Information measures for hierarchical and ontological data structures.

#### Information Content

```python
from metainformant.information import information_content

# Calculate information content of GO terms
term_frequencies = {"GO:0008150": 1000, "GO:0003674": 100}
ic = information_content(term_frequencies, "GO:0003674")
```

#### Semantic Similarity

```python
from metainformant.information import semantic_similarity, semantic_similarity_matrix

# Calculate similarity between terms
term_ic = {"A": 2.0, "B": 2.0, "C": 1.0}
hierarchy = {"A": {"C"}, "B": {"C"}}
similarity = semantic_similarity("A", "B", term_ic, hierarchy)

# Pairwise similarity matrix
terms = ["A", "B", "C"]
similarity_matrix = semantic_similarity_matrix(terms, term_ic, hierarchy)
```

### Analysis Functions (`analysis.py`)

High-level analysis workflows for biological sequences and data.

#### Information Profile

```python
from metainformant.information import information_profile

# Analyze sequence information content
sequences = ["ATCGATCG", "AAAAAA", "ATCGATCG"]
profile = information_profile(sequences, k=2)
print(f"Entropy: {profile['entropy']:.3f} bits")
print(f"Unique k-mers: {profile['unique_kmers']}")
print(f"Complexity: {profile['sequence_complexity']:.3f}")
```

#### Sequence Information Analysis

```python
from metainformant.information import analyze_sequence_information

# Comprehensive sequence analysis
sequence = "ATCGATCGATCG"
analysis = analyze_sequence_information(sequence, k_values=[1, 2, 3])
```

#### Sequence Comparison

```python
from metainformant.information import compare_sequences_information

# Compare two sequences
seq1 = "ATCGATCG"
seq2 = "ATCGATCG"
comparison = compare_sequences_information(seq1, seq2, k=2)
print(f"KL divergence: {comparison['kl_divergence']:.3f}")
print(f"Mutual information: {comparison['mutual_information']:.3f}")
```

### Network Integration (`networks.py`)

Information-theoretic analysis of biological networks.

#### Network Entropy

```python
from metainformant.information.networks import network_entropy
import networkx as nx

# Create network
G = nx.karate_club_graph()

# Calculate entropy of degree distribution
entropy = network_entropy(G)

# Or entropy of node attributes
entropy_attr = network_entropy(G, attribute="club")
```

#### Information Flow

```python
from metainformant.information.networks import information_flow

# Analyze information flow through network
flow = information_flow(G, source_nodes=["Mr. Hi"], target_nodes=list(G.nodes()))
print(f"Path length entropy: {flow['path_length_entropy']:.3f}")
```

#### Mutual Information Network

```python
from metainformant.information.networks import mutual_information_network

# Calculate MI matrix from network structure
mi_matrix = mutual_information_network(G)
```

### Visualization Integration (`visualization.py`)

Plotting functions for information-theoretic analysis.

#### Plot Entropy Distribution

```python
from metainformant.information.visualization import plot_entropy_distribution

# Plot distribution of entropy values
entropies = [2.5, 3.1, 2.8, 3.0, 2.9]
plot_entropy_distribution(entropies, output_path="output/information/entropy.png")
```

#### Plot Mutual Information Matrix

```python
from metainformant.information.visualization import plot_mutual_information_matrix
import numpy as np

# Plot MI matrix as heatmap
mi_matrix = np.random.rand(10, 10)
labels = [f"Gene_{i}" for i in range(10)]
plot_mutual_information_matrix(mi_matrix, labels, output_path="output/information/mi_matrix.png")
```

#### Plot Information Profile

```python
from metainformant.information import information_profile
from metainformant.information.visualization import plot_information_profile

# Analyze and visualize
profile = information_profile(sequences, k=2)
plot_information_profile(profile, output_path="output/information/profile.png")
```

## Integration with Other Modules

### With DNA Module

```python
from metainformant.dna import sequences
from metainformant.information import information_profile

# Analyze DNA sequence information
dna_seqs = sequences.read_fasta("data/sequences.fasta")
profile = information_profile(list(dna_seqs.values()), k=3)
```

### With Networks Module

```python
from metainformant.networks import create_network
from metainformant.information.networks import network_entropy, information_flow

# Create network and analyze information
network = create_network(["A", "B", "C"], directed=False)
network.add_edge("A", "B")
network.add_edge("B", "C")

entropy = network_entropy(network)
flow = information_flow(network)
```

### With Ontology Module

```python
from metainformant.ontology import go
from metainformant.information import information_content, semantic_similarity

# Analyze GO term information content
go_terms = go.load_go_terms("go.obo")
term_freqs = {term.id: term.frequency for term in go_terms}
ic = information_content(term_freqs, "GO:0008150")
```

### With Visualization Module

```python
from metainformant.information import analyze_sequence_information
from metainformant.information.visualization import plot_information_profile

# Analyze and visualize
analysis = analyze_sequence_information("ATCGATCGATCG")
profile = {"entropy": analysis["kmer_analyses"][1]["entropy"]}
plot_information_profile(profile)
```

## Mathematical Background

### Shannon Entropy

Shannon entropy measures the uncertainty in a probability distribution:
```
H(X) = -Σ p(x) × log₂(p(x))
```

### Mutual Information

Mutual information measures shared information between variables:
```
I(X; Y) = H(X) + H(Y) - H(X, Y) = H(X) - H(X|Y)
```

### Kullback-Leibler Divergence

KL divergence measures how different two distributions are:
```
D_KL(P||Q) = Σ p(x) × log(p(x) / q(x))
```

### Transfer Entropy

Transfer entropy measures information transfer between time series:
```
T(Y → X) = H(X_t | X_{t-1}) - H(X_t | X_{t-1}, Y_{t-1})
```

## Performance Features

- **Efficient Algorithms**: Optimized implementations for large datasets
- **Memory Conscious**: Handles large sequences and networks efficiently
- **Vectorized Operations**: NumPy-based implementations where applicable
- **Parallel Processing**: Supports parallel computation for multiple sequences

## Testing

Comprehensive tests ensure mathematical correctness:
- Verification against known analytical solutions
- Edge case handling (empty sequences, zero probabilities)
- Integration tests with other modules
- Performance benchmarks

## Dependencies

- **Core**: NumPy for numerical computations
- **Optional**: NetworkX for network analysis
- **Visualization**: Matplotlib for plotting (via visualization module)

### Continuous Information Theory (`continuous.py`)

Methods for analyzing continuous probability distributions.

#### Differential Entropy

```python
from metainformant.information import differential_entropy
import numpy as np

# Calculate differential entropy for continuous data
samples = np.random.normal(0, 1, 1000)
h = differential_entropy(samples, method="histogram")
```

#### Continuous Mutual Information

```python
from metainformant.information import mutual_information_continuous

x = np.random.randn(1000)
y = x + np.random.randn(1000) * 0.1  # Strong correlation
mi = mutual_information_continuous(x, y)
```

### Estimation Methods (`estimation.py`)

Entropy and mutual information estimation with bias correction.

#### Entropy Estimation

```python
from metainformant.information import entropy_estimator

counts = {"A": 50, "T": 30, "G": 20}
# Plug-in estimator
h_plugin = entropy_estimator(counts, method="plugin")
# Miller-Madow bias correction
h_mm = entropy_estimator(counts, method="miller_madow", bias_correction=True)
```

#### Bias Correction

```python
from metainformant.information import bias_correction

# Correct entropy estimate for small sample bias
corrected = bias_correction(estimate=1.5, n_samples=100, n_bins=10, measure="entropy")
```

### Additional Syntactic Methods

#### Jensen-Shannon Divergence

```python
from metainformant.information import jensen_shannon_divergence

p = [0.5, 0.3, 0.2]
q = [0.4, 0.4, 0.2]
js = jensen_shannon_divergence(p, q)  # Symmetrized KL divergence
```

#### Rényi Entropy

```python
from metainformant.information import renyi_entropy

probs = [0.25, 0.25, 0.25, 0.25]
# Collision entropy (α=2)
h_2 = renyi_entropy(probs, alpha=2.0)
# α=1 equals Shannon entropy
h_1 = renyi_entropy(probs, alpha=1.0)
```

#### Tsallis Entropy

```python
from metainformant.information import tsallis_entropy

probs = [0.5, 0.5]
h_q = tsallis_entropy(probs, q=2.0)  # Non-extensive entropy
```

#### Normalized Mutual Information

```python
from metainformant.information import normalized_mutual_information

x = [0, 1, 0, 1]
y = [0, 1, 0, 1]
nmi = normalized_mutual_information(x, y, method="arithmetic")  # Normalized to [0, 1]
```

### Workflow Functions (`workflows.py`)

High-level batch processing and workflow functions.

#### Batch Entropy Analysis

```python
from metainformant.information import batch_entropy_analysis

sequences = ["ATCG", "ATCG", "AAAA", "GCTA"]
results = batch_entropy_analysis(sequences, k=2, output_dir="output/information")
```

#### Information Workflow

```python
from metainformant.information import information_workflow

sequences = ["ATCGATCG", "AAAA", "GCTAGCTA"]
results = information_workflow(sequences, k_values=[1, 2, 3], output_dir="output/information")
```

#### Dataset Comparison

```python
from metainformant.information import compare_datasets

dataset1 = ["ATCG", "ATCG"]
dataset2 = ["AAAA", "AAAA"]
comparison = compare_datasets(dataset1, dataset2, k=1, method="kl_divergence")
```

### Integration Functions (`integration.py`)

Wrappers for integrating with other modules.

#### DNA Integration

```python
from metainformant.information import dna_integration
from metainformant.dna import sequences

dna_seqs = sequences.read_fasta("data/sequences.fasta")
results = dna_integration(dna_seqs, k=2, analysis_type="profile")
```

#### RNA Integration

```python
from metainformant.information import rna_integration
import numpy as np

expression = np.random.randn(100, 50)  # 100 samples, 50 genes
results = rna_integration(expression, method="mutual_information")
```

#### Single-Cell Integration

```python
from metainformant.information import singlecell_integration

count_matrix = np.random.randint(0, 100, (100, 50))
cell_types = ["TypeA"] * 50 + ["TypeB"] * 50
results = singlecell_integration(count_matrix, cell_types=cell_types, method="cell_type_entropy")
```

#### Multi-Omics Integration

```python
from metainformant.information import multiomics_integration

genomics = np.random.randn(100, 50)
transcriptomics = np.random.randn(100, 50)
results = multiomics_integration(
    genomics_data=genomics,
    transcriptomics_data=transcriptomics,
    method="cross_platform_mi"
)
```

#### ML Integration

```python
from metainformant.information import ml_integration

X = np.random.randn(100, 50)
y = np.random.randint(0, 2, 100)
results = ml_integration(X, y, method="feature_mi")  # Feature selection using MI
```

### Enhanced Visualization

#### Rényi Spectrum

```python
from metainformant.information.visualization import plot_renyi_spectrum

probs = [0.25, 0.25, 0.25, 0.25]
plot_renyi_spectrum(probs, output_path="output/information/renyi_spectrum.png")
```

#### Information Network

```python
from metainformant.information.visualization import plot_information_network
import networkx as nx

G = nx.karate_club_graph()
mi_matrix = np.random.rand(len(G.nodes()), len(G.nodes()))
plot_information_network(G, mi_matrix=mi_matrix, output_path="output/information/network.png")
```

#### Entropy Landscape

```python
from metainformant.information.visualization import plot_entropy_landscape

entropy_data = np.random.rand(10, 10)
plot_entropy_landscape(entropy_data, output_path="output/information/landscape.png")
```

## API Reference

### Syntactic Information Theory

#### `shannon_entropy(probs: Sequence[float], base: float = 2.0) -> float`
Calculate Shannon entropy of a probability distribution.

#### `mutual_information(x: Sequence[Any], y: Sequence[Any], base: float = 2.0) -> float`
Calculate mutual information between two sequences.

#### `kl_divergence(p: Sequence[float], q: Sequence[float], base: float = 2.0) -> float`
Calculate Kullback-Leibler divergence.

#### `jensen_shannon_divergence(p: Sequence[float], q: Sequence[float], base: float = 2.0) -> float`
Calculate Jensen-Shannon divergence (symmetrized KL).

#### `renyi_entropy(probs: Sequence[float], alpha: float = 2.0, base: float = 2.0) -> float`
Calculate Rényi entropy of order α.

#### `tsallis_entropy(probs: Sequence[float], q: float = 2.0, base: float = 2.0) -> float`
Calculate Tsallis entropy (non-extensive).

#### `normalized_mutual_information(x: Sequence[Any], y: Sequence[Any], method: str = "arithmetic", base: float = 2.0) -> float`
Calculate normalized mutual information.

#### `transfer_entropy(x: Sequence[Any], y: Sequence[Any], lag: int = 1, base: float = 2.0) -> float`
Calculate transfer entropy for time series.

### Continuous Information Theory

#### `differential_entropy(samples: np.ndarray, method: str = "histogram", bins: int | None = None) -> float`
Calculate differential entropy for continuous data.

#### `mutual_information_continuous(x: np.ndarray, y: np.ndarray, method: str = "histogram", bins: int | None = None) -> float`
Calculate mutual information for continuous variables.

### Estimation Methods

#### `entropy_estimator(counts: dict[Any, int] | list[int], method: str = "plugin", bias_correction: bool = True) -> float`
Estimate entropy from counts with various methods.

#### `mutual_information_estimator(x: list[Any], y: list[Any], method: str = "plugin", bias_correction: bool = True) -> float`
Estimate mutual information with bias correction.

### Workflow Functions

#### `batch_entropy_analysis(sequences: list[str], k: int = 1, output_dir: Path | None = None) -> dict[str, Any]`
Perform entropy analysis on multiple sequences.

#### `information_workflow(sequences: list[str], k_values: list[int] | None = None, output_dir: Path | None = None) -> dict[str, Any]`
Complete information-theoretic workflow for sequence analysis.

## Performance Considerations

### Complexity

- **Shannon Entropy**: O(n) where n is number of symbols
- **Mutual Information**: O(n) where n is sequence length
- **KL Divergence**: O(n) where n is distribution size
- **Rényi/Tsallis Entropy**: O(n) where n is distribution size
- **Transfer Entropy**: O(n × m) where n is sequence length, m is lag
- **Network Entropy**: O(V) where V is number of nodes

### Memory Usage

- Most functions have O(n) memory complexity
- Large sequence analysis: Consider processing in batches
- Network analysis: Memory scales with number of edges
- Continuous methods: Histogram-based methods use O(bins) memory

### Best Practices

1. **Discretization**: For continuous data, choose appropriate bin counts
   - Too few bins: Loss of information
   - Too many bins: Overfitting and bias

2. **Bias Correction**: Always use bias correction for small samples
   - Sample size < 100: Use Miller-Madow or Chao-Shen
   - Sample size > 1000: Plug-in estimator usually sufficient

3. **K-mer Size**: For sequence analysis
   - k=1: Base composition
   - k=2-3: Short motifs
   - k>3: Long-range patterns (requires more data)

4. **Normalization**: Use normalized MI for comparing different datasets

5. **Continuous Data**: Prefer histogram method for most cases, KDE for smooth distributions

## Troubleshooting

### Common Issues

**Issue**: Entropy returns 0.0 for non-empty data
- **Cause**: All values are identical
- **Solution**: Check data diversity

**Issue**: MI returns negative values
- **Cause**: Numerical precision issues (should be non-negative)
- **Solution**: Use `max(0.0, mi)` or check for bugs

**Issue**: KL divergence returns infinity
- **Cause**: Q has zeros where P has positives
- **Solution**: Smooth distributions or check data alignment

**Issue**: Continuous entropy estimation fails
- **Cause**: Insufficient samples or inappropriate binning
- **Solution**: Increase sample size or adjust bin count

**Issue**: Network functions fail with ImportError
- **Cause**: NetworkX not installed
- **Solution**: Install NetworkX: `pip install networkx`

## Integration Tutorials

### Tutorial 1: DNA Sequence Analysis Workflow

```python
from metainformant.dna import sequences
from metainformant.information import information_workflow, information_report

# Load sequences
dna_seqs = sequences.read_fasta("data/sequences.fasta")
seq_list = list(dna_seqs.values())

# Run complete workflow
results = information_workflow(seq_list, k_values=[1, 2, 3], output_dir="output/information")

# Generate report
information_report(results, output_path="output/information/report.md", format="markdown")
```

### Tutorial 2: RNA Expression Analysis

```python
from metainformant.information import rna_integration
import numpy as np

# Load expression data (samples x genes)
expression = np.load("data/expression.npy")

# Calculate entropy for each gene
results = rna_integration(expression, method="entropy")

# Find high-entropy genes (diverse expression)
high_entropy_genes = [
    g for g in results["gene_entropies"] if g["entropy"] > results["mean_entropy"]
]
```

### Tutorial 3: Feature Selection using Information Theory

```python
from metainformant.information import ml_integration
import numpy as np

# Load features and labels
X = np.load("data/features.npy")
y = np.load("data/labels.npy")

# Select features using mutual information
results = ml_integration(X, y, method="feature_mi")

# Get top features
top_features = results["top_features"]
X_selected = X[:, top_features]
```

### Tutorial 4: Multi-Omics Information Integration

```python
from metainformant.information import multiomics_integration
import numpy as np

# Load multi-omics data
genomics = np.load("data/genomics.npy")
transcriptomics = np.load("data/transcriptomics.npy")
proteomics = np.load("data/proteomics.npy")

# Cross-platform information analysis
results = multiomics_integration(
    genomics_data=genomics,
    transcriptomics_data=transcriptomics,
    proteomics_data=proteomics,
    method="platform_entropy"
)

# Compare information content across platforms
print(f"Genomics entropy: {results['genomics_entropy']:.3f}")
print(f"Transcriptomics entropy: {results['transcriptomics_entropy']:.3f}")
print(f"Proteomics entropy: {results['proteomics_entropy']:.3f}")
```

### Tutorial 5: Network Information Analysis

```python
from metainformant.networks import create_network
from metainformant.information.networks import network_entropy, mutual_information_network
from metainformant.information.visualization import plot_mi_network

# Create or load network
network = create_network(["A", "B", "C", "D"], directed=False)
network.add_edge("A", "B")
network.add_edge("B", "C")
network.add_edge("C", "D")

# Calculate network entropy
entropy = network_entropy(network)

# Calculate MI matrix
mi_matrix = mutual_information_network(network)

# Visualize
plot_mi_network(mi_matrix, labels=["A", "B", "C", "D"], output_path="output/information/mi_net.png")
```

## Integration with Other Modules

### With RNA Module

```python
from metainformant.information import rna_integration
from metainformant.rna import workflow

# Analyze expression data entropy
expression_data = workflow.extract_expression_data("data/rna_counts.csv")
results = rna_integration(expression_data, method="mutual_information")
```

### With Single-Cell Module

```python
from metainformant.information import singlecell_integration
from metainformant.singlecell import preprocessing

# Load and preprocess single-cell data
count_matrix = preprocessing.load_count_matrix("data/counts.h5ad")
cell_types = preprocessing.get_cell_types(count_matrix)

# Information analysis
results = singlecell_integration(count_matrix, cell_types=cell_types, method="cell_type_entropy")
```

### With Multi-Omics Module

```python
from metainformant.information import multiomics_integration
from metainformant.multiomics import integration

# Load multi-omics dataset
omics_data = integration.load_multiomics_data("data/omics_data.h5")
results = multiomics_integration(
    genomics_data=omics_data["genomics"],
    transcriptomics_data=omics_data["transcriptomics"],
    method="cross_platform_mi"
)
```

### With ML Module

```python
from metainformant.information import ml_integration
from metainformant.ml import features

# Feature selection using information theory
X, y = features.load_dataset("data/ml_data.csv")
results = ml_integration(X, y, method="feature_mi")

# Use top features for model training
top_features = results["top_features"]
X_selected = X[:, top_features]
```

### With Quality Module

```python
from metainformant.quality import fastq
from metainformant.information import shannon_entropy_from_counts

# Analyze quality score entropy
quality_scores = fastq.read_quality_scores("data/reads.fastq")
from collections import Counter
counts = Counter(quality_scores.flatten())
entropy = shannon_entropy_from_counts(counts)
```

### With Simulation Module

```python
from metainformant.simulation import sequences
from metainformant.information import information_profile

# Analyze simulated sequence entropy
sim_seqs = sequences.generate_dna_sequences(n=100, length=1000)
profile = information_profile(sim_seqs, k=3)

# Compare with real data
real_profile = information_profile(real_sequences, k=3)
```

## References

- Shannon, C. E. (1948). A mathematical theory of communication. Bell System Technical Journal, 27(3), 379-423.
- Cover, T. M., & Thomas, J. A. (2006). Elements of Information Theory. John Wiley & Sons.
- Schreiber, T. (2000). Measuring information transfer. Physical Review Letters, 85(2), 461.
- Kullback, S., & Leibler, R. A. (1951). On information and sufficiency. The Annals of Mathematical Statistics, 22(1), 79-86.
- Rényi, A. (1961). On measures of entropy and information. Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics and Probability, 1, 547-561.
- Tsallis, C. (1988). Possible generalization of Boltzmann-Gibbs statistics. Journal of Statistical Physics, 52(1-2), 479-487.
- Lin, J. (1991). Divergence measures based on the Shannon entropy. IEEE Transactions on Information Theory, 37(1), 145-151.

This module provides comprehensive information-theoretic tools for analyzing biological data across multiple scales and domains.

