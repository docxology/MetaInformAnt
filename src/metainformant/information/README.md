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

## References

- Shannon, C. E. (1948). A mathematical theory of communication. Bell System Technical Journal, 27(3), 379-423.
- Cover, T. M., & Thomas, J. A. (2006). Elements of Information Theory. John Wiley & Sons.
- Schreiber, T. (2000). Measuring information transfer. Physical Review Letters, 85(2), 461.
- Kullback, S., & Leibler, R. A. (1951). On information and sufficiency. The Annals of Mathematical Statistics, 22(1), 79-86.

This module provides comprehensive information-theoretic tools for analyzing biological data across multiple scales and domains.

