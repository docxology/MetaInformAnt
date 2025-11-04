# Information Theory Workflows

Step-by-step workflow documentation for common information-theoretic analysis tasks.

## Workflow 1: Sequence Information Analysis

### Objective
Analyze information content and complexity of biological sequences.

### Steps

1. **Load Sequences**
```python
from metainformant.dna import sequences
dna_seqs = sequences.read_fasta("data/sequences.fasta")
```

2. **Calculate Information Profile**
```python
from metainformant.information import information_profile
profile = information_profile(list(dna_seqs.values()), k=2)
```

3. **Analyze Individual Sequences**
```python
from metainformant.information import analyze_sequence_information
for seq_id, seq in dna_seqs.items():
    analysis = analyze_sequence_information(seq, k_values=[1, 2, 3])
    print(f"{seq_id}: entropy={analysis['kmer_analyses'][1]['entropy']:.3f}")
```

4. **Visualize Results**
```python
from metainformant.information.visualization import plot_information_profile
plot_information_profile(profile, output_path="output/information/profile.png")
```

## Workflow 2: Network Information Analysis

### Objective
Analyze information flow and structure in biological networks.

### Steps

1. **Load Network**
```python
from metainformant.networks import create_network
network = create_network(["A", "B", "C", "D"], directed=False)
network.add_edge("A", "B")
network.add_edge("B", "C")
```

2. **Calculate Network Entropy**
```python
from metainformant.information.networks import network_entropy
entropy = network_entropy(network)
```

3. **Analyze Information Flow**
```python
from metainformant.information.networks import information_flow
flow = information_flow(network, source_nodes=["A"], target_nodes=["D"])
```

4. **Calculate MI Matrix**
```python
from metainformant.information.networks import mutual_information_network
mi_matrix = mutual_information_network(network)
```

5. **Visualize Network**
```python
from metainformant.information.visualization import plot_mi_network
plot_mi_network(mi_matrix, labels=["A", "B", "C", "D"], output_path="output/information/network.png")
```

## Workflow 3: Multi-Omics Information Integration

### Objective
Integrate information measures across multiple omics platforms.

### Steps

1. **Load Multi-Omics Data**
```python
import numpy as np
genomics = np.load("data/genomics.npy")
transcriptomics = np.load("data/transcriptomics.npy")
proteomics = np.load("data/proteomics.npy")
```

2. **Calculate Platform Entropy**
```python
from metainformant.information import multiomics_integration
results = multiomics_integration(
    genomics_data=genomics,
    transcriptomics_data=transcriptomics,
    proteomics_data=proteomics,
    method="platform_entropy"
)
```

3. **Cross-Platform Analysis**
```python
mi_results = multiomics_integration(
    genomics_data=genomics,
    transcriptomics_data=transcriptomics,
    method="cross_platform_mi"
)
```

4. **Compare Information Content**
```python
print(f"Genomics entropy: {results['genomics_entropy']:.3f}")
print(f"Transcriptomics entropy: {results['transcriptomics_entropy']:.3f}")
print(f"Proteomics entropy: {results['proteomics_entropy']:.3f}")
```

## Workflow 4: Feature Selection using Information Theory

### Objective
Select informative features for machine learning using mutual information.

### Steps

1. **Load Data**
```python
import numpy as np
X = np.load("data/features.npy")
y = np.load("data/labels.npy")
```

2. **Calculate Feature MI**
```python
from metainformant.information import ml_integration
results = ml_integration(X, y, method="feature_mi")
```

3. **Select Top Features**
```python
top_features = results["top_features"][:50]  # Top 50 features
X_selected = X[:, top_features]
```

4. **Train Model with Selected Features**
```python
from sklearn.ensemble import RandomForestClassifier
model = RandomForestClassifier()
model.fit(X_selected, y)
```

## Workflow 5: Batch Processing Workflow

### Objective
Process multiple datasets and generate comprehensive reports.

### Steps

1. **Prepare Data**
```python
from pathlib import Path
fasta_files = list(Path("data").glob("*.fasta"))
```

2. **Batch Analysis**
```python
from metainformant.information import batch_entropy_analysis
from metainformant.dna import sequences

all_results = {}
for fasta_file in fasta_files:
    seqs = sequences.read_fasta(str(fasta_file))
    results = batch_entropy_analysis(list(seqs.values()), k=2)
    all_results[fasta_file.stem] = results
```

3. **Generate Report**
```python
from metainformant.information import information_report
information_report(
    all_results,
    output_path="output/information/batch_report.md",
    format="markdown"
)
```

## Workflow 6: Continuous Data Analysis Workflow

### Objective
Analyze continuous biological data using differential entropy.

### Steps

1. **Load Continuous Data**
```python
import numpy as np
samples = np.random.normal(0, 1, 1000)
```

2. **Estimate Differential Entropy**
```python
from metainformant.information import differential_entropy
h = differential_entropy(samples, method="histogram", bins=20)
```

3. **Compare with Estimation Methods**
```python
from metainformant.information.continuous import entropy_estimation
h_plugin = entropy_estimation(samples, method="plugin")
h_mm = entropy_estimation(samples, method="miller_madow")
print(f"Plugin: {h_plugin:.3f}, Miller-Madow: {h_mm:.3f}")
```

## Workflow 7: Semantic Information Analysis

### Objective
Analyze semantic similarity and information content for ontological data.

### Steps

1. **Load Annotations**
```python
gene_annotations = {
    "gene1": {"GO:0008150", "GO:0003674"},
    "gene2": {"GO:0008150", "GO:0005524"},
}
```

2. **Calculate Information Content**
```python
from metainformant.information import information_content_from_annotations
term_ic = information_content_from_annotations(gene_annotations)
```

3. **Calculate Semantic Similarity**
```python
from metainformant.information import semantic_similarity_matrix
terms = list(term_ic.keys())
similarity_matrix = semantic_similarity_matrix(terms, term_ic)
```

4. **Visualize Similarity Network**
```python
from metainformant.information.visualization import plot_semantic_similarity_network
plot_semantic_similarity_network(
    similarity_matrix,
    terms,
    output_path="output/information/semantic_network.png"
)
```

## Workflow 8: Complete Information Analysis Pipeline

### Objective
End-to-end information-theoretic analysis from raw data to visualization.

### Steps

1. **Data Loading and Preprocessing**
```python
from metainformant.dna import sequences
from metainformant.information import information_workflow

seqs = sequences.read_fasta("data/sequences.fasta")
seq_list = list(seqs.values())
```

2. **Comprehensive Analysis**
```python
results = information_workflow(
    seq_list,
    k_values=[1, 2, 3],
    output_dir="output/information"
)
```

3. **Generate Report**
```python
from metainformant.information import information_report
information_report(
    results,
    output_path="output/information/full_report.md",
    format="markdown"
)
```

4. **Visualization**
```python
from metainformant.information.visualization import (
    plot_information_profile,
    plot_entropy_distribution,
)

# Plot profile
profile = results["profiles"][1]
plot_information_profile(profile, output_path="output/information/profile.png")

# Plot entropy distribution
entropies = [r["entropy"] for r in results.get("sequences", [])]
if entropies:
    plot_entropy_distribution(entropies, output_path="output/information/entropy_dist.png")
```

These workflows provide step-by-step guides for common information-theoretic analysis tasks in biological data.

