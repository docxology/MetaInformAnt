# Information Theory Workflows

Step-by-step workflow documentation for common information-theoretic analysis tasks.

All import paths and result keys below match the current `metainformant.information` API.

## Workflow 1: Sequence Information Analysis

### Objective
Analyze information content and complexity of biological sequences.

### Steps

1. **Load Sequences**
```python
from metainformant.dna.sequence.core import read_fasta
dna_seqs = read_fasta("data/sequences.fasta")
```

2. **Calculate Information Profile**
```python
from metainformant.information.metrics.analysis import information_profile
profile = information_profile(list(dna_seqs.values()), k=2)
print(f"Mean entropy: {profile['statistics']['mean_entropy']:.3f} bits")
```

3. **Analyze Individual Sequences**
```python
from metainformant.information.metrics.analysis import analyze_sequence_information
for seq_id, seq in dna_seqs.items():
    analysis = analyze_sequence_information(seq, k_values=[1, 2, 3])
    entropy = analysis["k_mer_analysis"]["k1"]["entropy"]
    print(f"{seq_id}: 1-mer entropy={entropy:.3f}")
```

## Workflow 2: Network Information Analysis

### Objective
Analyze information flow and structure in biological networks.

### Steps

1. **Build a Network**
```python
import networkx as nx
G = nx.Graph()
G.add_edges_from([("A", "B"), ("B", "C"), ("C", "D")])
```

2. **Calculate Network Entropy**
```python
from metainformant.information.integration.networks import network_entropy
entropy = network_entropy(G)
```

3. **Analyze Information Flow**
```python
from metainformant.information.integration.networks import information_flow
flow = information_flow(G, source_nodes=["A"], target_nodes=["D"], steps=50)
print(f"Path length entropy: {flow['path_length_entropy']:.3f} bits")
```

4. **Build an Information Flow Network from Time Series**
```python
from metainformant.information.network_info.information_flow import information_flow_network
time_series = {"A": [...], "B": [...], "C": [...], "D": [...]}  # equal-length series
result = information_flow_network(time_series, method="transfer_entropy", threshold=0.05)
print(result["edges"])  # significant directed edges
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
from metainformant.information.integration import multiomics_integration
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
from metainformant.information.integration import ml_integration
results = ml_integration(X, y, method="feature_mi")
```

3. **Select Top Features**
```python
top_features = results["top_features"]  # list of {"index": int, "mi": float}
top_indices = [f["index"] for f in top_features[:50]]  # top 50 features
X_selected = X[:, top_indices]
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
from metainformant.dna.sequence.core import read_fasta
from metainformant.information.workflow import batch_entropy_analysis

all_results = {}
for fasta_file in fasta_files:
    seqs = read_fasta(str(fasta_file))
    results = batch_entropy_analysis(list(seqs.values()), k=2)
    all_results[fasta_file.stem] = results
```

3. **Generate Report**
```python
from metainformant.information.workflow import information_report
information_report(
    results,
    output_path="docs/information/batch_report.md",
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
from metainformant.information.metrics.core.continuous import differential_entropy
h = differential_entropy(samples, method="histogram", bins=20)
```

3. **Compare with Estimation Methods**
```python
from metainformant.information.metrics.core.continuous import entropy_estimation
h_hist = entropy_estimation(samples, method="histogram")
h_knn = entropy_estimation(samples, method="knn")
print(f"Histogram: {h_hist:.3f}, k-NN: {h_knn:.3f}")
```

Note: the histogram/kde/knn methods estimate *differential* (continuous) entropy.
Discrete estimators such as plugin/Miller-Madow live in
`metainformant.information.metrics.core.estimation` and operate on counts.

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
from metainformant.information.metrics.advanced.semantic import (
    information_content_from_annotations,
)
# Per-term IC: returns -log2(fraction of genes annotated with the term)
ic_bp = information_content_from_annotations(gene_annotations, "GO:0008150")
```

3. **Calculate Semantic Similarity**
```python
from metainformant.information.metrics.advanced.semantic import semantic_similarity_matrix
from metainformant.information.metrics.advanced.semantic import information_content

# Build per-term IC values and a term -> parent-terms hierarchy
term_frequencies = {"GO:0008150": 2, "GO:0003674": 1, "GO:0005524": 1}
term_ic = {t: information_content(term_frequencies, t) for t in term_frequencies}
hierarchy = {"GO:0003674": {"GO:0008150"}, "GO:0005524": {"GO:0008150"}}

terms = list(term_ic)
similarity_matrix = semantic_similarity_matrix(terms, term_ic, hierarchy)
```

## Workflow 8: Complete Information Analysis Pipeline

### Objective
End-to-end information-theoretic analysis from raw data to report.

### Steps

1. **Data Loading and Preprocessing**
```python
from metainformant.dna.sequence.core import read_fasta
from metainformant.information.workflow import information_workflow

seqs = read_fasta("data/sequences.fasta")
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
from metainformant.information.workflow import information_report
information_report(
    results,
    output_path="docs/information/full_report.md",
    format="markdown"
)
```

4. **Inspect Aggregate Results**
```python
print(results["aggregate_results"]["k_mer_statistics"]["k1"]["mean_entropy"])
print(results["workflow_status"], f"{results['processing_time']:.2f}s")
```

These workflows provide step-by-step guides for common information-theoretic analysis tasks in biological data.
