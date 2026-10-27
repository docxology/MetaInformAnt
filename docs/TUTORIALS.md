# METAINFORMANT Tutorials

Complete end-to-end tutorials for using METAINFORMANT across all biological domains. Each tutorial includes data preparation, analysis workflow, and result interpretation.

## Quick Start Tutorials

### 1. DNA Sequence Analysis Tutorial

Learn basic DNA sequence analysis including composition, motifs, and alignments.

#### Data Preparation

```python
# Create sample DNA sequences
sequences = {
    "gene1": "ATCGATCGATCGATCGATCG",
    "gene2": "GCTAGCTAGCTAGCTAGCTA",
    "gene3": "ATATATATATATATATATAT"
}
```

#### Basic Sequence Analysis

```python
from metainformant.dna.sequence import composition, core as sequence_core

# Analyze sequence composition
for name, seq in sequences.items():
    gc_content = composition.gc_content(seq)
    length = sequence_core.sequence_length(seq)
    print(f"{name}: GC={gc_content:.2f}, Length={length}")

# Find motifs
motifs = ["ATCG", "GCTA"]
for name, seq in sequences.items():
    for motif in motifs:
        positions = sequence_core.find_motifs(seq, [motif])
        print(f"{name} - {motif}: {positions}")
```

#### Advanced Analysis

```python-snippet
# Calculate pairwise distances
from metainformant.dna.alignment import distances

seq_list = list(sequences.values())
distance_matrix = []
for i, seq1 in enumerate(seq_list):
    row = []
    for j, seq2 in enumerate(seq_list):
        if i == j:
            row.append(0.0)
        else:
            dist = distances.jukes_cantor_distance(seq1, seq2)
            row.append(dist)
    distance_matrix.append(row)

print("Distance matrix:")
for i, row in enumerate(distance_matrix):
    print(f"Seq{i+1}: {row}")
```

### 2. Population Genetics Analysis Tutorial

Analyze genetic variation in populations using coalescent theory and neutrality tests.

#### Sample Data

```python
# Simulated population data
sequences = [
    "ATCGATCGATCG",  # Sample 1
    "ATCGATCGATCG",  # Sample 2
    "ATCGATCGATCG",  # Sample 3
    "ATCGATCGATCG",  # Sample 4
    "ATCGATCGATCG",  # Sample 5
]
```

#### Coalescent Analysis

```python-snippet
from metainformant.math.population_genetics import coalescent

# Calculate coalescent parameters
n_samples = len(sequences)
effective_size = 10000
mutation_rate = 0.0001

theta = coalescent.expected_pairwise_diversity(effective_size, mutation_rate)  # 4·Ne·μ
total_length = coalescent.expected_total_branch_length(n_samples, effective_size)
time_to_mrca = coalescent.expected_time_to_mrca(n_samples, effective_size)

print(f"Population mutation rate (θ): {theta:.2f}")
print(f"Expected total branch length: {total_length:.2f}")
print(f"Expected time to MRCA: {time_to_mrca:.2f}")

# Calculate Tajima's D (would use real sequence data):
# mean pairwise differences (pi), segregating sites (S) and sample size (n)
# tajima_d = coalescent.tajimas_D(pi, S, n_samples)
```

#### Neutrality Tests

```python
from metainformant.dna import population

# Calculate population genetics statistics
pi_diversity = population.nucleotide_diversity(sequences)
segregating_sites = population.segregating_sites(sequences)

print(f"Nucleotide diversity (π): {pi_diversity:.4f}")
print(f"Segregating sites (S): {segregating_sites}")

# Tajima's D test
tajima_d = population.tajimas_d(sequences)
print(f"Tajima's D: {tajima_d:.4f}")

# Interpret results
if tajima_d < -2.0:
    print("Evidence of positive selection or population expansion")
elif tajima_d > 2.0:
    print("Evidence of balancing selection or population contraction")
else:
    print("Neutral evolution")
```

### 3. GWAS Analysis Tutorial

Complete genome-wide association study workflow from variant data to results.

#### Data Preparation

```python
# Sample genotype and phenotype data
import numpy as np

# Simulate genotype matrix (100 samples, 50 SNPs)
n_samples, n_snps = 100, 50
genotypes = np.random.randint(0, 3, size=(n_samples, n_snps))  # 0, 1, 2

# Simulate phenotype (quantitative trait)
phenotype = genotypes[:, 0] * 0.5 + np.random.normal(0, 0.1, n_samples)  # SNP 0 has effect
```

#### GWAS Analysis

```python-snippet
from metainformant.gwas.analysis import association

# Perform association tests
results = []
for snp_idx in range(n_snps):
    genotypes_snp = genotypes[:, snp_idx]

    # Linear regression association test
    result = association.association_test_linear(
        genotypes=genotypes_snp.tolist(),
        phenotypes=phenotype.tolist(),
        covariates=None  # No covariates in this simple example
    )

    results.append({
        "snp": f"SNP_{snp_idx}",
        "p_value": result["p_value"],
        "beta": result["beta"],
        "se": result["se"]
    })

# Sort by p-value
results.sort(key=lambda x: x["p_value"])

# Show top associations
print("Top GWAS associations:")
for i, result in enumerate(results[:5]):
    print(f"  {result['snp']}: p={result['p_value']:.2e}, beta={result['beta']:.3f}")
```

#### Visualization

```python
from metainformant.core import io
from metainformant.gwas import visualization

# Prepare data for Manhattan plot (records use chrom/pos/p_value keys)
manhattan_data = [
    {
        "chrom": "1",
        "pos": i * 1000,  # Position in bp
        "p_value": result["p_value"],
        "variant_id": result["snp"]
    }
    for i, result in enumerate(results)
]

# Create Manhattan plot (the output directory must exist; the function returns a matplotlib Figure)
io.ensure_directory("output/gwas")
output_path = "output/gwas/manhattan.png"
fig = visualization.manhattan_plot(manhattan_data, output_path)
print(f"Manhattan plot saved: {output_path}")
```

### 4. RNA-seq Analysis Tutorial

Complete RNA-seq analysis pipeline from raw reads to gene expression quantification.

#### Configuration Setup

```python
from metainformant.rna.engine.workflow import AmalgkitWorkflowConfig
from pathlib import Path

# Configure RNA-seq analysis
config = AmalgkitWorkflowConfig(
    work_dir=Path("output/rna_analysis"),
    threads=4,
    species_list=["Drosophila_melanogaster"]
)

print("RNA-seq workflow configured:")
print(f"Working directory: {config.work_dir}")
print(f"Species: {config.species_list}")
print(f"Threads: {config.threads}")
```

#### Quality Control

```python
from pathlib import Path

from metainformant.quality.io import fastq

# Analyze FASTQ quality (assuming you have FASTQ files)
fastq_files = ["data/rna/sample1_R1.fastq.gz", "data/rna/sample1_R2.fastq.gz"]

for fastq_file in fastq_files:
    if not Path(fastq_file).exists():
        print(f"FASTQ file not found: {fastq_file}")
        continue

    quality_report = fastq.analyze_fastq_quality(fastq_file)
    stats = quality_report["basic_statistics"]
    print(f"FASTQ Quality for {fastq_file}:")
    print(f"  Total reads: {stats['total_reads']}")
    print(f"  Mean quality: {stats['mean_quality']:.2f}")
    print(f"  GC content: {stats['mean_gc']:.2f}%")
```

#### Expression Quantification

```python-snippet
from metainformant.rna.amalgkit import quant

# Quantify expression with the amalgkit quant step (kallisto on the downloaded reads).
# Needs the amalgkit CLI and prior metadata/getfastq steps, so this is a demonstration of the API.
# Further options include index_dir, metadata and quant_backend.
quant_params = {
    "out_dir": "output/rna_analysis/amalgkit",
    "threads": 4,
}

# result = quant(quant_params)  # returns a subprocess.CompletedProcess
print("Gene expression quantification would run here")
print("Output: per-sample quantification files under out_dir")
```

### 5. Information Theory Analysis Tutorial

Apply information-theoretic methods to analyze sequence complexity and patterns.

#### Sequence Complexity Analysis

```python-snippet
from metainformant.information.metrics.analysis import analysis

# Analyze sequence information content
test_sequences = [
    "AAAAAAAAAAAAAAAAAAAAAAAA",  # Low complexity
    "ATCGATCGATCGATCGATCGATCG",  # Periodic
    "ATCGTAGCTAGCTAGCTAGCTAGC",  # Random-like
]

print("Sequence Information Analysis:")
for i, seq in enumerate(test_sequences):
    info_result = analysis.analyze_sequence_information(seq, k_values=[1, 2, 3])

    print(f"\nSequence {i+1}: {seq[:20]}...")
    for k in [1, 2, 3]:
        if f"k{k}" in info_result["k_mer_analysis"]:
            analysis_k = info_result["k_mer_analysis"][f"k{k}"]
            print(f"  k={k}: entropy={analysis_k['entropy']:.3f}, "
                  f"unique_kmers={analysis_k['n_unique']}")
```

#### Batch Processing

```python-snippet
from metainformant.information.workflow import workflows

# Batch entropy analysis
sequences = [
    "ATCGATCGATCG",
    "GCTAGCTAGCTA",
    "ATATATATATAT",
    "CGCGCGCGCGCG"
]

# Methods default to ["plugin", "miller_madow"]
batch_results = workflows.batch_entropy_analysis(sequences, k=2)

print("\nBatch Entropy Analysis Results:")
print(f"Total sequences analyzed: {batch_results['batch_info']['n_sequences']}")
plugin_stats = batch_results["summary"]["entropy_statistics"]["plugin"]
print(f"Mean entropy: {plugin_stats['mean']:.3f}")
print(f"Entropy range: {plugin_stats['min']:.3f} - {plugin_stats['max']:.3f}")

# Detailed results per sequence
for seq_result in batch_results["sequence_results"]:
    plugin = seq_result["methods"]["plugin"]
    print(f"Sequence {seq_result['sequence_index']}: "
          f"entropy={plugin['entropy']:.3f}, "
          f"length={seq_result['sequence_length']}")
```

### 6. Multi-Omics Integration Tutorial

Combine data from multiple biological layers for comprehensive analysis.

#### Data Integration Setup

```python
from metainformant import multiomics

# Example multi-omics datasets
# In practice, these would be real biological data

# Simulated gene expression data
expression_data = {
    "sample1": {"GENE1": 10.5, "GENE2": 8.2, "GENE3": 12.1},
    "sample2": {"GENE1": 9.8, "GENE2": 8.9, "GENE3": 11.5},
    "sample3": {"GENE1": 11.2, "GENE2": 7.5, "GENE3": 12.8}
}

# Simulated methylation data
methylation_data = {
    "sample1": {"CGI1": 0.85, "CGI2": 0.72},
    "sample2": {"CGI1": 0.82, "CGI2": 0.68},
    "sample3": {"CGI1": 0.88, "CGI2": 0.75}
}

# Simulated protein data
protein_data = {
    "sample1": {"PROT1": 15.2, "PROT2": 8.7},
    "sample2": {"PROT1": 14.8, "PROT2": 9.1},
    "sample3": {"PROT1": 16.1, "PROT2": 8.4}
}
```

#### Integration Analysis

```python-snippet
import pandas as pd
from metainformant.multiomics.analysis import integration

# Each layer is a samples x features DataFrame
layers = {
    "expression": pd.DataFrame.from_dict(expression_data, orient="index"),
    "methylation": pd.DataFrame.from_dict(methylation_data, orient="index"),
    "proteomics": pd.DataFrame.from_dict(protein_data, orient="index"),
}

# Combine multi-omics data (only samples present in every layer are kept)
integrated_data = integration.integrate_omics_data(data=layers)

# Perform joint PCA: returns (sample embeddings, per-layer loadings, explained variance ratio)
embeddings, loadings, explained_variance = integration.joint_pca(integrated_data, n_components=3)

print("Multi-omics Integration Results:")
print(f"Joint embedding shape: {embeddings.shape}")
print("Available data types:", list(loadings.keys()))
```

### 7. Machine Learning Tutorial

Apply machine learning methods to biological classification and prediction tasks.

#### Classification Example

```python
from metainformant.ml.features import features
from metainformant.ml.models import classification
import numpy as np

# Create sample biological features and labels
np.random.seed(42)
n_samples = 200
n_features = 50

# Simulated gene expression features
X = np.random.randn(n_samples, n_features)

# Binary classification labels (e.g., disease vs healthy)
y = np.random.randint(0, 2, n_samples)

print(f"Dataset: {n_samples} samples, {n_features} features")
print(f"Class distribution: {np.bincount(y)}")
```

#### Model Training and Evaluation

```python
# Split data for training/testing
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42, stratify=y
)

# Train classifier (method: "rf", "gb" or "lr")
classifier = classification.create_biological_classifier("rf", random_state=42)
classifier.fit(X_train, y_train)

# Cross-validation (stratified k-fold; one entry per metric, each with finite per-fold scores)
cv_results = classification.cross_validate_biological(
    X_train, y_train, method="rf", cv_folds=5, random_state=42
)

print("Cross-validation Results:")
for metric, scores in cv_results["cross_validation"].items():
    print(f"  {metric}: {scores['mean']:.3f} (+/- {scores['std']:.3f})")

# Test predictions
predictions = classifier.predict(X_test)
probabilities = classifier.predict_proba(X_test)

accuracy = np.mean(predictions == y_test)
print(f"Test accuracy: {accuracy:.3f}")
```

`cross_validate_biological` reports `accuracy`, `precision_weighted`, `recall_weighted` and `f1_weighted`, plus `roc_auc` when there are two classes; each entry holds `mean`, `std` and the per-fold `scores`. A metric that cannot be computed (for example `cv_folds` larger than the sample count of every class) raises instead of producing a `None` or NaN entry, so the loop above never meets a missing metric.

To compare several methods, `classification.compare_classifiers(X_train, y_train, methods=["rf", "gb", "lr"], cv_folds=5)` returns the per-method results under `comparison` and the winner under `best_method`. A method that fails is reported as `{"error": message}` and is never selected; `best_method` is `None` when every method failed.

#### Feature Selection

```python
# Select important features (returns the reduced matrix and the selected column indices)
X_selected, selected_indices = features.select_features_univariate(
    X_train, y_train, method="mutual_info_classif", k=20
)

print(f"Selected {len(selected_indices)} features out of {n_features}")
print("Feature selection method: mutual information")
```

### 8. Network Analysis Tutorial

Analyze biological networks including protein-protein interactions and pathways.

#### Network Construction

```python
from metainformant.networks.analysis import community, graph

# Create sample protein interaction network
edges = [
    ("PROT1", "PROT2"),
    ("PROT1", "PROT3"),
    ("PROT2", "PROT3"),
    ("PROT2", "PROT4"),
    ("PROT3", "PROT5"),
    ("PROT4", "PROT5"),
]

# Build network
G = graph.create_network(edges, directed=False)

print("Protein Interaction Network:")
print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")
print(f"Connected components: {len(graph.get_connected_components(G))}")
```

#### Network Analysis

```python
# Calculate centrality measures (a dict of {measure: {node: value}})
centralities = graph.centrality_measures(G)
degree_cent = centralities["degree"]
betweenness_cent = centralities["betweenness"]

print("\nTop nodes by degree centrality:")
for node, centrality in sorted(degree_cent.items(), key=lambda x: x[1], reverse=True)[:3]:
    print(f"  {node}: {centrality:.3f}")

print("\nTop nodes by betweenness centrality:")
for node, centrality in sorted(betweenness_cent.items(), key=lambda x: x[1], reverse=True)[:3]:
    print(f"  {node}: {centrality:.3f}")
```

#### Community Detection

```python
# Detect communities: returns {node: community_id}
# (Louvain needs the optional python-louvain package; without it detection falls back to greedy modularity)
node_communities = community.detect_communities(G, method="louvain")

communities = {}
for node, community_id in node_communities.items():
    communities.setdefault(community_id, []).append(node)

print(f"\nDetected {len(communities)} communities:")
for community_id, members in sorted(communities.items()):
    print(f"  Community {community_id + 1}: {sorted(members)}")
```

## Advanced Tutorials

### 9. Custom Workflow Development

Learn to build custom analysis pipelines using METAINFORMANT components.

#### Workflow Orchestration

Workflow steps are plain functions referenced as `module:function`. A step receives its `params` as keyword arguments, plus the result of every step it `depends_on` as a keyword argument named after that upstream step. Put the step functions in a module that Python can import (for example `my_pipeline.py` next to the script you run, or anywhere on `PYTHONPATH`):

```python-snippet
# File: my_pipeline.py
from metainformant.core import io
from metainformant.dna.sequence import core
from metainformant.information.workflow import workflows


def load_sequences(path: str) -> dict[str, str]:
    """Step 1: read a FASTA file."""
    return core.read_fasta(path)


def entropy_report(load_sequences: dict[str, str], k: int = 2) -> dict:
    """Step 2: receives the result of the ``load_sequences`` step as ``load_sequences``."""
    return workflows.batch_entropy_analysis(list(load_sequences.values()), k=k)


def save_report(entropy_report: dict, path: str) -> str:
    """Step 3: write the report to JSON."""
    io.dump_json(entropy_report, path)
    return path
```

Then describe the steps in a configuration and run it. `run_config_based_workflow` takes the path of a configuration file:

```python-snippet
from metainformant.core import io
from metainformant.core.execution.workflow import run_config_based_workflow

# Define custom workflow configuration: "steps" maps step names to a function reference,
# its parameters and the steps it depends on
workflow_config = {
    "steps": {
        "load_sequences": {
            "function": "my_pipeline:load_sequences",
            "params": {"path": "data/sequences.fasta"},
        },
        "entropy_report": {
            "function": "my_pipeline:entropy_report",
            "params": {"k": 2},
            "depends_on": ["load_sequences"],
        },
        "save_report": {
            "function": "my_pipeline:save_report",
            "params": {"path": "output/custom_workflow/entropy_report.json"},
            "depends_on": ["entropy_report"],
        },
    }
}

# Execute workflow
io.dump_json(workflow_config, "config/custom_workflow.json")
results = run_config_based_workflow("config/custom_workflow.json")
print(f"Workflow completed: {results['success']}")
print(f"Execution order: {results['execution_order']}")
```

The result dictionary also holds each step's output under `results["results"]`, any failures under `results["errors"]`, and `total_duration`. When a step fails, the steps that depend on it are skipped.

### 10. Performance Optimization Tutorial

Techniques for optimizing METAINFORMANT performance with large datasets.

#### Parallel Processing

```python-snippet
from metainformant.core.execution import parallel
import time

# Example: Parallel entropy calculation
def calculate_entropy(seq):
    """Calculate entropy for a single sequence."""
    from metainformant.information.metrics.analysis.analysis import analyze_sequence_information
    # Restrict the methods to entropy: the default "patterns" method scans repeats and palindromes
    # and is very slow on long sequences
    report = analyze_sequence_information(seq, k_values=[1, 2], methods=["entropy"])
    return report["k_mer_analysis"]["k1"]["entropy"]

def main():
    # Large dataset simulation
    sequences = [f"ATCG{i}" * 100 for i in range(1000)]  # 1000 sequences

    # Sequential processing
    start_time = time.time()
    sequential_results = [calculate_entropy(seq) for seq in sequences]
    sequential_time = time.time() - start_time

    # Parallel processing: process_map for CPU-bound work (the function must be defined at module
    # level); use parallel.thread_map for I/O-bound work such as downloads
    start_time = time.time()
    parallel_results = parallel.process_map(calculate_entropy, sequences, max_workers=4)
    parallel_time = time.time() - start_time

    print("Performance Results:")
    print(f"Sequential time: {sequential_time:.2f}s")
    print(f"Parallel time: {parallel_time:.2f}s")
    print(f"Speedup: {sequential_time/parallel_time:.1f}x")
    print(f"Results identical: {sequential_results == parallel_results}")

if __name__ == "__main__":  # required where worker processes are spawned (macOS, Windows)
    main()
```

#### Caching for Repeated Computations

```python-snippet
from functools import lru_cache
from metainformant.math.population_genetics import coalescent
import time

# tajima_constants is not cached internally; memoize it when sample sizes repeat
cached_tajima_constants = lru_cache(maxsize=None)(coalescent.tajima_constants)

# Demonstrate caching benefits
sample_sizes = [10, 50, 100, 10, 50, 100]  # Repeated calculations

print("Tajima constants calculation (with caching):")
start_time = time.time()

for n in sample_sizes:
    constants = cached_tajima_constants(n)
    print(f"n={n}: a1={constants['a1']:.2f}")

total_time = time.time() - start_time
print(f"Total time: {total_time:.3f}s")
print(f"(Note: {cached_tajima_constants.cache_info().hits} of {len(sample_sizes)} calls were served from the cache)")
```

## Best Practices

### 1. Data Management
- Always validate input data before analysis
- Use appropriate file formats for your data type
- Implement proper error handling for missing files

### 2. Performance Optimization
- Use parallel processing for independent operations
- Leverage caching for expensive repeated calculations
- Vectorize operations with NumPy when possible

### 3. Result Interpretation
- Understand the biological meaning of statistical results
- Validate results with known biological expectations
- Use appropriate significance thresholds

### 4. Reproducibility
- Save workflow configurations
- Document analysis parameters
- Use version control for analysis scripts

### 5. Error Handling
- Implement proper exception handling
- Provide informative error messages
- Log analysis progress and issues

## Getting Help

### Common Issues and Solutions

**Import Errors**: Ensure METAINFORMANT is properly installed with all dependencies:
```bash
uv pip install metainformant[scientific,ml,networks]
```

**Memory Issues**: For large datasets, use chunked processing or increase system memory.

**Slow Performance**: Enable parallel processing and caching where available.

**File Format Issues**: Check documentation for supported formats and conversion tools.

### Documentation Resources
- **[Documentation Index](index.md)**: Entry point for all user documentation
- **[Documentation Guide](DOCUMENTATION_GUIDE.md)**: Navigation and conventions
- **[FAQ](FAQ.md)**: Common issues and solutions
- **[Examples](../examples/README.md)**: Example code and usage patterns

---

These tutorials provide a foundation for using METAINFORMANT effectively across all biological domains. Start with the basic tutorials and progress to advanced topics as needed.
