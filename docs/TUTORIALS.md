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
from metainformant import dna

# Analyze sequence composition
for name, seq in sequences.items():
    gc_content = dna.composition.gc_content(seq)
    length = dna.sequences.sequence_length(seq)
    print(f"{name}: GC={gc_content:.2f}, Length={length}")

# Find motifs
motifs = ["ATCG", "GCTA"]
for name, seq in sequences.items():
    for motif in motifs:
        positions = dna.sequences.find_motifs(seq, [motif])
        print(f"{name} - {motif}: {positions}")
```

#### Advanced Analysis

```python
# Calculate pairwise distances
from metainformant.dna import distances

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

```python
from metainformant.math import coalescent

# Calculate coalescent parameters
summary = coalescent.CoalescentSummary(
    sample_size=len(sequences),
    effective_population_size=10000,
    mutation_rate=0.0001
)

print(f"Population mutation rate (θ): {summary.theta():.2f}")
print(f"Expected total branch length: {summary.total_branch_length():.2f}")
print(f"Expected time to MRCA: {summary.time_to_mrca():.2f}")

# Calculate Tajima's D (would use real sequence data)
# tajima_d = summary.tajimas_D(num_segregating_sites)
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

```python
from metainformant.gwas import association

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
    print(".2e")
```

#### Visualization

```python
from metainformant.gwas import visualization

# Prepare data for Manhattan plot
manhattan_data = [
    {
        "CHROM": "1",
        "POS": i * 1000,  # Position in kb
        "p_value": result["p_value"],
        "ID": result["snp"]
    }
    for i, result in enumerate(results)
]

# Create Manhattan plot
output_path = "output/gwas/manhattan.png"
plot_result = visualization.manhattan_plot(manhattan_data, output_path)
print(f"Manhattan plot saved: {plot_result['status']}")
```

### 4. RNA-seq Analysis Tutorial

Complete RNA-seq analysis pipeline from raw reads to gene expression quantification.

#### Configuration Setup

```python
from metainformant.rna.workflow import AmalgkitWorkflowConfig
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
from metainformant import quality

# Analyze FASTQ quality (assuming you have FASTQ files)
fastq_files = ["data/rna/sample1_R1.fastq.gz", "data/rna/sample1_R2.fastq.gz"]

for fastq_file in fastq_files:
    try:
        quality_report = quality.fastq.assess_quality(fastq_file)
        print(f"FASTQ Quality for {fastq_file}:")
        print(f"  Total reads: {quality_report['total_reads']}")
        print(f"  Mean quality: {quality_report['mean_quality']:.2f}")
        print(f"  GC content: {quality_report['gc_content']:.2f}%")
    except FileNotFoundError:
        print(f"FASTQ file not found: {fastq_file}")
```

#### Expression Quantification

```python
from metainformant.rna.engine.workflow_steps import quant

# Quantify gene expression (would use real BAM files)
# This is a demonstration of the API

quant_config = {
    "bam_files": ["data/rna/aligned/sample1.bam"],
    "gtf_file": "data/rna/annotation.gtf",
    "output_dir": "output/rna_analysis/quant",
    "threads": 4
}

# quant.run_quantification(quant_config)
print("Gene expression quantification would run here")
print("Output: gene expression matrix, quality metrics")
```

### 5. Information Theory Analysis Tutorial

Apply information-theoretic methods to analyze sequence complexity and patterns.

#### Sequence Complexity Analysis

```python
from metainformant.information import analysis

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
        if k in info_result["kmer_analyses"]:
            analysis_k = info_result["kmer_analyses"][k]
            print(f"  k={k}: entropy={analysis_k['entropy']:.3f}, "
                  f"unique_kmers={analysis_k['unique_kmers']}")
```

#### Batch Processing

```python
from metainformant.information import workflows

# Batch entropy analysis
sequences = [
    "ATCGATCGATCG",
    "GCTAGCTAGCTA",
    "ATATATATATAT",
    "CGCGCGCGCGCG"
]

batch_results = workflows.batch_entropy_analysis(sequences, k=2)

print("\nBatch Entropy Analysis Results:")
print(f"Total sequences analyzed: {batch_results['summary']['total_sequences']}")
print(f"Mean entropy: {batch_results['summary']['mean_entropy']:.3f}")
print(f"Entropy range: {batch_results['summary']['min_entropy']:.3f} - "
      f"{batch_results['summary']['max_entropy']:.3f}")

# Detailed results per sequence
for seq_result in batch_results["sequences"]:
    print(f"Sequence {seq_result['index']}: "
          f"entropy={seq_result['entropy']:.3f}, "
          f"length={seq_result['length']}")
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

```python
# Combine multi-omics data
integrated_data = {
    "expression": expression_data,
    "methylation": methylation_data,
    "proteomics": protein_data
}

# Perform joint PCA
from metainformant.multiomics import integration

try:
    joint_pca_result = integration.joint_pca(integrated_data, n_components=3)

    print("Multi-omics Integration Results:")
    print(f"Integrated PCA components shape: {joint_pca_result['expression'].shape}")
    print("Available data types:", list(joint_pca_result.keys()))

except Exception as e:
    print(f"Integration would require real multi-omics data: {e}")
```

### 7. Machine Learning Tutorial

Apply machine learning methods to biological classification and prediction tasks.

#### Classification Example

```python
from metainformant import ml
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

# Train classifier
classifier = ml.classification.train_classifier(X_train, y_train, method="rf")

# Cross-validation
cv_results = ml.classification.cross_validate_classifier(
    classifier, X_train, y_train, cv=5
)

print("Cross-validation Results:")
for metric, value in cv_results.items():
    print(f"  {metric}: {value:.3f}")

# Test predictions
predictions, probabilities = ml.classification.predict_with_confidence(
    classifier, X_test
)

accuracy = np.mean(predictions == y_test)
print(f"Test accuracy: {accuracy:.3f}")
```

#### Feature Selection

```python
# Select important features
selected_features, feature_mask = ml.features.select_features(
    X_train, y_train, method="mutual_info", k=20
)

print(f"Selected {len(selected_features[1])} features out of {n_features}")
print(f"Feature selection method: mutual information")
```

### 8. Network Analysis Tutorial

Analyze biological networks including protein-protein interactions and pathways.

#### Network Construction

```python
from metainformant import networks

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
G = networks.graph.create_network(edges, directed=False)

print("Protein Interaction Network:")
print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")
print(f"Connected components: {networks.graph.connected_components_count(G)}")
```

#### Network Analysis

```python
# Calculate centrality measures
degree_cent = networks.centrality.degree_centrality(G)
betweenness_cent = networks.centrality.betweenness_centrality(G)

print("\nTop nodes by degree centrality:")
for node, centrality in sorted(degree_cent.items(), key=lambda x: x[1], reverse=True)[:3]:
    print(f"  {node}: {centrality:.3f}")

print("\nTop nodes by betweenness centrality:")
for node, centrality in sorted(betweenness_cent.items(), key=lambda x: x[1], reverse=True)[:3]:
    print(f"  {node}: {centrality:.3f}")
```

#### Community Detection

```python
# Detect communities
communities = networks.community.louvain_communities(G)

print(f"\nDetected {len(communities)} communities:")
for i, community in enumerate(communities):
    print(f"  Community {i+1}: {sorted(community)}")
```

## Advanced Tutorials

### 9. Custom Workflow Development

Learn to build custom analysis pipelines using METAINFORMANT components.

#### Workflow Orchestration

```python
from pathlib import Path
from metainformant.core.workflow import run_config_based_workflow

# Define custom workflow configuration
workflow_config = {
    "name": "custom_biology_analysis",
    "steps": [
        {
            "name": "data_ingestion",
            "function": "metainformant.core.io.load_json",
            "params": {"path": "data/input.json"},
            "output_key": "input_data"
        },
        {
            "name": "sequence_analysis",
            "function": "metainformant.dna.sequences.read_fasta",
            "params": {"path": "data/sequences.fasta"},
            "output_key": "sequences"
        },
        {
            "name": "information_analysis",
            "function": "metainformant.information.analysis.batch_entropy_analysis",
            "params": {
                "sequences": "{sequences}",  # Reference previous step output
                "k": 2
            },
            "output_key": "entropy_results"
        },
        {
            "name": "results_save",
            "function": "metainformant.core.io.dump_json",
            "params": {
                "obj": "{entropy_results}",
                "path": "output/analysis_results.json"
            }
        }
    ],
    "output_dir": "output/custom_workflow"
}

# Execute workflow
results = run_config_based_workflow(workflow_config)
print(f"Workflow completed: {results['status']}")
```

### 10. Performance Optimization Tutorial

Techniques for optimizing METAINFORMANT performance with large datasets.

#### Parallel Processing

```python
from metainformant.core import parallel
import time

# Example: Parallel entropy calculation
def calculate_entropy(seq):
    """Calculate entropy for a single sequence."""
    from metainformant.information.analysis import analyze_sequence_information
    return analyze_sequence_information(seq, k_values=[1, 2])["kmer_analyses"][1]["entropy"]

# Large dataset simulation
sequences = [f"ATCG{i}" * 100 for i in range(1000)]  # 1000 sequences

# Sequential processing
start_time = time.time()
sequential_results = [calculate_entropy(seq) for seq in sequences]
sequential_time = time.time() - start_time

# Parallel processing
start_time = time.time()
parallel_results = parallel.run_parallel(calculate_entropy, sequences, max_workers=4)
parallel_time = time.time() - start_time

print("Performance Results:")
print(f"Sequential time: {sequential_time:.2f}s")
print(f"Parallel time: {parallel_time:.2f}s")
print(f"Speedup: {sequential_time/parallel_time:.1f}x")
print(f"Results identical: {sequential_results == parallel_results}")
```

#### Caching for Repeated Computations

```python
from metainformant.math import coalescent
import time

# Demonstrate caching benefits
sample_sizes = [10, 50, 100, 10, 50, 100]  # Repeated calculations

print("Tajima constants calculation (with caching):")
start_time = time.time()

for n in sample_sizes:
    constants = coalescent.tajima_constants(n)
    print(f"n={n}: a1={constants['a1']:.2f}")

total_time = time.time() - start_time
print(f"Total time: {total_time:.3f}s")
print("(Note: First calculations of each n are cached for reuse)")
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
