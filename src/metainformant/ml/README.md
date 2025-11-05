# Machine Learning Module

The `ml` module provides statistical and machine learning methods for biological data analysis, including classification, regression, feature selection, and model validation.

## Overview

This module offers a comprehensive toolkit for applying machine learning techniques to biological datasets, with emphasis on biological interpretation and robust validation.

## Submodules

### Classification (`classification.py`)
Supervised learning methods for biological prediction tasks.

**Key Features:**
- Binary and multi-class classification
- Biological sequence classification
- Expression-based phenotype prediction
- Model interpretability tools

**Usage:**
```python
from metainformant.ml import BiologicalClassifier, evaluate_classifier

# Create and train classifier
classifier = BiologicalClassifier(algorithm="random_forest", random_state=42)
classifier.fit(X_train, y_train)

# Make predictions
predictions = classifier.predict(X_test)

# Get feature importance
importance = classifier.get_feature_importance()
top_features = np.argsort(importance)[-10:]

# Evaluate performance
metrics = evaluate_classifier(classifier, X_test, y_test)
```

### Regression (`regression.py`)
Continuous trait prediction and modeling.

**Key Features:**
- Linear and non-linear regression
- Regularization methods (Lasso, Ridge, Elastic Net)
- Survival analysis
- Time-series prediction

**Usage:**
```python
from metainformant.ml import BiologicalRegressor, evaluate_regressor

# Create and train regressor
regressor = BiologicalRegressor(algorithm="linear", random_state=42)
regressor.fit(X_train, y_train)

# Make predictions
predictions = regressor.predict(X_test)

# Evaluate performance
metrics = evaluate_regressor(regressor, X_test, y_test)
```

### Feature Selection (`features.py`)
Dimensionality reduction and feature importance analysis.

**Key Features:**
- Univariate feature selection
- Recursive feature elimination
- L1-based selection
- Biological feature ranking

**Usage:**
```python
from metainformant.ml import (
    select_features_univariate,
    select_features_recursive,
    select_features_stability,
    biological_feature_ranking
)

# Select features using univariate tests
selected_X, selected_indices = select_features_univariate(X, y, k=100, method="f_score")

# Recursive feature elimination
selected_X_rfe, selected_indices_rfe = select_features_recursive(X, y, k=50)

# Feature ranking
rankings = biological_feature_ranking(X, y)
```

### Model Validation (`validation.py`)
Comprehensive model assessment and validation.

**Key Features:**
- Cross-validation strategies
- Bootstrap resampling
- Permutation testing
- Model comparison and selection

**Usage:**
```python
from metainformant.ml import (
    cross_validate,
    cross_validate_biological,
    bootstrap_validate,
    train_test_split,
    k_fold_split
)

# Cross-validation for biological data
cv_scores = cross_validate_biological(classifier, X, y, cv=5)

# Bootstrap validation
bootstrap_results = bootstrap_validate(classifier, X, y, n_bootstrap=100)

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
```

### Dimensionality Reduction (`dimensionality.py`)
Manifold learning and dimensionality reduction techniques.

**Key Features:**
- Principal Component Analysis (PCA)
- t-SNE and UMAP
- Non-negative Matrix Factorization (NMF)
- Independent Component Analysis (ICA)

**Usage:**
```python
from metainformant.ml import (
    reduce_dimensions_pca,
    reduce_dimensions_umap,
    reduce_dimensions_tsne,
    biological_embedding
)

# PCA dimensionality reduction
X_reduced, components, variance = reduce_dimensions_pca(X, n_components=50)

# UMAP embedding
X_umap = reduce_dimensions_umap(X, n_neighbors=15, n_components=2)

# Biological embedding
embedding = biological_embedding(X, method="umap")
```

## Biological Applications

### Sequence-Based Prediction
```python
from metainformant.dna import sequences
from metainformant.ml import BiologicalClassifier

# Classify sequences by function
sequences = sequences.read_fasta("sequences.fasta")
# Extract features and labels as needed
classifier = BiologicalClassifier(algorithm="random_forest")
classifier.fit(X_train, y_train)
```

### Expression-Based Phenotype Prediction
```python
from metainformant.rna import workflow
from metainformant.ml import BiologicalRegressor

# Predict continuous phenotypes
# Get expression and phenotype data
regressor = BiologicalRegressor(algorithm="ridge")
regressor.fit(expression_train, phenotypes_train)
predictions = regressor.predict(expression_test)
```

### Network-Based Learning
```python
from metainformant.networks import ppi
from metainformant.ml import BiologicalClassifier

# Use network features for prediction
# Extract network features and train classifier
classifier = BiologicalClassifier(algorithm="random_forest")
classifier.fit(network_features_train, labels_train)
```

## Performance Features

- **Scalability**: Efficient algorithms for large biological datasets
- **Parallel Processing**: Multi-core support for computationally intensive tasks
- **Memory Optimization**: Streaming processing for large feature matrices
- **GPU Support**: CUDA acceleration where applicable

## Model Interpretability

### Feature Importance
```python
from metainformant.ml import biological_feature_ranking

# Analyze feature contributions
rankings = biological_feature_ranking(X, y)
# Access feature importance from classifier
importance = classifier.feature_importance_
```

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import sequences, composition
from metainformant.ml import BiologicalClassifier, kmer_frequencies
import numpy as np

# Sequence classification using k-mer features
dna_seqs = sequences.read_fasta("sequences.fasta")

# Extract k-mer features from sequences
def extract_kmer_features(seq_dict, k=3):
    features = []
    for seq_id, seq in seq_dict.items():
        kmers = composition.kmer_frequencies(seq, k=k)
        # Convert to feature vector
        feature_vec = [kmers.get(kmer, 0) for kmer in sorted(kmers.keys())]
        features.append(feature_vec)
    return np.array(features)

X = extract_kmer_features(dna_seqs)
# Create labels (e.g., functional categories)
y = np.array([0, 1, 0, 1, 0])  # Binary classification labels

# Train classifier
classifier = BiologicalClassifier(algorithm="random_forest", random_state=42)
classifier.fit(X, y)
predictions = classifier.predict(X)
```

### With RNA Module
```python
from metainformant.rna import workflow
from metainformant.ml import BiologicalRegressor, select_features_univariate

# Expression-based phenotype prediction
# Load expression data from RNA workflow
expression_data = workflow.extract_expression("expression.tsv")
phenotypes = np.array([...])  # Continuous trait values

# Feature selection for expression-based prediction
X_selected, selected_indices = select_features_univariate(
    expression_data.values,
    phenotypes,
    k=100,
    method="f_score"
)

# Train regression model
regressor = BiologicalRegressor(algorithm="ridge", random_state=42)
regressor.fit(X_selected, phenotypes)
predictions = regressor.predict(X_selected)
```

### With GWAS Module
```python
from metainformant.gwas import parse_vcf_full, apply_qc_filters
from metainformant.ml import BiologicalClassifier, cross_validate_biological

# Variant effect prediction
vcf_data = parse_vcf_full("variants.vcf")
qc_variants = apply_qc_filters(vcf_data, maf_threshold=0.05)

# Extract variant features (e.g., allele frequencies, functional annotations)
variant_features = extract_variant_features(qc_variants)
# Create labels (e.g., pathogenic vs benign)
variant_labels = np.array([0, 1, 0, 1, 1])  # Binary labels

# Train variant effect classifier
classifier = BiologicalClassifier(algorithm="gradient_boosting", random_state=42)
classifier.fit(variant_features, variant_labels)

# Cross-validate model
cv_results = cross_validate_biological(classifier, variant_features, variant_labels, cv=5)
```

### With Single-Cell Module
```python
from metainformant.singlecell import load_count_matrix, compute_pca, leiden_clustering
from metainformant.ml import BiologicalClassifier, reduce_dimensions_umap

# Cell type classification from single-cell data
sc_data = load_count_matrix("singlecell_counts.h5ad")
sc_data = compute_pca(sc_data, n_components=50)

# Use clustering results as training labels
sc_data = leiden_clustering(sc_data, resolution=1.0)
cell_types = sc_data.obs["cluster"].values  # Cluster assignments

# Extract features from PCA or raw counts
X_features = sc_data.obsm["X_pca"]  # PCA features
# Or use raw expression
# X_features = sc_data.X

# Train cell type classifier
classifier = BiologicalClassifier(algorithm="random_forest", random_state=42)
classifier.fit(X_features, cell_types)

# Visualize with UMAP
X_umap = reduce_dimensions_umap(X_features, n_neighbors=15, n_components=2)
```

### With Multiomics Module
```python
from metainformant.multiomics import MultiOmicsData, joint_pca
from metainformant.ml import BiologicalClassifier, select_features_recursive

# Multi-omic feature integration for classification
omics_data = MultiOmicsData(
    genomics=genomics_df,
    transcriptomics=transcriptomics_df,
    proteomics=proteomics_df
)

# Joint PCA for dimensionality reduction
embeddings, loadings, variance = joint_pca(omics_data, n_components=50)

# Use joint embeddings as features
X = embeddings
y = np.array([...])  # Class labels

# Recursive feature elimination for multi-omic data
X_selected, selected_indices = select_features_recursive(X, y, k=30)

# Train multi-omic classifier
classifier = BiologicalClassifier(algorithm="support_vector", random_state=42)
classifier.fit(X_selected, y)
predictions = classifier.predict(X_selected)
```

### With Networks Module
```python
from metainformant.networks import create_network, centrality_measures, detect_communities
from metainformant.ml import BiologicalClassifier

# Network-based learning using network topology features
ppi_network = load_string_interactions("interactions.tsv")
network = ppi_network.create_network()

# Extract network features (centrality measures, community membership)
centralities = centrality_measures(network)
communities = detect_communities(network, method="louvain")

# Create feature matrix from network properties
def extract_network_features(network, protein_ids):
    features = []
    for pid in protein_ids:
        node_features = [
            centralities["degree"].get(pid, 0),
            centralities["betweenness"].get(pid, 0),
            centralities["closeness"].get(pid, 0),
            communities.get(pid, -1)
        ]
        features.append(node_features)
    return np.array(features)

network_features = extract_network_features(network, protein_ids)
labels = np.array([...])  # Functional labels

# Train network-based classifier
classifier = BiologicalClassifier(algorithm="random_forest", random_state=42)
classifier.fit(network_features, labels)
```

### With Visualization Module
```python
from metainformant.ml import biological_feature_ranking, reduce_dimensions_umap
from metainformant.visualization import scatter_plot, heatmap

# Visualize feature importance
importance = biological_feature_ranking(X, y)
# Create sorted importance plot
sorted_importance = sorted(importance.items(), key=lambda x: x[1], reverse=True)[:20]
features, scores = zip(*sorted_importance)
ax = scatter_plot(range(len(features)), scores, 
                  xlabel="Feature Rank", ylabel="Importance Score",
                  title="Top 20 Feature Importance")

# Visualize dimensionality reduction
X_umap = reduce_dimensions_umap(X, n_components=2)
ax = scatter_plot(X_umap[:, 0], X_umap[:, 1],
                  xlabel="UMAP 1", ylabel="UMAP 2",
                  title="UMAP Embedding of Feature Space")

# Visualize feature importance heatmap
importance_matrix = np.array([importance[f] for f in feature_names]).reshape(10, 10)
ax = heatmap(importance_matrix, title="Feature Importance Heatmap")
```

### With Statistical Analysis
```python
from metainformant.ml import cross_validate_biological, bootstrap_validate

# Statistical validation of ML results
cv_results = cross_validate_biological(classifier, X, y, cv=5)
print(f"Mean CV accuracy: {np.mean(cv_results['test_score']):.3f}")
print(f"Std CV accuracy: {np.std(cv_results['test_score']):.3f}")

# Bootstrap validation for robust performance estimation
bootstrap_results = bootstrap_validate(classifier, X, y, n_bootstrap=100)
print(f"Bootstrap mean accuracy: {bootstrap_results['mean_score']:.3f}")
print(f"Bootstrap 95% CI: {bootstrap_results['ci_95']}")
```

## Testing and Validation

- **Algorithm Correctness**: Verification against known implementations
- **Biological Relevance**: Testing with real biological datasets
- **Performance Benchmarking**: Scalability testing
- **Robustness**: Edge case and error handling validation

## Dependencies

- **Core**: scikit-learn for ML algorithms
- **Optional**: xgboost, lightgbm for gradient boosting
- **Deep Learning**: tensorflow/pytorch for neural networks (optional)

This module provides a complete machine learning toolkit tailored for biological data analysis and interpretation.
