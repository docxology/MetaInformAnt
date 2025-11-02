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

### With Visualization Module
```python
from metainformant.ml import biological_feature_ranking
from metainformant.visualization import plots

# Visualize feature importance
importance = biological_feature_ranking(X, y)
# Visualize using plots module
```

### With Statistical Analysis
```python
from metainformant.ml import cross_validate_biological

# Statistical validation of ML results
cv_results = cross_validate_biological(classifier, X, y, cv=5)
# Extract scores from cv_results as needed
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
