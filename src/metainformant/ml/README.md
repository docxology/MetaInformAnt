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
from metainformant.ml import classification

# Train classifier
model = classification.train_classifier(features, labels, method="random_forest")

# Make predictions
predictions = classification.predict(model, test_features)

# Evaluate performance
metrics = classification.evaluate_classification(predictions, test_labels)
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
from metainformant.ml import regression

# Train regression model
model = regression.train_regressor(features, targets, method="xgboost")

# Make predictions
predictions = regression.predict(model, test_features)

# Evaluate performance
metrics = regression.evaluate_regression(predictions, test_targets)
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
from metainformant.ml import features

# Select important features
selected = features.select_features(features, labels, k=100)
importance_scores = features.feature_importance(model, features)

# Biological interpretation
biological_features = features.interpret_features(selected_features, annotation)
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
from metainformant.ml import validation

# Cross-validation
cv_scores = validation.cross_validate(model, features, labels, cv=5)

# Bootstrap confidence intervals
bootstrap_results = validation.bootstrap_validation(model, features, labels)

# Model comparison
best_model = validation.compare_models([model1, model2, model3], features, labels)
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
from metainformant.ml import dimensionality

# Reduce dimensions
pca_result = dimensionality.pca(features, n_components=50)
umap_result = dimensionality.umap(features, n_neighbors=15)
```

## Biological Applications

### Sequence-Based Prediction
```python
from metainformant.dna import sequences
from metainformant.ml import classification

# Classify sequences by function
sequences = sequences.read_fasta("sequences.fasta")
features = sequences.extract_features(sequences)
labels = sequences.get_functional_labels(sequences)

model = classification.train_classifier(features, labels)
```

### Expression-Based Phenotype Prediction
```python
from metainformant.rna import workflow
from metainformant.ml import regression

# Predict continuous phenotypes
expression = workflow.extract_expression_patterns(rna_data)
phenotypes = get_phenotype_data()

model = regression.train_regressor(expression, phenotypes)
predictions = regression.predict(model, test_expression)
```

### Network-Based Learning
```python
from metainformant.networks import ppi
from metainformant.ml import classification

# Use network features for prediction
network_features = ppi.extract_network_features(interaction_network)
model = classification.train_classifier(network_features, labels)
```

## Performance Features

- **Scalability**: Efficient algorithms for large biological datasets
- **Parallel Processing**: Multi-core support for computationally intensive tasks
- **Memory Optimization**: Streaming processing for large feature matrices
- **GPU Support**: CUDA acceleration where applicable

## Model Interpretability

### Feature Importance
```python
from metainformant.ml import features

# Analyze feature contributions
importance = features.permutation_importance(model, features, labels)
biological_interpretation = features.interpret_importance(importance, annotation)
```

### Model Explanation
```python
from metainformant.ml import explain

# Explain individual predictions
explanation = explain.explain_prediction(model, instance, features)
shap_values = explain.shap_analysis(model, features)
```

## Integration with Other Modules

### With Visualization Module
```python
from metainformant.ml import features
from metainformant.visualization import heatmap

# Visualize feature importance
importance = features.feature_importance(model, features)
heatmap(importance.reshape(1, -1), title="Feature Importance")
```

### With Statistical Analysis
```python
from metainformant.ml import validation
from metainformant.math import statistics

# Statistical validation of ML results
cv_results = validation.cross_validate(model, features, labels)
p_values = statistics.permutation_test(cv_results)
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
