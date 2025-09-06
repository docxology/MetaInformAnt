"""Machine learning tools for bioinformatics applications.

This module provides specialized ML algorithms and tools for biological data:
- Feature selection for high-dimensional biological data
- Classification and regression for biological datasets
- Dimensionality reduction techniques
- Cross-validation and model evaluation
- Ensemble methods for biological predictions
"""

from .classification import (
    BiologicalClassifier,
    cross_validate_biological,
    evaluate_classifier,
    train_ensemble_classifier,
)
from .dimensionality import biological_embedding, reduce_dimensions_pca, reduce_dimensions_tsne, reduce_dimensions_umap
from .features import (
    biological_feature_ranking,
    select_features_recursive,
    select_features_stability,
    select_features_univariate,
)
from .regression import BiologicalRegressor, evaluate_regressor, train_ensemble_regressor
from .validation import bootstrap_validate, cross_validate, k_fold_split, learning_curve, train_test_split

__all__ = [
    # Feature selection
    "select_features_univariate",
    "select_features_recursive",
    "select_features_stability",
    "biological_feature_ranking",
    # Classification
    "BiologicalClassifier",
    "train_ensemble_classifier",
    "evaluate_classifier",
    "cross_validate_biological",
    # Regression
    "BiologicalRegressor",
    "train_ensemble_regressor",
    "evaluate_regressor",
    # Dimensionality reduction
    "reduce_dimensions_pca",
    "reduce_dimensions_umap",
    "reduce_dimensions_tsne",
    "biological_embedding",
    # Validation
    "train_test_split",
    "k_fold_split",
    "cross_validate",
    "bootstrap_validate",
    "learning_curve",
]
