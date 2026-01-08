"""Feature selection and engineering module.

This module provides tools for dimensionality reduction and feature selection
specifically tailored for biological data.
"""

from .dimensionality import (
    biological_embedding,
    reduce_dimensions_pca,
    reduce_dimensions_tsne,
    reduce_dimensions_umap,
)
from .features import (
    biological_feature_ranking,
    select_features_recursive,
    select_features_stability,
    select_features_univariate,
)

__all__ = [
    "biological_embedding",
    "reduce_dimensions_pca",
    "reduce_dimensions_tsne",
    "reduce_dimensions_umap",
    "biological_feature_ranking",
    "select_features_recursive",
    "select_features_stability",
    "select_features_univariate",
]
