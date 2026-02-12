# Feature Selection and Dimensionality Reduction

The features module provides feature selection and dimensionality reduction methods designed for biological data analysis. It includes univariate statistical tests, model-based selection, and manifold learning methods (PCA, t-SNE, UMAP, ICA).

## Key Concepts

### Feature Selection

Feature selection identifies the most informative features (genes, biomarkers, etc.) from high-dimensional biological data. Three main approaches are provided:

- **Univariate selection**: Ranks features independently using statistical tests (ANOVA F-test, chi-squared, mutual information)
- **Model-based selection**: Uses feature importance from trained models (Random Forest, Lasso) to select features
- **Recursive elimination**: Iteratively removes the least important features, retraining the model at each step

### Dimensionality Reduction

Dimensionality reduction projects high-dimensional data into a lower-dimensional space while preserving structure:

- **PCA**: Linear projection preserving maximum variance. Best for visualization of global structure and initial dimensionality reduction before clustering.
- **t-SNE**: Non-linear embedding optimized for preserving local neighborhood structure. Best for visualization of clusters in 2D/3D.
- **UMAP**: Non-linear embedding preserving both local and global structure. Faster than t-SNE and scales better.
- **ICA**: Independent Component Analysis for separating mixed signals into independent sources.

## Function Reference

### Feature Selection

```python
def select_features_univariate(
    X: np.ndarray,
    y: np.ndarray,
    method: str = "f_classif",
    k: int | str = "all",
    **kwargs: Any,
) -> tuple[np.ndarray, np.ndarray]
```

Select features using univariate statistical tests. Methods: `"f_classif"` (ANOVA F-test), `"chi2"` (chi-squared), `"mutual_info_classif"`, `"f_regression"`, `"mutual_info_regression"`. Returns (X_selected, selected_indices).

### Dimensionality Reduction

```python
def pca_reduction(
    X: np.ndarray,
    n_components: int | None = None,
    scale_data: bool = True,
    random_state: int | None = None,
    **kwargs: Any,
) -> tuple[np.ndarray, PCA]
```

PCA dimensionality reduction. Optionally standardizes data first. Returns (transformed_data, fitted_pca_model).

```python
def tsne_reduction(
    X: np.ndarray,
    n_components: int = 2,
    perplexity: float = 30.0,
    random_state: int | None = None,
    **kwargs: Any,
) -> tuple[np.ndarray, Any]
```

t-SNE embedding for visualization. Returns (embedded_data, tsne_model).

```python
def umap_reduction(
    X: np.ndarray,
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    random_state: int | None = None,
    **kwargs: Any,
) -> tuple[np.ndarray, Any]
```

UMAP embedding. Requires `umap-learn` package. Returns (embedded_data, umap_model).

```python
def ica_reduction(
    X: np.ndarray,
    n_components: int | None = None,
    random_state: int | None = None,
    **kwargs: Any,
) -> tuple[np.ndarray, Any]
```

ICA for independent component extraction. Returns (transformed_data, ica_model).

## Usage Examples

```python
from metainformant.ml.features import features, dimensionality
import numpy as np

# Univariate feature selection
X = np.random.randn(100, 500)  # 100 samples, 500 features
y = np.random.randint(0, 3, 100)  # 3-class labels
X_selected, indices = features.select_features_univariate(X, y, method="f_classif", k=50)
print(f"Selected {X_selected.shape[1]} features")

# PCA reduction
X_pca, pca_model = dimensionality.pca_reduction(X, n_components=20, scale_data=True)
print(f"Explained variance: {sum(pca_model.explained_variance_ratio_):.2%}")

# t-SNE for visualization
X_tsne, _ = dimensionality.tsne_reduction(X_pca, n_components=2, perplexity=30)

# UMAP embedding (requires umap-learn)
X_umap, _ = dimensionality.umap_reduction(X_pca, n_components=2, n_neighbors=15)

# ICA
X_ica, ica_model = dimensionality.ica_reduction(X, n_components=10)
```

## Configuration

- **Environment prefix**: `ML_`
- **Required**: numpy, scikit-learn
- **Optional**: umap-learn (for UMAP)
- All methods accept a `random_state` parameter for reproducibility
- PCA scaling is enabled by default (recommended for gene expression data)

## Related Modules

- `ml.models` -- Classification and regression using selected features
- `ml.evaluation` -- Cross-validation to evaluate feature selection quality
- `ml.interpretability` -- Feature importance methods (SHAP, permutation importance)
- `ml.automl` -- Automated feature preprocessing pipeline
