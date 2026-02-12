# Model Evaluation

The evaluation module provides comprehensive model validation tools designed for biological data, including stratified train/test splitting, k-fold cross-validation, bootstrap validation, and learning curve analysis.

## Key Concepts

### Stratified Splitting

Biological classification tasks often involve imbalanced classes (rare disease subtypes, uncommon cell types). Stratified splitting ensures each fold or split maintains the original class distribution. The module automatically applies stratification when fewer than 20 unique target values are detected.

### Cross-Validation Strategies

- **Stratified K-Fold**: Preserves class proportions in each fold. Default for classification tasks.
- **Standard K-Fold**: For regression tasks where stratification is not applicable.
- **Leave-One-Out**: Maximum data utilization for small sample sizes (common in biological studies).

### Bootstrap Validation

Resampling with replacement to estimate model performance confidence intervals. Particularly useful for small biological datasets where cross-validation may have high variance.

### Learning Curves

Plot model performance as a function of training set size to diagnose overfitting, underfitting, and determine whether more samples would improve performance.

## Function Reference

### train_test_split_biological

```python
def train_test_split_biological(
    X: np.ndarray,
    y: np.ndarray,
    test_size: float = 0.2,
    stratify: np.ndarray | None = None,
    random_state: int | None = None,
    **kwargs: Any,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
```

Split biological data into train/test sets with automatic stratification for classification targets. Returns (X_train, X_test, y_train, y_test).

### cross_validate_biological

```python
def cross_validate_biological(
    model: Any,
    X: np.ndarray,
    y: np.ndarray,
    cv: int = 5,
    scoring: str | list[str] | None = None,
    stratified: bool = True,
    random_state: int | None = None,
) -> dict[str, Any]
```

Perform cross-validation with biological data considerations. Supports multiple scoring metrics simultaneously. Returns fold scores, mean, std, and per-fold details.

### bootstrap_validation

```python
def bootstrap_validation(
    model: Any,
    X: np.ndarray,
    y: np.ndarray,
    n_iterations: int = 100,
    sample_fraction: float = 0.8,
    random_state: int | None = None,
) -> dict[str, Any]
```

Bootstrap validation with confidence intervals. Returns metric distributions, mean, std, and 95% confidence intervals.

### learning_curve_analysis

```python
def learning_curve_analysis(
    model: Any,
    X: np.ndarray,
    y: np.ndarray,
    train_sizes: list[float] | None = None,
    cv: int = 5,
    scoring: str = "accuracy",
) -> dict[str, Any]
```

Compute learning curves. Default train sizes are [0.1, 0.2, 0.4, 0.6, 0.8, 1.0]. Returns training and validation scores at each size for diagnosing bias-variance tradeoffs.

### permutation_test

```python
def permutation_test(
    model: Any,
    X: np.ndarray,
    y: np.ndarray,
    n_permutations: int = 100,
    scoring: str = "accuracy",
) -> dict[str, Any]
```

Permutation test to assess whether model performance is significantly better than chance. Returns observed score, null distribution, and p-value.

## Usage Examples

```python
from metainformant.ml.evaluation.validation import (
    train_test_split_biological,
    cross_validate_biological,
    bootstrap_validation,
    learning_curve_analysis,
)
from metainformant.ml.models.classification import BiologicalClassifier

# Stratified train/test split
X_train, X_test, y_train, y_test = train_test_split_biological(
    X, y, test_size=0.2, random_state=42
)

# Cross-validation
model = BiologicalClassifier(algorithm="random_forest", random_state=42)
cv_results = cross_validate_biological(
    model.model, X, y, cv=5,
    scoring=["accuracy", "f1_weighted"],
)
print(f"CV Accuracy: {cv_results['mean_accuracy']:.3f}")

# Bootstrap validation for confidence intervals
boot = bootstrap_validation(model.model, X, y, n_iterations=200)
print(f"Accuracy: {boot['mean']:.3f} (95% CI: {boot['ci_lower']:.3f}-{boot['ci_upper']:.3f})")

# Learning curve to check for overfitting
lc = learning_curve_analysis(model.model, X, y, cv=5, scoring="accuracy")
# lc['train_sizes'], lc['train_scores'], lc['validation_scores']

# Permutation test for statistical significance
perm = permutation_test(model.model, X, y, n_permutations=100)
print(f"p-value: {perm['p_value']:.4f}")
```

## Configuration

- **Environment prefix**: `ML_`
- **Required**: numpy, scikit-learn
- Automatic stratification threshold: 20 unique target values
- Bootstrap uses out-of-bag samples for evaluation by default
- Learning curves use 5-fold cross-validation at each training size

## Related Modules

- `ml.models` -- Classification and regression models to evaluate
- `ml.features` -- Feature selection (evaluate with cross-validation)
- `ml.interpretability` -- Model explanation methods
- `ml.automl` -- Automated model selection with built-in evaluation
