# Evaluation

Model validation and cross-validation utilities for biological data analysis.

## Contents

| File | Purpose |
|------|---------|
| `validation.py` | Stratified train/test splits, k-fold CV, permutation importance, bootstrap validation, learning curves |

## Key Functions

| Function | Description |
|----------|-------------|
| `train_test_split_biological` | Stratified splitting for imbalanced biological classes |
| `cross_validation_scores` | K-fold cross-validation of an estimator with one or several metrics (stratified for classification-like targets) |
| `cross_validate` | Single-metric cross-validation summary (`scores`, `mean_score`, `std_score`) |
| `permutation_importance_biological` | Feature importance via permutation |
| `bootstrap_validate` | Bootstrap validation of a `(X_train, y_train, X_test, y_test)` model function (`mean_mse`, `std_mse`, per-iteration results) |
| `validate_model_stability` | Accuracy (or R-squared) spread across bootstrap samples |
| `learning_curve` | Train/validation scores at increasing training sizes |

`cross_validate_biological`, which cross-validates a classifier chosen by name, lives in `metainformant.ml.models.classification` together with `compare_classifiers`. Unlike `cross_validation_scores`, which logs a warning and returns an empty array for a metric that could not be computed, `cross_validate_biological` raises: the original exception propagates with the failing metric named in an exception note, and non-finite fold scores raise a `ValueError`, so no metric entry is ever `None` or NaN. `compare_classifiers` reports a failed method as `{"error": message}`, never selects it as `best_method`, and returns `best_method=None` when every method failed.

## Usage

```python
import numpy as np
from sklearn.linear_model import LogisticRegression
from metainformant.ml.evaluation.validation import (
    train_test_split_biological,
    cross_validation_scores,
)

rng = np.random.default_rng(0)
X = rng.normal(size=(100, 6))
y = (X[:, 0] + rng.normal(scale=0.5, size=100) > 0).astype(int)
model = LogisticRegression(max_iter=500)

X_train, X_test, y_train, y_test = train_test_split_biological(
    X, y, test_size=0.2, stratify=y, random_state=42
)
cv_scores = cross_validation_scores(model, X, y, cv=5, scoring=["accuracy", "f1"])
print(cv_scores["f1"].mean())
```
