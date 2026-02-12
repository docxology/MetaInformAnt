# Evaluation

Model validation and cross-validation utilities for biological data analysis.

## Contents

| File | Purpose |
|------|---------|
| `validation.py` | Stratified train/test splits, k-fold CV, permutation importance, bootstrap validation |

## Key Functions

| Function | Description |
|----------|-------------|
| `train_test_split_biological` | Stratified splitting for imbalanced biological classes |
| `cross_validate_biological` | K-fold cross-validation with biological-aware stratification |
| `permutation_importance_biological` | Feature importance via permutation testing |
| `bootstrap_evaluate` | Bootstrap confidence intervals for model metrics |

## Usage

```python
from metainformant.ml.evaluation.validation import (
    train_test_split_biological,
    cross_validate_biological,
)

X_train, X_test, y_train, y_test = train_test_split_biological(
    X, y, test_size=0.2, stratify=y, random_state=42
)
cv_results = cross_validate_biological(model, X, y, cv=5, scoring="f1")
```
