# Models

Classification and regression model wrappers for biological data analysis.

## Contents

| File | Purpose |
|------|---------|
| `classification.py` | `BiologicalClassifier` — ensemble classification (RF, SVM, GBM, voting) with evaluation |
| `regression.py` | `BiologicalRegressor` — regression models (linear, RF, SVR, elastic net) for trait prediction |

## Key Classes

| Class | Description |
|-------|-------------|
| `BiologicalClassifier` | Sklearn wrapper supporting algorithm-based or model-based initialization, CV scoring |
| `BiologicalRegressor` | Sklearn wrapper for quantitative trait prediction with R2, MAE, RMSE evaluation |

## Key Functions

| Function | Description |
|----------|-------------|
| `classification.cross_validate_biological()` | Cross-validate a classifier chosen by name (`'rf'`, `'gb'`, `'lr'`, `'ensemble'`); `random_state` seeds the folds and the model, and `parallel_jobs` governs the cross-validation and the ensemble unless an explicit `n_jobs` keyword is given (`parallel_jobs=1` starts no worker processes) |
| `classification.compare_classifiers()` | Run `cross_validate_biological` per method; a failed method is reported as `{"error": message}` and never selected |
| `regression.analyze_prediction_uncertainty()` | Bootstrap prediction intervals; `random_state` seeds the bootstrap samples (drawn from a local generator) and the regressors, and the global NumPy random state is untouched |

## Usage

```python
from metainformant.ml.models.classification import BiologicalClassifier
from metainformant.ml.models.regression import BiologicalRegressor

clf = BiologicalClassifier(algorithm="random_forest", random_state=42)
clf.fit(X_train, y_train)
report = clf.evaluate(X_test, y_test)

reg = BiologicalRegressor(algorithm="linear", random_state=42)
reg.fit(X_train, y_train)
metrics = reg.evaluate(X_test, y_test)
```
