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
