# Classification and Regression Models

The models module provides wrapper classes for classification and regression tasks optimized for biological data analysis. Both `BiologicalClassifier` and `BiologicalRegressor` provide a unified interface for training, prediction, evaluation, and cross-validation.

## Key Concepts

### BiologicalClassifier

A wrapper for classification models with built-in evaluation metrics. Supports algorithm-based initialization (specify an algorithm name) or model-based initialization (pass a pre-built scikit-learn model).

Supported algorithms:
- `"random_forest"`: RandomForestClassifier
- `"gradient_boosting"`: GradientBoostingClassifier
- `"logistic_regression"`: LogisticRegression (max_iter=1000)
- `"svm"`: SVC with probability=True

### BiologicalRegressor

A wrapper for regression models with evaluation metrics tailored for quantitative biological predictions.

Supported algorithms:
- `"linear"`: LinearRegression
- `"random_forest"`: RandomForestRegressor
- `"gradient_boosting"`: GradientBoostingRegressor
- `"ridge"`: Ridge regression
- `"lasso"`: Lasso regression
- `"svm"`: Support Vector Regression (SVR)

### Evaluation Metrics

**Classification**: accuracy, precision, recall, F1 score, ROC AUC, confusion matrix, classification report.

**Regression**: R-squared, mean squared error, mean absolute error, median absolute error, explained variance score.

## Class Reference

### BiologicalClassifier

```python
class BiologicalClassifier:
    def __init__(
        self,
        model: Any = None,           # Pre-built sklearn model
        model_type: str = "unknown",  # Model type name
        algorithm: str | None = None, # Algorithm name
        random_state: int | None = None,
        **params: Any,
    )

    def fit(self, X, y, feature_names=None) -> "BiologicalClassifier"
    def predict(self, X) -> np.ndarray
    def predict_proba(self, X) -> np.ndarray
    def evaluate(self, X, y) -> dict[str, Any]
    def cross_validate(self, X, y, cv=5) -> dict[str, Any]
    def feature_importances(self) -> dict[str, float] | None
```

### BiologicalRegressor

```python
class BiologicalRegressor:
    def __init__(
        self,
        model: Any = None,
        model_type: str = "unknown",
        algorithm: str | None = None,
        random_state: int | None = None,
        **params: Any,
    )

    def fit(self, X, y, feature_names=None) -> "BiologicalRegressor"
    def predict(self, X) -> np.ndarray
    def evaluate(self, X, y) -> dict[str, Any]
    def cross_validate(self, X, y, cv=5) -> dict[str, Any]
    def feature_importances(self) -> dict[str, float] | None
```

## Usage Examples

```python
from metainformant.ml.models.classification import BiologicalClassifier
from metainformant.ml.models.regression import BiologicalRegressor
import numpy as np

# Classification with Random Forest
clf = BiologicalClassifier(algorithm="random_forest", random_state=42, n_estimators=100)
clf.fit(X_train, y_train, feature_names=gene_names)

predictions = clf.predict(X_test)
probabilities = clf.predict_proba(X_test)

# Evaluate on test set
metrics = clf.evaluate(X_test, y_test)
print(f"Accuracy: {metrics['accuracy']:.3f}, AUC: {metrics['roc_auc']:.3f}")

# Cross-validation
cv_results = clf.cross_validate(X, y, cv=5)
print(f"CV accuracy: {cv_results['mean_accuracy']:.3f} +/- {cv_results['std_accuracy']:.3f}")

# Feature importances
importances = clf.feature_importances()
top_features = sorted(importances.items(), key=lambda x: x[1], reverse=True)[:10]

# Regression with Gradient Boosting
reg = BiologicalRegressor(algorithm="gradient_boosting", random_state=42)
reg.fit(X_train, y_continuous_train)
reg_metrics = reg.evaluate(X_test, y_continuous_test)
print(f"R2: {reg_metrics['r2']:.3f}, MAE: {reg_metrics['mae']:.3f}")

# Using a pre-built model
from sklearn.ensemble import ExtraTreesClassifier
custom_model = ExtraTreesClassifier(n_estimators=200)
clf = BiologicalClassifier(model=custom_model, model_type="extra_trees")
clf.fit(X_train, y_train)
```

## Configuration

- **Environment prefix**: `ML_`
- **Required**: numpy, scikit-learn
- All models support `random_state` for reproducibility
- Classifiers default to probability-enabled predictions (required for ROC AUC)
- Regressors compute R-squared, MSE, MAE, and explained variance by default

## Related Modules

- `ml.features` -- Feature selection before model training
- `ml.evaluation` -- Extended validation methods (bootstrap, learning curves)
- `ml.interpretability` -- SHAP, LIME, and permutation importance for model explanation
- `ml.automl` -- Automated model selection and hyperparameter optimization
