# AutoML: Automated Model Selection and Optimization

The AutoML module provides automated hyperparameter tuning and model selection for biological machine learning workflows. It includes random search, Bayesian optimization with a Gaussian process surrogate, exhaustive grid search, automatic model selection, and preprocessing pipelines.

## Key Concepts

### Random Search

Randomly samples hyperparameter configurations from specified distributions. More efficient than grid search for high-dimensional parameter spaces, as it explores diverse regions of the search space.

### Bayesian Optimization

Uses a Gaussian process surrogate model to intelligently select the next hyperparameter configuration to evaluate. Balances exploitation (refining known good regions) and exploration (testing unexplored regions) via an acquisition function. More sample-efficient than random search.

### Grid Search

Exhaustive evaluation of all combinations in a discrete hyperparameter grid. Guarantees finding the best configuration within the grid but scales poorly with the number of parameters.

### Model Selection

Automatically evaluates multiple model types (Random Forest, Gradient Boosting, Logistic Regression, SVM, etc.) with default hyperparameters and selects the best-performing model via cross-validation.

### Auto Preprocessing

Automatic preprocessing pipeline that handles missing values, scaling, encoding, and optional feature selection based on data characteristics.

## Function Reference

### random_search

```python
def random_search(
    model: Any,
    X: Any,
    y: Any,
    param_distributions: dict[str, list[Any]],
    n_iter: int = 50,
    cv: int = 5,
    scoring: str | Callable | None = None,
    random_state: int | None = None,
) -> dict[str, Any]
```

Perform random hyperparameter search. Returns `best_params`, `best_score`, `all_results` (list of param/score pairs), and `elapsed_time`.

### bayesian_optimization

```python
def bayesian_optimization(
    model: Any,
    X: Any,
    y: Any,
    param_space: dict[str, tuple[float, float]],
    n_iterations: int = 30,
    n_initial: int = 5,
    cv: int = 5,
    scoring: str | Callable | None = None,
    random_state: int | None = None,
) -> dict[str, Any]
```

Bayesian optimization with a Gaussian process surrogate. `param_space` maps parameter names to (min, max) bounds for continuous parameters. Returns `best_params`, `best_score`, `convergence_history`, and `surrogate_model`.

### grid_search

```python
def grid_search(
    model: Any,
    X: Any,
    y: Any,
    param_grid: dict[str, list[Any]],
    cv: int = 5,
    scoring: str | Callable | None = None,
) -> dict[str, Any]
```

Exhaustive grid search over all parameter combinations. Returns `best_params`, `best_score`, `all_results`, and `n_combinations`.

### model_selection

```python
def model_selection(
    X: Any,
    y: Any,
    task: str = "classification",
    cv: int = 5,
    scoring: str | None = None,
    random_state: int | None = None,
) -> dict[str, Any]
```

Evaluate multiple model types and select the best. Task must be `"classification"` or `"regression"`. Returns `best_model`, `best_model_name`, `best_score`, and `model_scores` (all models ranked).

### auto_preprocess

```python
def auto_preprocess(
    X: Any,
    y: Any | None = None,
    scale: bool = True,
    handle_missing: bool = True,
    feature_selection: bool = False,
    n_features: int | None = None,
) -> dict[str, Any]
```

Automatic data preprocessing pipeline. Returns `X_processed`, `y_processed` (if label encoding applied), `transformers` (fitted scaler, imputer, etc.), and `feature_mask` (if selection applied).

## Usage Examples

```python
from metainformant.ml import (
    random_search, bayesian_optimization, grid_search,
    model_selection, auto_preprocess,
)
from sklearn.ensemble import RandomForestClassifier

# Random search for hyperparameter tuning
rf = RandomForestClassifier(random_state=42)
param_dists = {
    "n_estimators": [50, 100, 200, 500],
    "max_depth": [5, 10, 20, None],
    "min_samples_split": [2, 5, 10],
}
result = random_search(rf, X, y, param_dists, n_iter=50, cv=5)
print(f"Best params: {result['best_params']}, Score: {result['best_score']:.3f}")

# Bayesian optimization for continuous parameters
param_space = {
    "max_depth": (2.0, 30.0),
    "min_samples_split": (2.0, 20.0),
}
bayes_result = bayesian_optimization(rf, X, y, param_space, n_iterations=30)
print(f"Best score: {bayes_result['best_score']:.3f}")

# Grid search (exhaustive)
param_grid = {"n_estimators": [100, 200], "max_depth": [10, 20]}
grid_result = grid_search(rf, X, y, param_grid, cv=5)

# Automatic model selection
selection = model_selection(X, y, task="classification", cv=5)
print(f"Best model: {selection['best_model_name']} (score: {selection['best_score']:.3f})")
best_model = selection["best_model"]

# Auto preprocessing
prep = auto_preprocess(X, y, scale=True, handle_missing=True, feature_selection=True)
X_clean = prep["X_processed"]
```

## Configuration

- **Environment prefix**: `ML_`
- **Required**: numpy
- **Optional**: scikit-learn (for cross-validation and model training)
- Bayesian optimization uses a pure Python GP surrogate with no external dependencies
- Model selection evaluates Random Forest, Gradient Boosting, Logistic Regression, and SVM by default
- All methods support custom scoring functions

## Related Modules

- `ml.models` -- Model classes used in selection and tuning
- `ml.evaluation` -- Cross-validation methods used internally by AutoML
- `ml.features` -- Feature selection as part of preprocessing
- `ml.interpretability` -- Interpret selected models
