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

```python-snippet
def random_search(
    model_fn: Any,
    param_distributions: dict,
    X: Any,
    y: Any,
    n_iter: int = 50,
    cv: int = 5,
    metric: str = "accuracy",
    random_state: int | None = None,
) -> dict
```

Perform random hyperparameter search. `model_fn` is a callable factory that takes keyword arguments and returns a model with fit/predict methods. Each `param_distributions` value is a list (uniform choice) or a dict with `low`/`high` keys (uniform range), optionally `"log": True` for log-uniform sampling. Returns `best_params`, `best_score`, `all_results` (per-iteration param/score pairs), and `n_evaluations`.

### bayesian_optimization

```python-snippet
def bayesian_optimization(
    objective_fn: Any,
    param_space: dict,
    n_iter: int = 30,
    n_initial: int = 5,
    random_state: int | None = None,
) -> dict
```

Bayesian optimization with a Gaussian process (RBF kernel) surrogate and Expected Improvement acquisition. It does not take `X`/`y`: wrap your data in `objective_fn`, a callable that takes a params dict and returns the score to maximize. `param_space` maps parameter names to `{"low": float, "high": float, "log": bool}` specs. Returns `best_params`, `best_score`, `history` (per-iteration records), and a `surrogate_model` summary.

### grid_search

```python-snippet
def grid_search(
    model_fn: Any,
    param_grid: dict,
    X: Any,
    y: Any,
    cv: int = 5,
    metric: str = "accuracy",
) -> dict
```

Exhaustive grid search over all parameter combinations. `param_grid` maps parameter names to exact value lists. Returns `best_params`, `best_score`, and `all_results` (failed combinations are kept with an `error` entry).

### model_selection

```python-snippet
def model_selection(
    X: Any,
    y: Any,
    task: str = "classification",
    cv: int = 5,
) -> dict
```

Evaluate multiple model types (linear, tree-based, ensemble, KNN, SVM) and rank them by cross-validation score — accuracy for `"classification"`, R-squared for `"regression"`. Returns `best_model_type`, `rankings` (sorted `(model_type, score)` tuples), and `cv_results_per_model`.

### auto_preprocess

```python-snippet
def auto_preprocess(
    X: Any,
    y: Any | None = None,
) -> dict
```

Automatic preprocessing pipeline: detects column data types, imputes missing values, scales numeric features, and encodes categorical-like columns. Returns `X_processed`, `transformations_applied`, and `feature_info`.

## Usage Examples

```python-snippet
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from metainformant.ml.automl.optimization import (
    random_search, bayesian_optimization, grid_search,
    model_selection, auto_preprocess,
)

def make_model(**params):
    return RandomForestClassifier(random_state=42, **params)

# Random search: pass a model factory, not an instantiated model
param_dists = {
    "n_estimators": [50, 100, 200, 500],
    "max_depth": [5, 10, 20, None],
    "min_samples_split": [2, 5, 10],
}
result = random_search(make_model, param_dists, X, y, n_iter=50, cv=5, metric="accuracy")
print(f"Best params: {result['best_params']}, Score: {result['best_score']:.3f}")

# Bayesian optimization: wrap the data in an objective function
param_space = {
    "max_depth": {"low": 2.0, "high": 30.0},
    "min_samples_split": {"low": 2.0, "high": 20.0},
}
def objective(params):
    return cross_val_score(make_model(**params), X, y, cv=5).mean()

bayes_result = bayesian_optimization(objective, param_space, n_iter=30, n_initial=5)
print(f"Best score: {bayes_result['best_score']:.3f}")

# Grid search (exhaustive)
param_grid = {"n_estimators": [100, 200], "max_depth": [10, 20]}
grid_result = grid_search(make_model, param_grid, X, y, cv=5)

# Automatic model selection
selection = model_selection(X, y, task="classification", cv=5)
print(f"Best model type: {selection['best_model_type']}")
print(selection["rankings"])

# Auto preprocessing
prep = auto_preprocess(X, y)
X_clean = prep["X_processed"]
```

## Configuration

- **Required**: numpy
- **Optional**: scikit-learn (for cross-validation and model training)
- Bayesian optimization uses a pure Python GP surrogate with no external dependencies
- Model selection evaluates linear, tree-based, ensemble, KNN, and SVM model types (scikit-learn-backed when available)
- Scoring is selected with the `metric` argument (`"accuracy"`, `"mse"`, `"r2"`); there is no custom-callable scoring parameter

## Related Modules

- `ml.models` -- Model classes used in selection and tuning
- `ml.evaluation` -- Cross-validation methods used internally by AutoML
- `ml.features` -- Feature selection as part of preprocessing
- `ml.interpretability` -- Interpret selected models
