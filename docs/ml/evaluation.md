# Model Evaluation

The evaluation module provides model validation tools designed for biological data, including stratified train/test splitting, k-fold cross-validation, bootstrap validation, learning curves and permutation-based feature importance.

## Key Concepts

### Stratified Splitting

Biological classification tasks often involve imbalanced classes (rare disease subtypes, uncommon cell types). Stratified splitting ensures each fold or split maintains the original class distribution. The module automatically applies stratification when fewer than 20 unique target values are detected.

### Cross-Validation Strategies

- **Stratified K-Fold**: Preserves class proportions in each fold. Used whenever the target has fewer than 20 unique values (classification-like targets).
- **Standard K-Fold**: Used for targets with 20 or more unique values (regression).

### Bootstrap Validation

Resampling with replacement to check how stable model performance is. Each bootstrap iteration draws a sample with replacement, splits it into train and test parts and scores the model on the test part; the functions report the mean and spread over iterations. Particularly useful for small biological datasets where cross-validation may have high variance.

### Learning Curves

Model performance as a function of training set size, to diagnose overfitting and underfitting and to judge whether more samples would improve performance. `learning_curve` returns the scores; plotting them is left to the caller.

### Reproducibility

Every function that draws random numbers takes a `random_state` and draws from a generator of its own: the same seed gives the same result, and the global NumPy random state is neither seeded nor consumed. `random_state=None` draws fresh entropy and leaves the global stream alone. This covers the feature permutations of `permutation_importance_biological`, the bootstrap samples of `bootstrap_validate` and `validate_model_stability`, the training-set subsamples of `learning_curve`, and the shuffled folds of every cross-validation helper. The functions do not seed the models they evaluate, so give a stochastic model its own `random_state`.

## Function Reference

All functions live in `metainformant.ml.evaluation.validation` and need scikit-learn (they raise `ImportError` without it).

### train_test_split_biological

```python-snippet
def train_test_split_biological(
    X: np.ndarray,
    y: np.ndarray,
    test_size: float = 0.2,
    stratify: np.ndarray | None = None,
    random_state: int | None = None,
    **kwargs: Any,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]: ...
```

Split biological data into train/test sets with automatic stratification for classification targets. Returns (X_train, X_test, y_train, y_test).

### cross_validation_scores

```python-snippet
def cross_validation_scores(
    model: Any,
    X: np.ndarray,
    y: np.ndarray,
    cv: int = 5,
    scoring: str | List[str] = "accuracy",
    stratified: bool = True,
    random_state: int | None = None,
) -> Dict[str, np.ndarray]: ...
```

Cross-validate a scikit-learn estimator with one or several scoring metrics at once. Returns a dictionary that maps each metric name to the array of per-fold scores (an empty array for a metric that could not be computed). Stratified folds are used when `stratified=True` and the target has fewer than 20 unique values.

### cross_validate

```python-snippet
def cross_validate(
    model: Any = None,
    X: np.ndarray | None = None,
    y: np.ndarray | None = None,
    cv: int = 5,
    scoring: str = "accuracy",
    classifier_func: Any = None,
    cv_folds: int | None = None,
    random_state: int | None = None,
) -> Dict[str, Any]: ...
```

Single-metric cross-validation. With an estimator (or a factory that returns one) the result holds `scores`, `mean_score`, `std_score`, `cv_folds` and `scoring`. With `classifier_func(X_train, y_train, X_val, y_val)` returning a dict, the result holds `mean_accuracy`, `std_accuracy`, `fold_results` and `cv_folds`. `cv_folds` is an alias for `cv`.

### bootstrap_validate

```python-snippet
def bootstrap_validate(
    X: np.ndarray,
    y: np.ndarray,
    model_func: Callable[..., Any],
    n_bootstrap: int = 100,
    test_size: float = 0.2,
    random_state: int = 42,
) -> Dict[str, Any]: ...
```

Bootstrap validation of a function `model_func(X_train, y_train, X_test, y_test)` that returns either an array of predictions or a dict (optionally with `predictions`, `mse` and your own metrics). Returns `mean_mse`, `std_mse`, `bootstrap_results` (one entry per iteration), `n_bootstrap` and `test_size`. The module does not compute confidence intervals; derive them from `bootstrap_results` as shown below.

### validate_model_stability

```python-snippet
def validate_model_stability(
    X: np.ndarray,
    y: np.ndarray,
    model_factory: Callable[..., Any],
    n_bootstraps: int = 50,
    test_size: float = 0.2,
    random_state: int | None = None,
    **kwargs: Any,
) -> Dict[str, Any]: ...
```

Bootstrap stability of an estimator built by `model_factory(**kwargs)`. Scores are accuracy for targets with fewer than 20 unique values and R-squared otherwise. Returns `n_bootstraps`, `scores` (`mean`, `std`, `min`, `max`, `values`), `stability_metrics` (`coefficient_of_variation`, `stability_score`) and, when the model exposes `feature_importances_` or `coef_`, `feature_importance_stability`.

### learning_curve

```python-snippet
def learning_curve(
    X: np.ndarray,
    y: np.ndarray,
    model_factory: Any = None,
    train_sizes: np.ndarray | None = None,
    cv: int = 5,
    random_state: int = 42,
    classifier_func: Any = None,
    cv_folds: int | None = None,
) -> Dict[str, Any]: ...
```

Compute learning curves. `train_sizes` are fractions of each training fold; the default is `np.linspace(0.1, 1.0, 10)`. Returns `train_sizes`, `train_scores` and `validation_scores` (one list of per-fold scores for each training size) and `cv_folds`. Scores are accuracy for targets with fewer than 20 unique values and R-squared otherwise.

### permutation_importance_biological

```python-snippet
def permutation_importance_biological(
    model: Any,
    X: np.ndarray,
    y: np.ndarray,
    n_repeats: int = 10,
    random_state: int | None = None,
    scoring: str = "accuracy",
) -> Dict[str, Any]: ...
```

Feature importance by permutation: the drop in cross-validated score when one feature is shuffled. Returns `baseline_score`, `scoring_metric`, `n_repeats`, `feature_importances` (a list of `{"feature", "importance_mean", "importance_std", "rank"}`, best first) and `raw_scores`. It ranks features; it is not a significance test of the model.

### compare_validation_strategies

```python-snippet
def compare_validation_strategies(
    X: np.ndarray,
    y: np.ndarray,
    model_factory: Callable[..., Any],
    strategies: List[str] | None = None,
    random_state: int | None = None,
) -> Dict[str, Any]: ...
```

Score one model under several validation strategies. `model_factory()` returns a fresh model for each strategy, and `strategies` defaults to `["holdout", "5fold_cv", "10fold_cv", "stratified_5fold"]`:

- `"holdout"`: one stratified 80/20 split (`train_test_split_biological(test_size=0.2)`).
- `"<k>fold_cv"` (for example `"5fold_cv"` or `"10fold_cv"`): k-fold cross-validation with `cross_validation_scores(cv=k)`; the folds are stratified for classification-like targets.
- `"stratified_5fold"`: shuffled stratified 5-fold cross-validation (`cross_validation_scores(cv=5, stratified=True)`), mean accuracy over the folds.

Returns `comparison` (each strategy mapped to `{"strategy", "score"}`, or to `{"error": message}` when it failed; an unknown strategy name is reported this way too), `best_strategy` (the strategy with the highest score, `None` when every strategy failed) and `strategies_tested`. The score is the accuracy of a classification-like target (fewer than 20 unique values); the holdout of a target with 20 or more unique values reports R-squared.

### cross_validate_biological

`cross_validate_biological(X, y, method="rf", cv_folds=5, random_state=None, algorithm=None, parallel_jobs=-1, **kwargs)` lives in `metainformant.ml.models.classification`, next to the classifiers. It cross-validates a classifier chosen by name (`'rf'`, `'gb'`, `'lr'` or `'ensemble'`) over shuffled stratified folds.

- `random_state` seeds the folds and the model, for `method="ensemble"` too (an integer `random_state` works with every method).
- `parallel_jobs` is the number of parallel jobs for the cross-validation and, for `method="ensemble"`, for the voting ensemble itself. An explicit `n_jobs` keyword overrides it for the ensemble only. `parallel_jobs=1` starts no worker processes at all, so use it in restricted or resource-sensitive environments where process pools are unavailable.

### k_fold_split and biological_data_validator

`k_fold_split(X, y, k=5, stratify=True, random_state=42)` returns a list of `(train_indices, test_indices)` pairs. `biological_data_validator(X, y, checks=None)` checks for missing values, infinite values, constant features and class balance and returns the per-check results with an overall `passed` flag.

## Usage Examples

```python
import numpy as np
from sklearn.linear_model import LogisticRegression
from metainformant.ml.evaluation.validation import (
    train_test_split_biological,
    cross_validation_scores,
    compare_validation_strategies,
    bootstrap_validate,
    learning_curve,
    permutation_importance_biological,
)
from metainformant.ml.models.classification import BiologicalClassifier, cross_validate_biological

# Synthetic data: 120 samples, 8 features, two informative
rng = np.random.default_rng(0)
X = rng.normal(size=(120, 8))
y = (X[:, 0] + X[:, 1] + rng.normal(scale=0.5, size=120) > 0).astype(int)

# Stratified train/test split
X_train, X_test, y_train, y_test = train_test_split_biological(
    X, y, test_size=0.2, random_state=42
)

# Cross-validation (pass the underlying scikit-learn estimator)
model = BiologicalClassifier(algorithm="random_forest", random_state=42)
cv_results = cross_validation_scores(
    model.model, X, y, cv=5,
    scoring=["accuracy", "f1_weighted"],
    random_state=42,
)
print(f"CV Accuracy: {cv_results['accuracy'].mean():.3f}")

# The same model under several validation strategies (holdout, 5- and 10-fold, stratified 5-fold)
strategies = compare_validation_strategies(
    X, y, lambda: LogisticRegression(max_iter=500), random_state=42
)
print(strategies["best_strategy"], strategies["comparison"]["stratified_5fold"])

# A classifier chosen by name, here the voting ensemble; parallel_jobs=1 starts no worker processes
ensemble_cv = cross_validate_biological(
    X, y, method="ensemble", cv_folds=5, random_state=42, parallel_jobs=1
)
print(f"Ensemble accuracy: {ensemble_cv['accuracy']:.3f}")

# Bootstrap validation: model_func receives (X_train, y_train, X_test, y_test)
def fit_predict(X_tr, y_tr, X_te, y_te):
    predictions = LogisticRegression(max_iter=500).fit(X_tr, y_tr).predict(X_te)
    return {"predictions": predictions, "accuracy": float((predictions == y_te).mean())}

boot = bootstrap_validate(X, y, fit_predict, n_bootstrap=50, random_state=42)
accuracies = np.array([result["accuracy"] for result in boot["bootstrap_results"]])
print(
    f"Accuracy: {accuracies.mean():.3f} "
    f"(2.5-97.5th percentile: {np.percentile(accuracies, 2.5):.3f}-{np.percentile(accuracies, 97.5):.3f})"
)

# Learning curve to check for overfitting
lc = learning_curve(X, y, model_factory=lambda: LogisticRegression(max_iter=500), cv=5)
print(lc["train_sizes"])
print(np.mean(lc["train_scores"], axis=1), np.mean(lc["validation_scores"], axis=1))

# Permutation importance of each feature
importance = permutation_importance_biological(
    LogisticRegression(max_iter=500), X, y, n_repeats=3, random_state=42
)
print(importance["feature_importances"][0])  # most important feature
```

## Configuration

- **Required**: numpy, scikit-learn
- The functions are configured through their arguments; no environment variables are read
- Automatic stratification threshold: fewer than 20 unique target values
- Learning curves use 5-fold cross-validation and ten training sizes (0.1 to 1.0) by default
- `bootstrap_validate` and `validate_model_stability` evaluate each bootstrap sample on a random hold-out split of that sample

## Related Modules

- `ml.models` -- Classification and regression models to evaluate (`cross_validate_biological` there cross-validates a classifier chosen by name; see above)
- `ml.features` -- Feature selection (evaluate with cross-validation)
- `ml.interpretability` -- Model explanation methods
- `ml.automl` -- Automated model selection with built-in evaluation
