# Specification: evaluation

## 🎯 Scope
ML evaluation sub-package (validation, cross-validation).

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 2 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `validation.py`

### Behavior Contracts

#### `validation.py`
- `compare_validation_strategies(X, y, model_factory, strategies=None, random_state=None)` — `strategies` defaults to
  `["holdout", "5fold_cv", "10fold_cv", "stratified_5fold"]`. `"holdout"` is one stratified 80/20 split, `"<k>fold_cv"` is
  `cross_validation_scores(cv=k)`, and `"stratified_5fold"` is `cross_validation_scores(cv=5, stratified=True)`. A strategy
  that fails (an unknown name included) is reported as `{"error": message}` under `comparison`; `best_strategy` is the
  highest-scoring strategy, or `None` when every strategy failed.
- `permutation_importance_biological`, `validate_model_stability`, `bootstrap_validate` and `learning_curve` draw their
  random numbers from `numpy.random.RandomState(random_state)`, local to the call: the same seed gives the same
  result, `None` draws fresh entropy, and the global NumPy random state is neither seeded nor consumed. The shuffled
  folds of the cross-validation helpers come from scikit-learn splitters seeded with `random_state`. The evaluated
  models are not seeded by these functions.
- `cross_validate_biological` lives in `models.classification`, not here; see `models/README.md`.
