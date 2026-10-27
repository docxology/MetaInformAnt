# Specification: models

## 🎯 Scope
ML models sub-package (classification, regression).

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 3 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `classification.py`
- `regression.py`

### Behavior Contracts

#### `classification.py`
- `cross_validate_biological(X, y, method="rf", cv_folds=5, random_state=None, algorithm=None, parallel_jobs=-1, **kwargs)`
  — `method` is `'rf'`, `'gb'`, `'lr'` or `'ensemble'` (`algorithm` is an alias using the long names). The folds are
  `StratifiedKFold(cv_folds, shuffle=True, random_state=random_state)`. For `'ensemble'`, `random_state` goes to
  `train_ensemble_classifier` as its own argument and not also through the member options, so a seeded ensemble run
  works. `parallel_jobs` is the `n_jobs` of the cross-validation and, for `'ensemble'`, of the voting ensemble unless
  `n_jobs` is passed in `**kwargs`; `parallel_jobs=1` starts no worker processes.

#### `regression.py`
- `analyze_prediction_uncertainty(X, y, X_test, method="rf", n_bootstraps=100, random_state=None)` draws its bootstrap
  samples from `numpy.random.RandomState(random_state)` and trains each regressor with the same `random_state`; the same
  seed gives the same intervals, and the global NumPy random state is neither seeded nor consumed.
