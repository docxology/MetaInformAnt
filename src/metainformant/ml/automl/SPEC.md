# Specification: automl

## 🎯 Scope
Automated machine learning (AutoML) subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 2 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `optimization.py`

### Behavior Contracts

#### `optimization.py`
- `random_search(model_fn, param_distributions, X, y, n_iter=50, cv=5, metric="accuracy", random_state=None)` and
  `bayesian_optimization(objective_fn, param_space, n_iter=30, n_initial=5, random_state=None)` take
  `random_state: int | random.Random | None`. Every random draw (the parameter samples; the initial points and the 200
  expected-improvement candidates per iteration) comes from `ml._numeric.resolve_python_rng(random_state)`: the same
  seed gives the same search, `None` draws fresh entropy, and the global `random` state is neither seeded nor consumed.
