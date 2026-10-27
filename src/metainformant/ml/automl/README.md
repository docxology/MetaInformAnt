# AutoML

Automated machine learning capabilities including hyperparameter optimization, model selection, and preprocessing for bioinformatics ML pipelines.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `optimization` module |
| `optimization.py` | Random search, Bayesian optimization, grid search, model selection, auto-preprocessing |

## Key Functions

| Function | Description |
|----------|-------------|
| `random_search()` | Hyperparameter optimization via random sampling (`random_state` seeds the parameter draws) |
| `bayesian_optimization()` | Bayesian optimization with Gaussian process surrogate (`random_state` seeds the initial points and the expected-improvement candidates) |
| `grid_search()` | Exhaustive grid search over parameter space |
| `model_selection()` | Automatic model selection across multiple algorithm families |
| `auto_preprocess()` | Automatic preprocessing pipeline construction |

## Usage

```python
from metainformant.ml.automl import optimization

# model_fn builds a model from keyword arguments; param_distributions maps each name to a list of values
# or to a {"low": ..., "high": ...} range
best = optimization.random_search(model_fn, param_distributions, X, y, n_iter=50, random_state=0)

# objective_fn takes a params dict and returns the score to maximise (it holds the data itself);
# param_space maps each name to {"low": ..., "high": ..., "log": False}
best = optimization.bayesian_optimization(objective_fn, param_space, n_iter=30, random_state=0)

selected = optimization.model_selection(X, y, task="classification")
prep = optimization.auto_preprocess(X, y)  # a dict; the processed matrix is prep["X_processed"]
```

## Reproducibility

`random_search` and `bayesian_optimization` take `random_state`, an integer seed or a `random.Random`. Their randomness
comes from the generator that `metainformant.ml._numeric.resolve_python_rng` builds from it: the same seed gives the
same search, `None` draws fresh entropy, and the global `random` state is neither seeded nor consumed. The seed does
not reach the models themselves; give a stochastic model its own `random_state`.
