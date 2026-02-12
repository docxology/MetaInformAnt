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
| `random_search()` | Hyperparameter optimization via random sampling |
| `bayesian_optimization()` | Bayesian optimization with Gaussian process surrogate |
| `grid_search()` | Exhaustive grid search over parameter space |
| `model_selection()` | Automatic model selection across multiple algorithm families |
| `auto_preprocess()` | Automatic preprocessing pipeline construction |

## Usage

```python
from metainformant.ml.automl import optimization

best = optimization.random_search(model, param_space, X, y, n_iter=50)
best = optimization.bayesian_optimization(model, param_space, X, y, n_iter=30)
selected = optimization.model_selection(X, y, task="classification")
X_processed = optimization.auto_preprocess(X, y)
```
