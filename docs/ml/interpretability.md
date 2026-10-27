# Model Interpretability

The interpretability module provides model explanation methods for understanding which features drive predictions in biological machine learning models. It implements permutation importance, Kernel SHAP, LIME, partial dependence, feature interactions and attention-weight extraction, plus interpretability-based feature selection (Boruta, recursive elimination, stability selection, mutual information). Inputs can be NumPy arrays or plain nested lists, and the methods are implemented in pure Python on top of your model's `predict`.

## Key Concepts

### Permutation Importance

Measures feature importance by randomly shuffling each feature and measuring the decrease in model performance. Features whose shuffling causes large performance drops are important. Model-agnostic -- works with any model that has a `predict` method.

### Kernel SHAP

Approximation of Shapley values using a weighted linear regression approach. SHAP values explain the contribution of each feature to each individual prediction, decomposing the prediction into a sum of feature effects.

### LIME (Local Interpretable Model-agnostic Explanations)

Explains individual predictions by fitting a local interpretable model (a weighted linear regression) around the point of interest. Perturbs the input and observes how predictions change.

### Partial Dependence

Shows the marginal effect of one feature (`partial_dependence`, with individual conditional expectation curves) or two features (`feature_interaction`) on the predicted outcome, averaging over all other features. Reveals the functional relationship between a feature and the model's predictions.

### Feature Interaction

Measures the interaction between a pair of features as the variance of the part of their joint partial dependence that is not explained by the two one-feature partial dependences. Identifies features that have synergistic effects on predictions.

### Reproducibility

Permutation importance, Kernel SHAP, LIME, Boruta and stability selection are randomised and take a `random_state`: an integer seed or a `random.Random`. Their randomness (the shuffles, the sampled coalitions, the perturbed samples, the shadow features and the subsamples) comes from the generator that `metainformant.ml._numeric.resolve_python_rng` builds from it. The same seed gives the same result, `random_state=None` draws fresh entropy, and the global `random` state is neither seeded nor consumed. A seed passed to Boruta, `0` included, makes the selection reproducible. Partial dependence and feature interaction are deterministic.

## Function Reference

The explainers live in `metainformant.ml.interpretability.explainers` and the feature selection methods in `metainformant.ml.interpretability.feature_selection`; neither is re-exported from `metainformant.ml`.

### compute_permutation_importance

```python-snippet
def compute_permutation_importance(
    model: Any,
    X: Any,
    y: Any,
    n_repeats: int = 10,
    metric: str = "accuracy",
    random_state: int | random.Random | None = None,
) -> dict: ...
```

Compute permutation importance for all features of a fitted model. `metric` is `'accuracy'`, `'mse'` (reported as negative MSE, so higher is better) or `'r2'`. Returns a dictionary with `importances_mean` and `importances_std` (one value per feature), `feature_names` (`"feature_0"`, `"feature_1"`, ...) and `baseline_score`.

### compute_shap_values_kernel

```python-snippet
def compute_shap_values_kernel(
    predict_fn: Any,
    X: Any,
    n_samples: int = 100,
    background: Any | None = None,
    random_state: int | random.Random | None = None,
) -> dict: ...
```

Compute Kernel SHAP values for each row of `X`. `predict_fn` is a callable that takes a 2D array and returns one prediction per row (for a classifier pass a wrapper such as `lambda A: model.predict_proba(A)[:, 1]`). `n_samples` is the number of coalition samples per instance and `background` defaults to the mean of `X`; `random_state` seeds the sampled coalitions. Returns `shap_values` (list of lists, n_instances x n_features), `expected_value` (base prediction) and `feature_names`.

### compute_lime_explanation

```python-snippet
def compute_lime_explanation(
    predict_fn: Any,
    instance: list[float],
    feature_names: list[str],
    n_samples: int = 1000,
    n_features: int = 10,
    kernel_width: float | None = None,
    random_state: int | random.Random | None = None,
) -> dict: ...
```

Generate a LIME explanation for a single instance (a plain list of feature values; `feature_names` is required). `random_state` seeds the perturbed samples. Returns `coefficients` (a list of `(feature_name, coefficient)` tuples for the top `n_features`), `intercept`, `local_prediction`, `r_squared` (fidelity of the local model) and `feature_contributions` (every feature name mapped to its coefficient).

### partial_dependence

```python-snippet
def partial_dependence(
    model: Any,
    X: Any,
    feature: int,
    n_grid: int = 50,
) -> dict: ...
```

Compute partial dependence of a single feature. Returns `grid_values` (grid of feature values), `pdp_mean` (average prediction at each grid point), `pdp_std` and `ice_curves` (individual conditional expectation curves for up to 50 samples).

### feature_interaction

```python-snippet
def feature_interaction(
    model: Any,
    X: Any,
    feature_i: int,
    feature_j: int,
    n_grid: int = 20,
) -> dict: ...
```

Measure the interaction of one pair of features. Returns `grid_i`, `grid_j`, `pdp_values` (the n_grid x n_grid two-feature partial dependence) and `interaction_strength` (variance of the interaction term). To compare several pairs, call it once per pair.

### compute_attention_weights

```python-snippet
def compute_attention_weights(
    model: Any,
    X: Any,
) -> dict | None: ...
```

Extract attention weights from a model that provides them, either through a `get_attention_weights(X)` method or an `attention_weights` attribute. Returns `attention_matrix` and `layer_attentions`, or `None` when the model exposes neither (this is the case for scikit-learn models).

### Feature Selection (interpretability-based)

```python-snippet
def boruta_selection(
    X: Any,
    y: Any,
    max_iter: int = 100,
    alpha: float = 0.05,
    random_state: int | random.Random | None = None,
) -> dict: ...

def recursive_elimination(
    model: Any,
    X: Any,
    y: Any,
    n_features: int = 10,
    step: int = 1,
    cv: int = 5,
) -> dict: ...

def stability_selection(
    X: Any,
    y: Any,
    n_bootstrap: int = 100,
    threshold: float = 0.6,
    alpha: float = 0.1,
    random_state: int | random.Random | None = None,
) -> dict: ...

def mutual_information_selection(
    X: Any,
    y: Any,
    n_features: int = 20,
    discrete_target: bool = True,
    n_bins: int = 10,
) -> dict: ...
```

Feature selection methods based on interpretability. Boruta compares each feature with shuffled shadow features and returns `selected`, `tentative`, `rejected` (feature indices) and `importance_history`. Recursive elimination (the only one that takes a model, which needs `fit` and `predict`) returns `selected_features`, `ranking` and `cv_scores`. Stability selection fits L1-penalized models on random half-subsamples and returns `selected_features`, `selection_probabilities` and `threshold_used`. Mutual information selection returns `selected_features`, `mi_scores` and `ranking`.

## Usage Examples

```python
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from metainformant.ml.interpretability.explainers import (
    compute_permutation_importance, compute_shap_values_kernel,
    compute_lime_explanation, partial_dependence, feature_interaction,
)
from metainformant.ml.interpretability.feature_selection import (
    boruta_selection, stability_selection,
)

# Synthetic data: the first two of six features drive the label
rng = np.random.default_rng(0)
X = rng.normal(size=(150, 6))
y = (2 * X[:, 0] + X[:, 1] + rng.normal(scale=0.5, size=150) > 0).astype(int)
X_train, X_test, y_train, y_test = X[:100], X[100:], y[:100], y[100:]
model = RandomForestClassifier(n_estimators=30, random_state=0).fit(X_train, y_train)

# Permutation importance
perm_imp = compute_permutation_importance(model, X_test, y_test, n_repeats=20, random_state=0)
ranking = np.argsort(perm_imp["importances_mean"])[::-1]
for idx in ranking[:5]:
    print(f"{perm_imp['feature_names'][idx]}: importance={perm_imp['importances_mean'][idx]:.4f}")

# SHAP values (predict_fn returns the probability of class 1)
predict_fn = lambda A: model.predict_proba(A)[:, 1]
shap = compute_shap_values_kernel(predict_fn, X_test[:10], n_samples=100, background=X_train[:50], random_state=0)
feature_effects = np.abs(np.array(shap["shap_values"])).mean(axis=0)
print(f"Top feature effect: {feature_effects.max():.4f}")

# LIME explanation for a single sample
feature_names = [f"gene_{i}" for i in range(X.shape[1])]
lime = compute_lime_explanation(predict_fn, X_test[0].tolist(), feature_names, n_features=5, random_state=0)
for name, weight in lime["coefficients"]:
    print(f"{name}: weight={weight:.4f}")

# Partial dependence plot data
pd_result = partial_dependence(model, X_test, feature=0)
# Plot pd_result["grid_values"] vs pd_result["pdp_mean"]

# Feature interactions: one call per pair of features
pairs = [(0, 1), (0, 2), (1, 2)]
strengths = {pair: feature_interaction(model, X_test, *pair, n_grid=10)["interaction_strength"] for pair in pairs}
for pair, score in sorted(strengths.items(), key=lambda item: item[1], reverse=True):
    print(f"Features {pair}: interaction={score:.4f}")

# Boruta feature selection
boruta = boruta_selection(X, y, max_iter=20, random_state=0)
selected = boruta["selected"]

# Stability selection
stable = stability_selection(X, y.astype(float), n_bootstrap=30, threshold=0.6, random_state=0)
print(selected, stable["selected_features"])
```

## Configuration

- **Required**: only the METAINFORMANT package itself; NumPy is optional (inputs may be NumPy arrays or nested lists)
- The methods are configured through their arguments; no environment variables are read
- `compute_permutation_importance`, `partial_dependence` and `feature_interaction` call `model.predict`
- The randomised methods draw from their own generator (see Reproducibility), so a seed makes a run repeatable without touching the global random state
- `compute_shap_values_kernel` and `compute_lime_explanation` take any `predict_fn` that maps a 2D array to one prediction per row, so `predict_proba` works through a small wrapper
- `recursive_elimination` needs a model with `fit` and `predict` (a `TypeError` is raised otherwise)

## Related Modules

- `ml.models` -- Models to interpret
- `ml.features` -- Feature selection using statistical methods
- `ml.evaluation` -- Cross-validation for evaluating feature selection
- `ml.automl` -- Automated model selection
