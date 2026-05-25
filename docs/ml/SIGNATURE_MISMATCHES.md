# Signature Mismatches — Technical Reference for Developers

Use this file to align documentation with source or vice versa.

## AutoML (`automl/optimization.py`)

### random_search

**Documented:**
```python-snippet
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

**Actual Source:**
```python-snippet
def random_search(
    model_fn: Any,                              # Callable factory, NOT instantiated model
    param_distributions: dict,                  # Supports list, dict with low/high/log
    X: Any,
    y: Any,
    n_iter: int = 50,
    cv: int = 5,
    metric: str = "accuracy",                   # Named 'metric' not 'scoring'
    random_state: int | None = None,
) -> dict
```

**Fix Options:**
1. Update docs to match source (recommended)
2. Add wrapper `def random_search(model, X, y, param_distributions, ...)` that internally calls source version with factory

### bayesian_optimization

**Documented:**
```python-snippet
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

**Actual Source:**
```python-snippet
def bayesian_optimization(
    objective_fn: Any,                          # Callable(params_dict) -> float score
    param_space: dict,                          # dict[name] = {'low': float, 'high': float, 'log': bool}
    n_iter: int = 30,                           # Not 'n_iterations'
    n_initial: int = 5,
    cv: int = 5,
    metric: str = "accuracy",                   # Named 'metric'
    random_state: int | None = None,
) -> dict
```

**Note:** Source version does NOT take X, y directly — you must wrap your data in an objective function.

**Fix:** Docs must be rewritten to show objective function pattern.

### grid_search

**Documented:**
```python-snippet
def grid_search(
    model: Any,
    X: Any,
    y: Any,
    param_grid: dict[str, list[Any]],
    cv: int = 5,
    scoring: str | Callable | None = None,
) -> dict[str, Any]
```

**Actual Source:**
```python-snippet
def grid_search(
    model_fn: Any,                              # Callable factory
    param_grid: dict,                           # Exact choices, not distributions
    X: Any,
    y: Any,
    cv: int = 5,
    metric: str = "accuracy",
) -> dict
```

### auto_preprocess

**Documented:**
```python-snippet
def auto_preprocess(
    X: Any,
    y: Any | None = None,
    scale: bool = True,
    handle_missing: bool = True,
    feature_selection: bool = False,
    n_features: int | None = None,
) -> dict[str, Any]
```

**Actual Source (as of line 201):**
```python-snippet
def auto_preprocess(
    X: Any,
    y: Any | None = None,
) -> dict:
    # Implementation truncated in read — check actual file
```

**Likely under-implemented relative to docs.**

---

## Evaluation (`evaluation/validation.py`)

### cross_validate_biological

**Documented:**
```python-snippet
def cross_validate_biological(
    model: Any,
    X: np.ndarray,
    y: np.ndarray,
    cv: int = 5,
    scoring: str | list[str] | None = None,
    stratified: bool = True,
    random_state: int | None = None,
) -> dict[str, Any]
```

**Actual Source:**
- Function is named `cross_validation_scores()` NOT `cross_validate_biological`
- Signature:
```python-snippet
def cross_validation_scores(
    model: Any,
    X: np.ndarray,
    y: np.ndarray,
    cv: int = 5,
    scoring: str | List[str] = "accuracy",
    stratified: bool = True,
    random_state: int | None = None,
) -> Dict[str, np.ndarray]:
```
**Returns:** Dict mapping metric name → array of fold scores (NOT named `mean_accuracy` etc.)

**Alternative source function:**
```python-snippet
def cross_validate(
    model: Any = None,
    X: np.ndarray = None,
    y: np.ndarray = None,
    cv: int = 5,
    scoring: str = "accuracy",
    classifier_func: Any = None,
    cv_folds: int | None = None,
    random_state: int | None = None,
) -> Dict[str, Any]:
```
This is a more general wrapper.

**Recommendation:** Either rename source to match docs or update docs to show actual function names.

### bootstrap_validation

**Documented:** Does not exist in source.

**Closest alternatives:**
- `bootstrap_validate(model_factory, X, y, n_bootstraps, test_size, random_state, **kwargs)`
- `validate_model_stability(model_factory, X, y, n_bootstraps, test_size, random_state, **kwargs)`

Both take `model_factory` (callable) not a fitted model.

### learning_curve_analysis

**Documented:** Does not exist.

**Source has:** `learning_curve(model, X, y, train_sizes=None, cv=5, scoring="accuracy")`

### permutation_test

**Documented:** Does not exist.

**Source does not have permutation test statistical significance function.** Only `permutation_importance_biological()` which computes feature importance, not test of model significance.

---

## Interpretability (`interpretability/explainers.py`, `feature_selection.py`)

### compute_permutation_importance

**Documented:**
```python-snippet
def compute_permutation_importance(
    model: Any,
    X: Any,
    y: Any,
    n_repeats: int = 10,
    scoring: Callable | None = None,      # Named 'scoring'
    random_state: int | None = None,
) -> dict[str, Any]
# Returns: importances, importances_std, ranking, scores_per_feature
```

**Actual Source:**
```python-snippet
def compute_permutation_importance(
    model: Any,
    X: Any,
    y: Any,
    n_repeats: int = 10,
    metric: str = "accuracy",            # Named 'metric' not 'scoring'
    random_state: int | None = None,
) -> dict
# Returns: importances_mean, importances_std, feature_names, baseline_score
# NO 'ranking' key, NO 'scores_per_feature' — different return structure
```

### compute_lime_explanation

**Documented:**
```python-snippet
def compute_lime_explanation(
    model: Any,
    X: Any,
    instance_index: int,
    n_features: int = 10,
    n_samples: int = 1000,
) -> dict[str, Any]
# Returns: feature_weights, intercept, local_prediction, actual_prediction, r_squared
```

**Actual Source (inferred from reading):** Signature and return format differ; LIME implementation uses custom pure-Python approach. Exact signature not fully captured in initial read — verify in source file line 600+.

### feature_interaction

**Documented:**
```python-snippet
def feature_interaction(
    model: Any,
    X: Any,
    feature_pairs: list[tuple[int, int]] | None = None,
) -> dict[str, Any]
```

**Source has additional parameters** (per grep). Check source for exact signature.

### compute_shap_values_kernel

**Close match** — source uses:
```python-snippet
def compute_shap_values_kernel(
    predict_fn: Any,                      # Not 'model' — requires predict callable
    X: Any,
    n_samples: int = 100,
    background: Any | None = None,
) -> dict
```

Docs show `model` as first param but source expects `predict_fn` (model.predict or model.predict_proba).

---

## Interpretability — Feature Selection (`feature_selection.py`)

### boruta_selection

**Documented:**
```python-snippet
def boruta_selection(model, X, y, max_iter=100, random_state=None) -> dict
```

**Actual Source:**
```python-snippet
def boruta_selection(
    X: Any,
    y: Any,
    max_iter: int = 100,
    alpha: float = 0.05,                  # Extra parameter not in docs
    random_state: int | None = None,
) -> dict:
```
**Does NOT take `model` parameter** — builds internal random forest. Docs incorrectly show `model` as first arg.

### recursive_elimination

**Documented generically:**
```python-snippet
def recursive_elimination(model, X, y, n_features=None, cv=5) -> dict
```

**Actual Source:**
```python-snippet
def recursive_elimination(
    model: Any,
    X: Any,
    y: Any,
    n_features: int = 10,
    step: int = 1,                         # Missing in docs
    cv: int = 5,
) -> dict
```
Additional `step` parameter controls elimination step size.

### mutual_information_selection

**Documented:** Simple signature
**Source:** Has more complex implementation with internal `_mi_continuous`, `_mi_discrete_continuous` helpers. Return format is dictionary.

---

## Models (`models/classification.py`)

### BiologicalClassifier.feature_importances()

**Documented:** Method name `feature_importances()`
**Actual:** Method is `get_feature_importance()` → returns `np.ndarray`, not dict.

**Recommendation:** Either rename method in source or update docs.

### train_classifier convenience function

**Documented in examples:** `classification.train_classifier(features, labels, method=...)`
**Actual:** No such function exists in `classification.py`.

**Alternatives:**
- `train_ensemble_classifier(X, y, n_estimators, ...)` — returns BiologicalClassifier
- `create_biological_classifier(method, **kwargs)` — factory, still need to `.fit()`
- Direct instantiation: `BiologicalClassifier(algorithm=method).fit(X, y)`

**Recommendation:** Add `train_classifier(X, y, method="rf", **kwargs)` convenience wrapper, or update all docs to use explicit class.

---

## LLM (`llm/ollama/`)

### Class Names

**Documented:** `LLMConfig`
**Actual:** `OllamaConfig`

**Fix:** Rename in docs or add alias `LLMConfig = OllamaConfig` in `llm/__init__.py`.

### Missing Methods on OllamaClient

All domain-specific methods are missing:

| Documented Method | Likely Implementation Location |
|---|---|
| `query(prompt, context)` | Could add wrapper around `generate()` |
| `query_biological(query, data, annotations)` | Needs implementation |
| `query_variants(query, vcf_data)` | Needs implementation |
| `generate_code(description, library)` | Needs implementation using prompt templates |
| `explain_code(code)` | Needs implementation |
| `generate_docs(module, format)` | Needs implementation |
| `generate_report(analysis_name, results)` | Needs implementation |
| `interpret_gwas(results)` | Needs implementation |
| `explain_de_results(results)` | Needs implementation |

**Recommendation:** Either implement these high-level helpers in `client.py` or remove their usage from documentation entirely.

---

## Return Value Mismatches

| Function | Docs Promise | Source Returns | Impact |
|---|---|---|---|
| `cross_validate_biological` | Named keys `mean_accuracy`, `std_accuracy` | Raw score arrays per metric | Need to compute mean/std from arrays |
| `permutation_importance` | `importances`, `ranking`, `scores_per_feature` | `importances_mean`, `importances_std`, `feature_names` | Key names differ; no ranking computed |
| `model_selection` | `best_model` (fitted model), `model_scores` (ranked) | Returns dict, check exact structure | Verify `best_model` is actually fitted model object |

---

## Import Path Issues

Docs sometimes use:

```python-snippet
from metainformant.ml import classification
```

But `ml/__init__.py` only exports submodules (`automl`, `deep_learning`, `evaluation`, `features`, `interpretability`, `llm`, `models`). It does NOT import functions from submodules into top-level.

To use classification:
```python
from metainformant.ml.models import classification
from metainformant.ml.models.classification import BiologicalClassifier
```

Not:
```python-snippet
from metainformant.ml import classification  # ✗ This fails unless added to __init__
```

**Fix**: Either update `ml/__init__.py` to re-export commonly-used functions or correct all import examples.

---

## Priorities for Fixing

**P0 (Critical — breaks all examples):**
1. Deep learning API completely missing
2. LLM high-level methods missing
3. `train_classifier()` missing from models
4. `cross_validate_biological` vs `cross_validation_scores` naming mismatch

**P1 (High — breaks specific examples):**
5. `bootstrap_validation`, `learning_curve_analysis`, `permutation_test` missing from evaluation
6. `permutation_importance` location and return format wrong
7. `feature_importances()` method name mismatch
8. `explain` module doesn't exist

**P2 (Medium — documentation gaps):**
9. Additional useful functions undocumented
10. Parameter name mismatches (`metric` vs `scoring`)
11. Import path clarity

---

End of technical reference.
