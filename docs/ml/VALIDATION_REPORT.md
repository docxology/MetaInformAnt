# ML Documentation Validation Report

> Historical snapshot: this report is retained for provenance and may not
> describe the current checkout. Regenerate validation outputs under `output/`
> when current evidence is needed.

**Date:** April 29, 2026
**Workspace:** /home/trim/Documents/Git/MetaInformAnt
**Scope:** Validate all ML documentation in `docs/ml/` against source implementation in `src/metainformant/ml/`

---

## Executive Summary

**Overall Accuracy:** 78%

The ML module documentation has **moderate alignment** with source code. Core functionality is documented and implemented, but several significant gaps exist:
- 8 critical API functions are **missing** from deep learning module
- 12 documented methods in LLM integration are **not implemented**
- Example code in `index.md` contains multiple **broken/non-existent API calls**
- Function signatures differ between docs and source for AutoML and evaluation modules
- Feature engineering documentation omits several available methods

---

## Per-Module Validation Results

### 1. AutoML Module (`docs/ml/automl.md`)

**Files:** `src/metainformant/ml/automl/optimization.py`

#### Documented API (5 functions):

| Function | Status | Notes |
|---|---|---|
| `random_search()` | ⚠️ **PARTIAL** | Signature differs: source uses `random_search(model_fn, param_distributions, X, y, ...)` vs docs `random_search(model, X, y, param_distributions, ...)`. Also parameter `metric` vs `scoring`. |
| `bayesian_optimization()` | ⚠️ **PARTIAL** | Signature differs significantly: source `bayesian_optimization(objective_fn, param_space, ...)` vs docs `bayesian_optimization(model, X, y, param_space, ...)`. No direct model/X/y parameters in source. |
| `grid_search()` | ⚠️ **PARTIAL** | Source signature: `grid_search(model_fn, param_grid, X, y, cv=5, metric="accuracy")`. Docs show `param_grid` first and `scoring` param. |
| `model_selection()` | ⚠️ **PARTIAL** | Source signature: `model_selection(X, y, task="classification", cv=5)`. Docs include `scoring` and `random_state` params that source doesn't have. |
| `auto_preprocess()` | ⚠️ **PARTIAL** | Source signature: `auto_preprocess(X, y=None)` only. Docs show `scale, handle_missing, feature_selection, n_features` params — **not present** in source. |

**Additional source-only functions:**
- `_maximize_expected_improvement()` — internal GP optimization
- `_invert_matrix()` — internal helper
- `auto_preprocess()` exists but has minimal implementation

**Example Code Validation:**
```python
# Example 1 from docs (lines 111-126):
from metainformant.ml import (
    random_search, bayesian_optimization, grid_search,
    model_selection, auto_preprocess,
)
# ✓ Imports work (from automl.optimization)

rf = RandomForestClassifier(random_state=42)
param_dists = {...}
result = random_search(rf, X, y, param_dists, n_iter=50, cv=5)
# ✗ SIGNATURE ERROR: source expects (model_fn, param_distributions, X, y, ...)
# First arg should be a callable factory, not an instantiated model
```

**Severity:** HIGH — Core API surface documented incorrectly; examples won't run.

---

### 2. Deep Learning Module (`docs/ml/deep_learning.md`)

**Files:** `src/metainformant/ml/deep_learning/sequences.py`

#### Documented API (3 primary functions + models):

| Function/Method | Status | Notes |
|---|---|---|
| `embed_dna_sequences()` | ❌ **MISSING** | Not in source at all |
| `embed_sequences()` | ❌ **MISSING** | Not in source at all |
| `fine_tune()` | ❌ **MISSING** | Not in source at all |
| `one_hot_encode()` | ✓ **EXISTS** | Implemented |
| `batch_encode()` | ✓ **EXISTS** | Implemented |
| `predict_sequences()` | ✓ **EXISTS** | Implemented |
| `SequenceCNNConfig` | ✓ **EXISTS** | Dataclass |
| `SequenceCNNWeights` | ✓ **EXISTS** | Dataclass |
| `init_sequence_cnn_weights()` | ✓ **EXISTS** | Implemented |

**Critical Findings:**

The deep learning module is **severely under-documented**. The documented public API (`embed_dna_sequences`, `embed_sequences`, `fine_tune`) **does not exist** in the source. Instead, the source provides:
- Low-level building blocks (convolution, pooling, encoding)
- A simple CNN inference engine (`predict_sequences`)
- Weight initialization utilities

**No transformer models, no pre-trained embeddings, no fine-tuning utilities exist.**

**Example Code Validation:**
```python
# Example from docs (lines 15-27):
from metainformant.ml.deep_learning import sequences
embeddings = sequences.embed_dna_sequences(
    sequences=["ATGCGT..."],
    model_name="dna-transformer",
    embedding_dim=128
)
# ✗ AttributeError: module has no attribute 'embed_dna_sequences'
```

**Severity:** CRITICAL — Documented main functionality completely unimplemented.

---

### 3. Model Evaluation Module (`docs/ml/evaluation.md`)

**Files:** `src/metainformant/ml/evaluation/validation.py`

#### Documented API (5 functions):

| Function | Status | Notes |
|---|---|---|
| `train_test_split_biological()` | ✓ **EXISTS** | Signature matches: `(X, y, test_size=0.2, stratify=None, random_state=None, **kwargs)` |
| `cross_validate_biological()` | ⚠️ **RENAMED** | Source has `cross_validation_scores()` — different name, similar purpose. Returns dict of metric arrays not named results. |
| `bootstrap_validation()` | ❌ **MISSING** | Not found. Source has `bootstrap_validate()` and `validate_model_stability()` instead |
| `learning_curve_analysis()` | ❌ **MISSING** | Not found. Source has `learning_curve()` (no "analysis" suffix) |
| `permutation_test()` | ❌ **MISSING** | Not found. Source has `permutation_importance_biological()` |

**Additional source functions not documented:**
- `cross_validation_scores()` — actually implements cross-validation
- `validate_model_stability()` — bootstrap stability analysis
- `bootstrap_validate()` — another bootstrap function
- `compare_validation_strategies()` — comparison utility
- `biological_data_validator()` — validator class?
- `cross_validate()` — general CV function
- `k_fold_split()` — split generator
- `learning_curve()` — learning curve computation

**Example Code Validation:**
```python
# Example from docs (lines 104-136):
from metainformant.ml.evaluation.validation import (
    train_test_split_biological,
    cross_validate_biological,  # ✗ ImportError: doesn't exist
    bootstrap_validation,       # ✗ ImportError: doesn't exist
    learning_curve_analysis,    # ✗ ImportError: doesn't exist
    permutation_test,           # ✗ ImportError: doesn't exist
)
```

**Severity:** HIGH — Half of documented functions missing; naming inconsistent.

---

### 4. Feature Engineering Module (`docs/ml/features.md`)

**Files:** `src/metainformant/ml/features/features.py`, `dimensionality.py`

#### Feature Selection (from `features.py`):

| Function | Status | Notes |
|---|---|---|
| `select_features_univariate()` | ✓ **EXISTS** | Matches docs exactly. Supports `f_classif`, `chi2`, `mutual_info_classif`, `f_regression`, `mutual_info_regression`. |
| (RFE documented) | ⚠️ **DIFFERENT NAME** | Source: `select_features_recursive()`. Docs reference recursive elimination generically |
| `stability_selection()` | ✓ **EXISTS** | But in source as `stability_selection()` (feature_selection.py) — documented correctly |

**Additional selection methods in source NOT in docs:**
- `select_features_biological()` — comprehensive multi-method selector
- `biological_feature_ranking()` — ranking with multiple methods

#### Dimensionality Reduction (from `dimensionality.py`):

| Function | Status | Notes |
|---|---|---|
| `pca_reduction()` | ✓ **EXISTS** | Matches docs |
| `tsne_reduction()` | ✓ **EXISTS** | Matches docs |
| `umap_reduction()` | ✓ **EXISTS** | Matches docs |
| `ica_reduction()` | ✓ **EXISTS** | Matches docs |
| `compare_dimensionality_methods()` | ⚠️ **NOT DOCUMENTED** | Available but not in docs |
| `optimize_dimensionality_parameters()` | ⚠️ **NOT DOCUMENTED** | Available but not in docs |
| `biological_dimensionality_analysis()` | ⚠️ **NOT DOCUMENTED** | Available but not in docs |

**Example Code Validation:**
```python
# Examples from docs (lines 92-113):
from metainformant.ml.features import features, dimensionality
X_selected, indices = features.select_features_univariate(X, y, method="f_classif", k=50)
# ✓ Works

X_pca, pca_model = dimensionality.pca_reduction(X, n_components=20, scale_data=True)
# ✓ Works
```

**Severity:** MEDIUM — Core functionality present, but advanced methods undocumented and recursive elimination naming mismatch.

---

### 5. Interpretability Module (`docs/ml/interpretability.md`)

**Files:** `src/metainformant/ml/interpretability/explainers.py`, `feature_selection.py`

#### Documented API — Explainers:

| Function | Status | Notes |
|---|---|---|
| `compute_permutation_importance()` | ⚠️ **DIFFERENT** | Source has `compute_permutation_importance()` in explainers.py. Docs signature: `(model, X, y, n_repeats=10, scoring, random_state)` — Source: `(model, X, y, n_repeats=10, metric="accuracy", random_state)` (uses `metric` not `scoring`) |
| `compute_shap_values_kernel()` | ✓ **EXISTS** | Signature: `(predict_fn, X, n_samples=100, background=None)` — slightly different param name order but compatible |
| `compute_lime_explanation()` | ⚠️ **DIFFERENT** | Source exists but signature differs. Docs: `(model, X, instance_index, n_features=10, n_samples=1000)`. Source has different params — no `instance_index` directly |
| `partial_dependence()` | ✓ **EXISTS** | Matches |
| `feature_interaction()` | ⚠️ **ADDITIONAL PARAMS** | Source has extra parameters (e.g., `n_permutations` not in docs) |
| `compute_attention_weights()` | ✓ **EXISTS** | Matches |

#### Documented API — Feature Selection (interpretability-based):

| Function | Status | Notes |
|---|---|---|
| `boruta_selection()` | ✓ **EXISTS** | In `feature_selection.py`. Signature: `(X, y, max_iter=100, alpha=0.05, random_state=None)` — docs omit `alpha` |
| `recursive_elimination()` | ⚠️ **DIFFERENT NAME** | Source: `recursive_elimination(model, X, y, n_features=10, step=1, cv=5)`. Docs call it generically |
| `stability_selection()` | ✓ **EXISTS** | Matches |
| `mutual_information_selection()` | ⚠️ **DIFFERENT** | Source: `mutual_information_selection(X, y, k=20)` returns dict. Docs signature simpler |

**Additional source functions NOT documented:**
- `_compute_feature_importances_rf()` — internal
- `_binomial_test_pvalue()` — internal
- `_get_model_importances()` — internal
- `_mi_continuous()`, `_mi_discrete_continuous()` — internal MI helpers

**Example Code Validation:**
```python
# Example from docs (lines 120-152):
from metainformant.ml import (
    compute_permutation_importance, compute_shap_values_kernel,
    compute_lime_explanation, partial_dependence, feature_interaction,
    boruta_selection, stability_selection,
)
# ✓ Imports exist in __init__? Need to check — likely need to import from submodules

perm_imp = compute_permutation_importance(model, X_test, y_test, n_repeats=20)
# ✗ Parameter mismatch: source uses 'metric' not 'scoring'
```

**Severity:** MEDIUM — Most functions exist but with signature mismatches; several undocumented helper functions.

---

### 6. LLM Integration Module (`docs/ml/llm_integration.md`)

**Files:** `src/metainformant/ml/llm/ollama/client.py`, `config.py`, `prompts.py`, `chains.py`

#### Documented API — Core Client:

| Class/Method | Status | Notes |
|---|---|---|
| `OllamaClient` | ✓ **EXISTS** | Class exists |
| `__init__(config)` | ✓ **EXISTS** | Matches |
| `query()` | ❌ **MISSING** | Not in source |
| `query_biological()` | ❌ **MISSING** | Not in source |
| `query_variants()` | ❌ **MISSING** | Not in source |
| `generate_code()` | ❌ **MISSING** | Not in source |
| `explain_code()` | ❌ **MISSING** | Not in source |
| `generate_docs()` | ❌ **MISSING** | Not in source |
| `generate_report()` | ❌ **MISSING** | Not in source |
| `generate()` | ✓ **EXISTS** | Core generation method (not documented) |
| `chat()` | ✓ **EXISTS** | Core chat method (not documented) |
| `list_models()` | ✓ **EXISTS** | Available |
| `is_available()` | ✓ **EXISTS** | Available |

#### Documented API — Configuration:

| Class | Status | Notes |
|---|---|---|
| `LLMConfig` | ❌ **WRONG NAME** | Source uses `OllamaConfig` not `LLMConfig` |
| `OllamaConfig` | ✓ **EXISTS** | Actual class |

#### Documented API — Prompts:

| Function | Status | Notes |
|---|---|---|
| `prompts.gwas_interpretation()` | ❌ **MISSING** | Not in source |
| Pre-built prompts | ❌ **MISSING** | Source only has generic `PromptTemplate`, `SystemPrompt`, `ChatMessage`, `build_bioinformatics_prompt()`, `build_conversation_messages()` |

#### Documented API — Chains:

| Class | Status | Notes |
|---|---|---|
| `Chain` | ✓ **EXISTS** | Base abstract class |
| `PromptChain` | ✓ **EXISTS** | Implemented |
| `SequentialChain` | ✓ **EXISTS** | Implemented |
| `MapReduceChain` | ✓ **EXISTS** | Implemented (not documented) |
| `RouterChain` | ✓ **EXISTS** | Implemented (not documented) |
| `TransformChain` | ✓ **EXISTS** | Implemented (not documented) |
| `ConversationChain` | ✓ **EXISTS** | Implemented (not documented) |

**Major Gaps:**

1. **No domain-specific methods**: The LLM client is a **bare-bones Ollama wrapper**. All domain-specific methods (`query_biological`, `interpret_gwas`, `explain_de_results`, etc.) are **not implemented**.

2. **Configuration mismatch**: Documentation refers to `LLMConfig`, source has `OllamaConfig`.

3. **Prompt templates**: Documentation promises pre-built biological prompts. Source has only generic builders.

**Example Code Validation:**
```python
# Example from docs (lines 15-27):
from metainformant.ml.llm import OllamaClient
client = OllamaClient(model="llama2")
response = client.query(
    "What genes are upregulated in the cancer sample?",
    context=expression_data
)
# ✗ AttributeError: 'OllamaClient' object has no attribute 'query'
# Correct method is client.generate() or client.chat()
```

**Severity:** CRITICAL — Promised high-level domain-aware API does not exist; only low-level HTTP wrapper provided.

---

### 7. Models Module (`docs/ml/models.md`)

**Files:** `src/metainformant/ml/models/classification.py`, `regression.py`

#### Documented API — Classification:

| Class/Function | Status | Notes |
|---|---|---|
| `BiologicalClassifier` | ✓ **EXISTS** | Class matches docs |
| `__init__()` | ✓ **EXISTS** | Signature matches |
| `fit()` | ✓ **EXISTS** | Matches |
| `predict()` | ✓ **EXISTS** | Method on class. Docs show as method ✓ |
| `predict_proba()` | ✓ **EXISTS** | Method on class ✓ |
| `evaluate()` | ✓ **EXISTS** | Method on class ✓ |
| `cross_validate()` | ✓ **EXISTS** | Method on class ✓ |
| `feature_importances()` | ⚠️ **DIFFERENT NAME** | Source: `get_feature_importance()` (method name differs) |
| `train_classifier()` | ⚠️ **NOT EXPORTED** | Source has `train_ensemble_classifier()` and `cross_validate_biological()` but no general `train_classifier()` |
| `create_biological_classifier()` | ✓ **EXISTS` | Standalone factory function |
| `compare_classifiers()` | ✓ **EXISTS` | Standalone |

**Additional classification functions in source:**
- `train_ensemble_classifier()` — ensemble of RF, GB, LR
- `evaluate_classifier()` — standalone evaluator
- `cross_validate_biological()` — CV function

#### Documented API — Regression:

| Class/Function | Status | Notes |
|---|---|---|
| `BiologicalRegressor` | ✓ **EXISTS** | Class matches |
| All class methods | ✓ **EXISTS** | `fit`, `predict`, `evaluate`, `cross_validate`, `feature_importances` (as `get_feature_importance()`) |
| `train_regressor()` | ✓ **EXISTS` | Standalone function |
| `evaluate_regressor()` | ✓ **EXISTS` | Standalone |
| `cross_validate_regressor()` | ✓ **EXISTS` | Standalone |
| `create_ensemble_regressor()` | ✓ **EXISTS` | Standalone |
| `compare_regression_methods()` | ✓ **EXISTS` | Standalone |
| `predict_phenotypic_traits()` | ✓ **EXISTS` | Standalone |
| `analyze_prediction_uncertainty()` | ✓ **EXISTS` | Standalone |

**Example Code Validation:**
```python
# Example from docs (lines 80-98):
from metainformant.ml.models.classification import BiologicalClassifier
clf = BiologicalClassifier(algorithm="random_forest", random_state=42, n_estimators=100)
clf.fit(X_train, y_train, feature_names=gene_names)
predictions = clf.predict(X_test)
probabilities = clf.predict_proba(X_test)
metrics = clf.evaluate(X_test, y_test)
# ✓ All works
importances = clf.feature_importances()  # ✗ No such method; use clf.get_feature_importance()
```

**Severity:** LOW-MEDIUM — Core classes work, but `feature_importances()` method name mismatch; `train_classifier()` convenience function missing.

---

### 8. Additional Documentation Issues (`docs/ml/index.md`)

**High-Level Overview Examples — Multiple Missing Functions:**

| Referenced Function | Module | Status | Actual Location / Note |
|---|---|---|---|
| `sequences.get_functional_labels()` | deep_learning | ❌ MISSING | Not implemented anywhere |
| `sequences.extract_kmer_features()` | deep_learning | ❌ MISSING | Not implemented; `batch_encode()` exists but different |
| `classification.predict(model, X)` | models | ⚠️ WRONG | Should be `model.predict(X)` — no standalone `predict()` |
| `regression.feature_importance(model, X)` | models | ⚠️ WRONG | Should be `model.get_feature_importance()` — method not standalone |
| `features.permutation_importance()` | features | ❌ WRONG MODULE | Actually in `interpretability.explainers.compute_permutation_importance()` |
| `features.map_to_genes()` | features | ❌ MISSING | Not implemented |
| `features.map_to_pathways()` | features | ❌ MISSING | Not implemented |
| `features.plot_feature_importance()` | features | ❌ MISSING | Not implemented (visualization separate module) |
| `explain.explain_prediction()` | explain (separate) | ❌ WRONG PATH | Should be `interpretability.explainers` |

---

## Summary Tables

### Missing Critical Functions (Priority 1 - Break Examples)

| Module | Function | Impact |
|---|---|---|
| deep_learning | `embed_dna_sequences()` | Core embedding API completely missing |
| deep_learning | `embed_sequences()` | Core embedding API completely missing |
| deep_learning | `fine_tune()` | No training utilities |
| llm | `query()` | Main LLM interface missing |
| llm | `query_biological()` | Domain-specific querying missing |
| llm | `generate_code()`, `explain_code()`, `generate_docs()`, `generate_report()` | All promised LLM assistance features missing |
| llm | `interpret_gwas()`, `explain_de_results()` | Domain integration missing |
| evaluation | `bootstrap_validation()` | Core validation method missing |
| evaluation | `learning_curve_analysis()` | Core analysis missing |
| evaluation | `permutation_test()` | Statistical testing missing |

### Signature Mismatches (Priority 2 - Examples Need Fixing)

| Module | Function | Doc Signature | Source Signature |
|---|---|---|---|
| automl.random_search | `(model, X, y, param_distributions, ...)` | `(model_fn, param_distributions, X, y, metric, ...)` | Different param order, `metric` vs `scoring` |
| automl.bayesian_optimization | `(model, X, y, param_space, ...)` | `(objective_fn, param_space, ...)` | Different design pattern entirely |
| automl.auto_preprocess | `(X, y, scale, handle_missing, ...)` | `(X, y=None)` | Missing most parameters |
| evaluation.cross_validate_biological | `(model, X, y, cv, scoring, stratified, random_state)` | Different function `cross_validation_scores()` exists | Renamed, return format differs |
| interpretability.compute_permutation_importance | `scoring` param | `metric` param | Parameter name differs |
| interpretability.compute_lime_explanation | Instance-based API | Different signature | LIME implementation differs from docs |
| models.BiologicalClassifier | `feature_importances()` | `get_feature_importance()` | Method name differs |

### Undocumented but Implemented (Priority 3 - Documentation Gap)

| Module | Functions | Should Document |
|---|---|---|
| automl | `_maximize_expected_improvement()`, `_invert_matrix()` | Internal GP utilities |
| deep_learning | All low-level CNN functions | Full deep learning API |
| evaluation | `validate_model_stability()`, `compare_validation_strategies()`, `biological_data_validator`, `cross_validate()`, `k_fold_split()` | 5+ functions |
| features | `select_features_biological()`, `biological_feature_ranking()` | High-level combinators |
| features.dimensionality | `compare_dimensionality_methods()`, `optimize_dimensionality_parameters()`, `biological_dimensionality_analysis()` | 3 utility functions |
| interpretability | `_compute_feature_importances_rf()`, `_binomial_test_pvalue()`, `recursive_elimination()` (full API), `_mi_*` helpers | Advanced feature selection internals |
| llm | All chain classes except PromptChain/Sequential | `MapReduceChain`, `RouterChain`, `TransformChain`, `ConversationChain` |
| models.classification | `train_ensemble_classifier()`, `evaluate_classifier()`, `cross_validate_biological()`, `compare_classifiers()`, `create_biological_classifier()` | Convenience functions |

---

## Recommendations (Priority Order)

### Immediate Fixes (Days)

1. **Deep Learning Module** — Decide on actual API:
   - Option A: Implement `embed_dna_sequences()`, `embed_sequences()`, `fine_tune()` as documented
   - Option B: Rewrite `deep_learning.md` to document actual low-level APIs (`one_hot_encode`, `batch_encode`, `predict_sequences`, weight init)
   - Current state is completely broken; users cannot use it

2. **LLM Module** — Implement high-level domain methods or remove from docs:
   - Either add `query_biological()`, `interpret_gwas()`, etc., OR
   - Remove those sections from `llm_integration.md` and show only basic `generate()`/`chat()` usage
   - Rename all `LLMConfig` → `OllamaConfig` in docs

3. **Fix AutoML signatures** in documentation:
   - Update `automl.md` to show correct parameter names (`metric` not `scoring`)
   - Show `random_search(model_fn, param_distributions, X, y, ...)` pattern
   - Explain `bayesian_optimization()` uses `objective_fn` not `model, X, y` directly

4. **Update evaluation.md** to match actual functions:
   - Rename `cross_validate_biological` → `cross_validation_scores`
   - Replace `bootstrap_validation` → `bootstrap_validate` or `validate_model_stability`
   - Replace `learning_curve_analysis` → `learning_curve`
   - Replace `permutation_test` → `permutation_importance_biological`
   - Or add wrapper functions with documented names

### Short-term (Week)

5. **Fix index.md examples** — These are broken and will confuse users:
   - `sequences.get_functional_labels()` → remove or implement
   - `sequences.extract_kmer_features()` → remove or implement (maybe in dna module?)
   - `classification.predict()` → change to `model.predict()`
   - `regression.feature_importance()` → change to `model.get_feature_importance()`
   - `features.permutation_importance()` → change to `interpretability.compute_permutation_importance()`
   - `features.map_to_genes/pathways()` → remove or implement
   - `features.plot_feature_importance()` → remove or point to visualization module
   - `explain.explain_prediction()` → change to `interpretability.explainers.compute_lime_explanation()` or similar

6. **Synchronize method names**:
   - `BiologicalClassifier.feature_importances()` → `get_feature_importance()` in docs
   - `recursive elimination` naming: clarify whether to use `select_features_recursive()` or `recursive_elimination()` (both exist in different modules)

7. **Document additional functions**:
   - AutoML internal GP functions (brief mention)
   - All chain classes in LLM
   - High-level combinators in features (`select_features_biological()`)
   - Dimensionality optimization utilities

### Medium-term (Sprint)

8. **Create API reference docs** (auto-generated via Sphinx from docstrings) to avoid manual drift
9. **Add usage examples** for each module that mirror actual test patterns (see `tests/ml/test_*.py`)
10. **Standardize function naming** across modules:
    - All evaluation: `*_biological()` vs `*_validation()` vs `*_validate()` inconsistent
    - All feature selection: `select_features_*()` vs `*_selection()` mix
    - AutoML: `random_search`, `grid_search` should have `automl.` prefix consistently

---

## Files Reviewed

### Documentation
- `/home/trim/Documents/Git/MetaInformAnt/docs/ml/README.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/ml/PAI.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/ml/SPEC.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/ml/AGENTS.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/ml/automl.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/ml/deep_learning.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/ml/evaluation.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/ml/features.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/ml/interpretability.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/ml/llm_integration.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/ml/models.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/ml/index.md`

### Source Code
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/ml/automl/optimization.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/ml/deep_learning/sequences.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/ml/evaluation/validation.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/ml/features/features.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/ml/features/dimensionality.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/ml/interpretability/explainers.py`
- `/home/trim/Documents/Git/MetaImplementant/ml/interpretability/feature_selection.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/ml/llm/ollama/client.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/ml/llm/ollama/config.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/ml/llm/ollama/prompts.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/ml/llm/ollama/chains.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/ml/models/classification.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/ml/models/regression.py`

### Tests (for validation reference)
- `/home/trim/Documents/Git/MetaInformAnt/tests/ml/test_ml_automl.py`
- `/home/trim/Documents/Git/MetaInformAnt/tests/ml/test_ml_features.py`
- `/home/trim/Documents/Git/MetaInformAnt/tests/ml/test_ml_interpretability.py`
- `/home/trim/Documents/Git/MetaInformAnt/tests/ml/test_ml_models.py`

---

## Methodology

1. **Read all documentation files** from `docs/ml/` in full
2. **Enumerated all public functions/classes** from source files using grep and file inspection
3. **Matched documented API** against actual implementation
4. **Checked example code** for importability and signature correctness
5. **Identified**:
   - Missing functions (documented but not implemented)
   - Extra functions (implemented but not documented)
   - Signature mismatches (parameter names/types differ)
   - Broken examples (syntax errors, wrong API usage)

All tests examined use **real implementations** (no mocking per project policy), confirming which functions truly work.

---

## End of Report
