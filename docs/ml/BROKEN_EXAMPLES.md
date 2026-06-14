# Broken Examples — Specific Corrections for index.md

This file provides line-by-line fixes for the broken examples in `docs/ml/index.md`.

## 1. Deep Learning Example (Lines 197-208)

### Current (Broken):
```python-snippet
from metainformant.dna.sequence import core as sequences
from metainformant.ml import classification

# Extract features from DNA sequences
sequences = sequences.read_fasta("sequences.fasta")
sequence_features = sequences.extract_kmer_features(sequences, k=3)

# Classify sequences by function
functional_labels = sequences.get_functional_labels(sequences)
model = classification.train_classifier(sequence_features, functional_labels)

# Predict function for new sequences
new_features = sequences.extract_kmer_features(new_sequences, k=3)
predictions = classification.predict(model, new_features)
```

### Issues:
- `metainformant.dna.sequences` module — check if exists; likely in `dna.sequence` submodule
- `read_fasta()` — location uncertain; probably `dna.io.fasta`
- `extract_kmer_features()` — **NOT FOUND** in source anywhere
- `get_functional_labels()` — **NOT FOUND**
- `classification.train_classifier()` — **NOT FOUND** (source has `train_ensemble_classifier()` or create via `BiologicalClassifier`)
- `classification.predict()` — **NOT FOUND** (use model.predict())

### Recommended Rewrite (Option 1 — using actual deep learning):
```python
# Deep learning not fully implemented; skip or use placeholder
# Either remove this section or:
from metainformant.ml.deep_learning.sequences import batch_encode, predict_sequences
from metainformant.ml.models.classification import BiologicalClassifier

# Encode sequences
X = batch_encode(sequences_list, max_length=1000, alphabet="DNA")
# Need trained weights first; this is low-level
# Instead use classification directly on features
```

### Recommended Rewrite (Option 2 — use actual DNA module):
```python
# If k-mer features are in dna.sequence.kmer (not checked yet):
from metainformant.dna.sequence import kmer
from metainformant.ml.models.classification import BiologicalClassifier

# Extract k-mer features (if available)
# sequence_features = kmer.extract_kmer(sequences, k=3)  # Hypothetical

# Use placeholder since this feature doesn't exist yet
raise NotImplementedError("k-mer extraction not yet implemented")
```

---

## 2. Regression Example (Lines 212-231)

### Current (Broken):
```python-snippet
from metainformant.rna import workflow
from metainformant.ml import regression

# Load expression data
expression_data = workflow.extract_expression_patterns(rna_data)

# Predict continuous phenotypes
phenotype_values = load_continuous_phenotypes()
model = regression.train_regressor(
    expression_data,
    phenotype_values,
    method="xgboost",           # ✗ Not a supported method
    validation="cross_validation"  # ✗ Not a parameter
)

# Feature importance analysis
importance = regression.feature_importance(model, expression_data)
top_genes = importance.nlargest(50)
```

### Issues:
- `workflow.extract_expression_patterns()` — uncertain if exists
- `method="xgboost"` — source supports: `linear`, `rf`, `gb`, `ridge`, `lasso`, `svr`, `elasticnet`. **xgboost not in source**
- `validation="cross_validation"` — **NOT a parameter** of `train_regressor()`
- `regression.feature_importance()` — **NOT a function**; use `model.get_feature_importance()`

### Corrected Version:
```python
from metainformant.ml.models.regression import train_regressor

# Assuming expression_data and phenotype_values are numpy arrays
model = train_regressor(
    expression_data,
    phenotype_values,
    method="rf",  # or "gb", "linear", "ridge", etc.
    # No validation parameter; cross-validate separately
)

# Feature importance
importances = model.get_feature_importance()
# Assuming you have a pandas Series or numpy with gene names:
# top_genes_indices = np.argsort(importances)[-50:][::-1]
```

---

## 3. Network Features Example (Lines 234-250)

### Current (Broken):
```python-snippet
from metainformant.networks import ppi
from metainformant.ml import classification

# Use network features for prediction
network_features = ppi.extract_network_features(interaction_network)

# Combine with expression features
combined_features = pd.concat([expression_features, network_features], axis=1)

# Train integrated model
model = classification.train_classifier(
    combined_features,
    labels,
    feature_selection=True,      # ✗ Not a BiologicalClassifier parameter
    cross_validation=True        # ✗ Not a parameter
)
```

### Issues:
- `feature_selection=True` — `BiologicalClassifier.__init__()` doesn't accept this
- `cross_validation=True` — not a parameter
- Suggests pandas (`pd.concat`) — needs explicit import if used

### Corrected Version:
```python
from metainformant.ml.models.classification import BiologicalClassifier
import numpy as np

# Assuming network_features and expression_features are numpy arrays
combined_features = np.hstack([expression_features, network_features])

model = BiologicalClassifier(algorithm="random_forest", random_state=42)
model.fit(combined_features, labels)

# Cross-validate separately if needed:
from metainformant.ml.models.classification import cross_validate_biological
cv_results = cross_validate_biological(combined_features, labels, method="rf", cv_folds=5)
```

---

## 4. Feature Importance Analysis Example (Lines 264-275)

### Current (Broken):
```python
from metainformant.ml import features

# Analyze feature contributions
importance = features.permutation_importance(model, features, labels)

# Biological interpretation
gene_importance = features.map_to_genes(importance, gene_mapping)
pathway_importance = features.map_to_pathways(gene_importance, pathway_database)

# Visualize importance
features.plot_feature_importance(importance, top_n=20)
```

### Issues:
- `features.permutation_importance()` — **WRONG MODULE**. Should be `interpretability.explainers.compute_permutation_importance()`
- `features.map_to_genes()` — **NOT FOUND** anywhere in source
- `features.map_to_pathways()` — **NOT FOUND**
- `features.plot_feature_importance()` — **NOT FOUND** (plotting separate module)

### Corrected Version:
```python
from metainformant.ml.interpretability.explainers import compute_permutation_importance
import numpy as np

# Compute permutation importance
importance_result = compute_permutation_importance(
    model, X_test, y_test, n_repeats=20, metric="accuracy", random_state=42
)
importances = importance_result["importances_mean"]

# Manual gene mapping if you have gene names:
gene_importance = {gene_name: imp for gene_name, imp in zip(gene_names, importances)}

# Sorting for top genes
top_indices = np.argsort(importances)[::-1][:50]
top_genes = [(gene_names[i], importances[i]) for i in top_indices]

# Plotting — use visualization module or matplotlib:
# from metainformant.visualization import plot_bar
# plot_bar(top_genes[:20])
```

---

## 5. Model Explanation Example (Lines 278-293)

### Current (Broken):
```python-snippet
from metainformant.ml import explain

# Explain individual predictions
explanation = explain.explain_prediction(
    model,
    instance,
    features,
    method="shap"
)

print(f"Prediction: {explanation['prediction']}")
print(f"Top positive features: {explanation['top_positive']}")
print(f"Top negative features: {explanation['top_negative']}")
```

### Issues:
- The old `metainformant.ml` explain import path has no `explain` submodule.
- `explain.explain_prediction()` — **NOT FOUND**

### Corrected Version (using SHAP):
```python
from metainformant.ml.interpretability.explainers import compute_shap_values_kernel

# Get SHAP values for instance(s)
shap_result = compute_shap_values_kernel(
    model.predict if hasattr(model, "predict") else model.predict_proba,
    instance_array,  # 2D array with one or more instances
    n_samples=100,
    background=X_train[:50]  # optional background
)

# shap_result["shap_values"] is a list of lists (n_instances x n_features)
# For single instance i:
instance_shap = shap_result["shap_values"][0]
top_positive_idx = np.argsort(instance_shap)[-5:][::-1]
top_negative_idx = np.argsort(instance_shap)[:5]
```

### Corrected Version (using LIME):
```python
from metainformant.ml.interpretability.explainers import compute_lime_explanation

# For classification
lime_result = compute_lime_explanation(
    model, X_test, instance_index=0, n_features=10, n_samples=1000
)
# Returns: feature_weights, intercept, local_prediction, actual_prediction, r_squared
```

---

## 6. Multi-omics Example (Lines 300-321)

### Current (Seems OK but verify module exists):
```python-snippet
from metainformant.multiomics import integration  # ✗ Does multiomics module exist?
from metainformant.ml import classification

multiomics_data = integration.load_multiomics_data([...])
combined_features = integration.extract_multiomic_features(multiomics_data)
model = classification.train_classifier(
    combined_features, labels,
    method="ensemble",
    feature_selection=True,   # ✗ Not a param
    cross_validation=True     # ✗ Not a param
)
```

**Issues**: `multiomics` module may not exist in same structure; `train_classifier()` doesn't exist; invalid params.

---

## 7. Time Series Example (Lines 324-339)

### Current (Broken):
```python-snippet
from metainformant.ml import regression

time_series_data = load_time_series_expression()
time_features = regression.extract_time_features(time_series_data)  # ✗ NOT FOUND

model = regression.train_regressor(
    time_features, time_targets,
    method="lstm",      # ✗ NOT SUPPORTED: only 'rf', 'gb', 'linear', 'ridge', 'lasso', 'svr'
    time_window=10      # ✗ NOT A PARAMETER
)
```

---

## Summary of Broken References in index.md

| Line(s) | Broken Reference | Replacement |
|---|---|---|
| 201-208 | `sequences.get_functional_labels`, `extract_kmer_features` | **Not implemented** — remove or use placeholder |
| 203, 207 | `classification.train_classifier`, `classification.predict` | Use `BiologicalClassifier` class: `clf = BiologicalClassifier(...).fit(...); preds = clf.predict(X)` |
| 221 | `regression.train_regressor` with `method="xgboost"` | Use `train_regressor(..., method="rf")` |
| 229 | `regression.feature_importance` | `model.get_feature_importance()` |
| 268 | `features.permutation_importance` | `interpretability.explainers.compute_permutation_importance()` |
| 271-272 | `features.map_to_genes/pathways` | **Not implemented** — user must map manually |
| 275 | `features.plot_feature_importance` | **Not implemented** — use matplotlib or visualization module |
| 284 | `explain.explain_prediction` | **Module doesn't exist** — use `compute_shap_values_kernel()` or `compute_lime_explanation()` |
| 305-309 | `multiomics.integration` | Verify module exists; may need adjustment |
| 330 | `regression.extract_time_features` | **Not implemented** |
| 337 | `method="lstm"` | **Not supported** — only classical ML methods |

---

## Action Items for Documentation Team

1. **Rewrite or remove** the entire "With DNA Sequences" example (uses non-existent functions)
2. **Fix regression example**: correct method names, remove unsupported params
3. **Fix feature importance section**: point to `interpretability` module, remove `map_to_*` references
4. **Fix model explanation**: use `interpretability.explainers` instead of `explain`
5. **Time series**: either implement `extract_time_features()` and LSTM support, or remove example
6. **Multi-omics**: verify module exists before keeping example
7. **Standardize imports** throughout: `from metainformant.ml.models.classification import BiologicalClassifier`, not the old top-level classification alias.
