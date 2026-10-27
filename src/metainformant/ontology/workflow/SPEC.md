# Specification: workflow

## 🎯 Scope
Ontology workflow entry points.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 2 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `run_ontology.py`

### Behavior Contracts

#### `run_ontology.py`
- `run_ontology_analysis(config, phenotype, model, project_root)` reads its settings from `config["ontology"]` and its
  paths from `config["paths"]` (`results_dir`, default `output`). See `README.md` for the keys and their defaults.
- `ontology.random_state` — read by `_gsea_random_state`: unset or `None` gives `None` (unseeded permutations), any other
  value goes through `int()` (so `"11"` gives 11 and `0` is a valid seed). The result is passed to
  `gsea(..., random_state=...)`, so the same seed, ranking and gene sets reproduce the GSEA p-values, NES and FDR
  without touching the global `random` state.
- The GSEA step runs `min(ontology.gsea_permutations, 200)` permutations per gene set.
