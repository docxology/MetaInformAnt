# Specification: analysis

## 🎯 Scope
RNA analysis modules for expression analysis, QC, and validation.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 18 Python modules
- **Key Concepts**: Refer to Pydantic models in source. `statistics_contract.py` enforces the
  descriptive/inferential boundary: frozen `AnalysisProvenance` records (role-conditional
  multiplicity family/method and tested-feature count — `None`/`'not-applicable'` exactly for
  descriptive roles, declared exactly for inferential roles), fail-closed validation and rendering,
  `attrs["role"] = "descriptive"` markers on cross-species fingerprint outputs, a GATED BH-FDR
  inferential wrapper, and orthology-presence / species-tree invariants.

## 🔌 API Definition
### Exports
- `__init__.py`
- `across_species_orchestrator.py`
- `atlas_plots.py`
- `conservation_profiles.py`
- `cross_species.py`
- `expression.py`
- `expression_analysis.py`
- `expression_core.py`
- `ortholog_mapping.py`
- `protein_integration.py`
- `statistics_contract.py`
