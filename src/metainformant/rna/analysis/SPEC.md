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
  inferential wrapper, and orthology-presence / species-tree invariants. Records carry optional,
  fail-closed-when-declared reporting bindings (`data_root_snapshot_id`, `cohort_included_count` /
  `cohort_excluded_count`, `artifact_paths`, `metadata_harmonization_review`, `species_tree_source`,
  `species_tree_branch_length_scale`), a frozen `SensitivityAnalysis` registry (`name`,
  `varied_parameter`, `baseline_value`, `varied_values`, `expected_direction`, `notes`) validated by
  `validate_sensitivity_analysis()` and rendered as additive `analysis_provenance_sensitivity_*`
  lines, and the non-analysis roles `stopped`/`unavailable` that record a halted analysis without
  declaring result-implying fields. `cohort_accounting.py` derives the campaign cohort funnel from
  the progress DB and the per-species amalgkit configs: `build_cohort_funnel()` returns a frozen
  `FunnelReport` (stages `configured`, `with_progress`, `quantified_runs`, `failed_runs`,
  `excluded_runs`, `pending_runs`, `active_runs` in `STAGE_NAMES` order; durable failure-class
  `reason_codes`; source `db_path`/`config_dir`), rendered as additive `cohort_funnel_<stage>:
  <count>` lines by `render_funnel_lines()` and as a byte-deterministic two-column TSV by `to_tsv()`.

## 🔌 API Definition
### Exports
- `__init__.py`
- `across_species_orchestrator.py`
- `atlas_plots.py`
- `cohort_accounting.py`
- `conservation_profiles.py`
- `cross_species.py`
- `expression.py`
- `expression_analysis.py`
- `expression_core.py`
- `ortholog_mapping.py`
- `protein_integration.py`
- `statistics_contract.py`
