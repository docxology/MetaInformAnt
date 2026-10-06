# Specification: analysis

## 🎯 Scope
RNA analysis modules for expression analysis, QC, and validation.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: Expression, QC, comparative analysis, and statistics contracts
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
  `inferential_comparative.py` implements the gated predeclared inferential
  comparative analysis: study-aware OLS contrast fits, DerSimonian-Laird
  random-effects heterogeneity summaries, bootstrap confidence intervals, and
  BH-FDR exclusively via the gated `declared_inferential_bh_fdr` wrapper.
  `ortholog_mapping.py` additionally provides a fail-closed per-species
  retention audit (`audit_species_retention`, `DEFAULT_MIN_RETENTION`) and the
  versioned `MappingArtifactManifest` v1 with fail-closed write/read round-trips.
  `phylogenetic_comparative.py` implements phylogenetic comparative methods on
  contract-validated trees (rootedness caller-declared, fail-closed): Brownian
  tip covariance, Pagel's-lambda PGLS fitted by REML profile likelihood,
  diagnostics, seeded tree-uncertainty resampling, and deterministic Brownian
  trait simulation; inferential p-values are produced only inside validated
  fits, never for descriptive lanes.

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
- `inferential_comparative.py`
- `phylogenetic_comparative.py`
- `protein_integration.py`
- `statistics_contract.py`

## Numerical input and output contracts

Normalization requires finite nonnegative counts. Standard median-ratio
estimation uses features positive in every library. Differential expression
requires two or more samples per group and unique feature/sample IDs; p-value
adjustment rejects probabilities outside [0, 1] and infinities, while retaining
unscoreable NaNs. PCA rejects invalid dimensions, infinities, duplicated axes,
singleton samples, and absent variation. Missing-value mean imputation is an
explicit option whose policy/count accompany the result.

Within-species orchestration retains the complete expression-sample cohort,
requires explicit and unique metadata run IDs, normalizes library depth before
PCA by default, writes hash-bound method sidecars, and propagates calculation
failures. Convenience run-level DE remains exploratory until independent
biological units and a scientific design are established.

See the [Hymenoptera methods contract](../../../../docs/rna/HYMENOPTERA_METHODS.md).

## Project-adapter extraction contract

- Finalized matrix validation is reusable under `expression_io`; plain/gzip
  reads preserve unique axes across chunks and fail on malformed, nonfinite,
  negative or all-zero sample columns. Native profiles require two samples.
- `ortholog_diagnostics` owns membership cardinality, reconciled drop reasons
  and duplicate-evidence annotations, with precedence version retained explicitly.
- `ortholog_profiles` returns species mean-profile distances and overlap counts.
  It selects the first recorded transcript per mapping cell; insufficient overlap
  and constant profiles remain unavailable. This must not be substituted for
  per-gene correlations over aligned samples. Invalid declared matrix inputs fail.
- `counting_statistics` provides descriptive Wilson intervals; `statistics_io`
  parses predeclared provenance. Neither replaces inferential or biological gates.
- Parent modules must not import nested project scripts. Existing nested helper
  names may re-export parent methods while paths/manifests/artifact assembly stay
  project-specific. The project validator checks selected callable method owners
  and refuses shadowing local implementations or unavailable parent APIs.
