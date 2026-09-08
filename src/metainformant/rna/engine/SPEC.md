# Specification: engine

## 🎯 Scope
Current RNA/Amalgkit execution engine.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 18 Python modules
- **Key Concepts**: Refer to Pydantic models in source. `progress_db.py`
  owns the sample state machine (`pending → downloading → downloaded →
  quantifying → quantified`, with `failed` branches and provenance-audit
  `quarantined` transitions) and the
  failure-class contract: `classify_sample_error()` maps stored
  `samples.error` text via first-match order-sensitive markers
  (`SAMPLE_ERROR_CLASSES`) to `environment_write_denied`,
  `environment_missing_tool`, `transfer_all_sources_failed`,
  `extraction_timeout`, `quantification_timeout`, `quantification_failed`,
  `quantification_exception`, `unrecorded` (no error stored), or
  `unclassified` (no marker matched); the campaign status report surfaces
  the counts as `db_failure_classes`, and
  `metainformant.rna.analysis.cohort_accounting.build_cohort_funnel()`
  consumes the DB read-only for fail-closed cohort funnel accounting with
  durable failure-class reason codes. `streaming_orchestrator.py` runs the
  mandatory campaign preflight before discovery and supports the opt-in
  quant stall watchdog (`AMALGKIT_PIPELINE_QUANT_STALL_TIMEOUT_SECONDS`,
  default 0/disabled) around each quantification subprocess.

## 🔌 API Definition
### Exports
- `__init__.py`
- `discovery.py`
- `exclusions.py`
- `pipeline.py`
- `preflight.py`
- `progress_dashboard.py`
- `progress_db.py`
- `provenance.py`
- `raw_cleanup.py`
- `species.py`
- `sra_extraction.py`
- `streaming_orchestrator.py`
- `workflow.py`
- `workflow_cleanup.py`
- `workflow_core.py`
- `workflow_execution.py`
- `workflow_planning.py`
- `workflow_steps.py`
