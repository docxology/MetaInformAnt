# Specification: rna

## 🎯 Scope
Functionality for rna.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Orchestration Script

## 💾 Data Structures
- **Modules**: 34 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition

### Frozen acquisition

`acquisition.py` delegates `freeze`, `plan`, `estimate`, `quote-aws`, `local`/`worker`
and `aws` to the shared engine. `acquisition_worker.py` uses the same worker CLI.
AWS execution requires an explicit frozen configuration directory. Worker
prerequisites resolve native quantifier and batch identity before acquisition;
unsupported long-read reference bindings fail without changing method or cohort.
See [generic acquisition](../../docs/rna/GENERIC_ACQUISITION.md) for resource
controls, immutable admissions, deadline planning and gross-cost boundaries.

### Exports
- `_setup_utils.py`
- `_verify_utils.py`
- `adopt_existing_fastq.py`
- `adopt_existing_fastq_batch.py`
- `analyze_campaign_rate.py`
- `analyze_processing_times.py`
- `batch_genome_index.py`
- `check_environment.py`
- `check_pipeline_status.py`
- `check_tcc.py`
- ...
