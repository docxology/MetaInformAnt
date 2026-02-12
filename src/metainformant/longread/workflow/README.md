# Workflow

End-to-end pipeline orchestration for long-read sequencing analysis, providing DAG-based pipeline execution, pre-defined pipeline configurations, and structured report generation.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports orchestrator, orchestrator_core, pipeline_stages, pipelines, reporting |
| `orchestrator.py` | Re-export module aggregating public API from core and stages |
| `orchestrator_core.py` | LongReadOrchestrator class with DAG execution engine |
| `pipeline_stages.py` | Step builders for QC, assembly, methylation, and SV pipelines |
| `pipelines.py` | Pre-defined pipeline configurations with validation |
| `reporting.py` | QC, assembly, and methylation report generation (JSON/text/HTML) |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `LongReadOrchestrator` | Pipeline execution engine with dependency resolution |
| `PipelineStep` | Single step in a pipeline DAG |
| `PipelineResult` | Aggregated results from a pipeline run |
| `pipelines.get_qc_pipeline_config()` | QC pipeline configuration factory |
| `pipelines.get_assembly_pipeline_config()` | Assembly pipeline configuration factory |
| `pipelines.get_methylation_pipeline_config()` | Methylation pipeline configuration factory |
| `pipelines.get_sv_pipeline_config()` | SV calling pipeline configuration factory |
| `pipelines.validate_pipeline_config()` | Validate pipeline configuration parameters |
| `reporting.generate_qc_report()` | Generate structured QC report from pipeline results |
| `reporting.export_report()` | Export report as JSON, text, or HTML |

## Usage

```python
from metainformant.longread.workflow import orchestrator, pipelines, reporting

config = pipelines.get_qc_pipeline_config(min_length=1000, min_quality=7.0)
orch = orchestrator.LongReadOrchestrator(config=config, output_dir="output/lr")
result = orch.run_qc_pipeline(reads)
reporting.export_report(result, "output/lr/qc_report.json", fmt="json")
```
