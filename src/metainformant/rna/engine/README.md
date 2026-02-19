# Engine Module

Workflow execution, monitoring, and orchestration for RNA-seq pipelines.

## 📦 Components

| File | Purpose |
|------|---------|
| [`streaming_orchestrator.py`](streaming_orchestrator.py) | **Production pipeline** — multi-species, parallel, ENA-first |
| [`orchestrator.py`](orchestrator.py) | Per-sample `StreamingPipeline` (download→quant→cleanup) |
| [`orchestration.py`](orchestration.py) | Single-species workflow orchestration |
| [`orchestration_multi_species.py`](orchestration_multi_species.py) | Multi-species `PipelineOrchestrator` |
| [`workflow.py`](workflow.py) | Re-export hub for workflow_core, workflow_planning, workflow_execution |
| [`workflow_core.py`](workflow_core.py) | Config classes, validation, sample config creation |
| [`workflow_planning.py`](workflow_planning.py) | Workflow planning, step defaults, genome prep |
| [`workflow_execution.py`](workflow_execution.py) | Workflow execution engine, streaming mode |
| [`workflow_steps.py`](workflow_steps.py) | Individual step implementations |
| [`workflow_cleanup.py`](workflow_cleanup.py) | Disk cleanup, temp file removal |
| [`monitoring.py`](monitoring.py) | Real-time progress monitoring |
| [`discovery.py`](discovery.py) | Species and sample discovery via NCBI |
| [`pipeline.py`](pipeline.py) | Pipeline abstraction |
| [`progress_tracker.py`](progress_tracker.py) | Progress state persistence |
| [`sra_extraction.py`](sra_extraction.py) | SRA data extraction utilities |

## 🔑 Key Classes

### streaming_orchestrator.py

- `StreamingPipelineOrchestrator` - Multi-species end-to-end pipeline
  - `run_all()` - Process all species sequentially
  - `process_species()` - Download + quantify with ThreadPoolExecutor
  - `process_single_sample()` - Download → quant → cleanup per sample
  - ENA-first download with NCBI SRA fallback
  - Per-sample 2-hour timeout, automatic skip on failure
  - Downstream steps: merge → curate → sanity (30-min timeout)

### workflow.py

- `AmalgkitWorkflowConfig` - Load and validate YAML configs
- `execute_workflow()` - Main workflow entry point
- `_execute_streaming_mode()` - Disk-efficient chunked processing

### monitoring.py

- `WorkflowMonitor` - Real-time progress tracking
- `HeartbeatMonitor` - Process health checks

### progress_tracker.py

- `ProgressTracker` - Persistent progress state

## 🚀 Usage

```python
from metainformant.rna.engine.workflow import (
    AmalgkitWorkflowConfig,
    execute_workflow
)

config = AmalgkitWorkflowConfig.load("config/amalgkit/my_config.yaml")
result = execute_workflow(config)
```

## 🔗 Related

- [amalgkit/](../amalgkit/) - Tool wrappers
- [scripts/rna/](../../../../scripts/rna/) - CLI scripts
