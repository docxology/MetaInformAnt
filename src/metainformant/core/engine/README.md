# Engine

Pipeline orchestration engine that manages multi-phase workflows with sample-level state tracking, parallel execution, and progress monitoring.

## Contents

| File | Purpose |
|------|---------|
| `workflow_manager.py` | Pipeline manager with download, getfastq, and quantification phases |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `Stage` | Enum of pipeline stages (e.g., download, processing, complete) |
| `PipelineItem` | Tracks an individual sample through the pipeline |
| `PipelinePhase` | Defines a named phase with its execution function |
| `BasePipelineManager` | Abstract base for pipeline managers with phase registration |
| `WorkflowManager` | Concrete manager for RNA-seq download/quant workflows |
| `SampleStage` | Fine-grained sample lifecycle enum (queued, downloading, etc.) |
| `SampleState` | Mutable state object per sample with timestamps and error info |

## Usage

```python
from metainformant.core.engine.workflow_manager import WorkflowManager, PipelineItem

manager = WorkflowManager(config=workflow_config, output_dir="output/run1")
items = [PipelineItem(sample_id=sid, metadata=meta) for sid, meta in samples.items()]
manager.run(items)
```
