# Engine

Generic pipeline orchestration engine that runs work items through ordered phases with parallel execution and TUI progress.

## Contents

| File | Purpose |
|------|---------|
| `workflow_manager.py` | Domain-agnostic multi-phase pipeline manager |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `Stage` | Lifecycle enum for pipeline items (pending, running, done, failed) |
| `PipelineItem` | One work item with its metadata, stage and error |
| `PipelinePhase` | A named phase: handler callable, item filter and TUI color |
| `BasePipelineManager` | Runs phases in order over registered items; subclass it for a domain pipeline |

## Usage

```python
from metainformant.core.engine.workflow_manager import BasePipelineManager, PipelinePhase


def process(manager, items):
    for item in items:
        manager.mark_running(item)
        manager.mark_done(item)


manager = BasePipelineManager(phases=[PipelinePhase(name="process", handler=process)], max_threads=4)
manager.add_item("sample-1", {"path": "data/sample-1.tsv"})
results = manager.run()  # {item_id: succeeded}
```
