# Core Engine Module

The `metainformant.core.engine` module provides the high-level orchestration logic for multi-stage biological workflows.

## Overview

The primary component of this module is the `WorkflowManager`, which coordinates the execution of various analysis stages (download, extraction, quantification) while providing real-time progress monitoring.

## Key Components

### Workflow Manager
Coordinates samples and stages, ensuring robust execution with thread management and progress tracking.

## Usage Example

```python
from metainformant.core.engine.workflow_manager import WorkflowManager
from pathlib import Path

# Initialize manager
manager = WorkflowManager(config_path=Path("config.yaml"), max_threads=8)

# Add samples
manager.add_sample("sample1", "https://sra-download...", Path("data/sample1"))

# Run workflow
results = manager.run()
```
