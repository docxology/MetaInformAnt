# RNA Analysis Engine

## Overview

The **RNA Engine** is the central orchestration layer for the MetaInformAnt transcriptomics pipeline. It manages the execution of complex, multi-stage workflows (specifically Amalgkit) with a focus on **resilience**, **observability**, and **resource management**.

## Core Components

### 1. Workflow Orchestrator (`workflow.py`)

The intelligent heart of the engine.

* **Streaming Mode**: Automatically detects low-disk conditions and switches to chunked processing (e.g., download 5 -> quant 5 -> delete 5).
* **Resilience**: Implements fail-fast logic for disk space, "intelligent redo" to prevent data loss, and step-skipping for completed work.
* **Path Management**: Dynamically resolves paths for Amalgkit steps, handling custom output directories and cleanup.

### 2. Progress Tracking (`progress_tracker.py`)

maintains the state of the workflow.

* Tracks status of individual samples across steps.
* Persists state to prevent re-execution of successful steps.

### 3. Monitoring (`monitoring.py`)

Provides real-time visibility.

* **Heartbeats**: Logs progress and system vitals (disk space, memory) during long-running operations.
* **Process Watchers**: Monitors subprocess output for specific file patterns or directory growth.

### 4. Discovery (`discovery.py`)

Handles input identification.

* scans metadata and filesystems to identify valid samples for processing.

## Usage

The engine is typically invoked via the `run_workflow.py` script or directly through the `execute_workflow` function:

```python
from metainformant.rna.engine.workflow import load_workflow_config, execute_workflow

# Load validated configuration
config = load_workflow_config("config/my_species.yaml")

# Run robust, resumable workflow
execute_workflow(config)
```
