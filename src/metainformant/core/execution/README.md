# Core Execution Module

Workflow orchestration, parallel processing, and script discovery for MetaInformAnt.

## Purpose

This module provides:
- Config-driven workflow execution
- Parallel processing with `ParallelProcessor`
- Dynamic script and config discovery

## Key Components

| File | Description |
|------|-------------|
| [workflow.py](workflow.py) | `run_config_based_workflow`, `validate_config_file` |
| [parallel.py](parallel.py) | `ParallelProcessor`, `run_parallel` |
| [discovery.py](discovery.py) | `discover_functions`, `discover_configs` |

## Usage

```python
from metainformant.core.execution import run_config_based_workflow

results = run_config_based_workflow("config.yaml")
```

## Related Documentation

- **Parent**: [src/metainformant/core/README.md](../README.md)
- **SPEC**: [SPEC.md](SPEC.md)
- **AGENTS**: [AGENTS.md](AGENTS.md)
- **Engine Module**: [../engine/README.md](../engine/README.md)
