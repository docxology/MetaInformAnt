# RNA Module Technical Specification

## Architecture Overview

The `rna` module in METAINFORMANT provides a robust Python interface to the **Amalgkit** transcriptomic analysis toolkit.

**Core Philosophy: Thin CLI Wrapper**
To ensure 100% feature parity and reliability, the module acts as a "thin wrapper" around the official `amalgkit` command-line interface (CLI). It does *not* attempt to reimplement `amalgkit`'s internal logic in native Python, which would risk drift and maintenance burden.

### Components

1.  **`amalgkit.py` (The Core)**:
    *   Contains the `AmalgkitParams` dataclass.
    *   Implements wrapper functions for all amalgkit subcommands (`metadata`, `getfastq`, `quant`, etc.).
    *   Handles command construction, execution (via `subprocess`), logging, and error propagation.
    *   Manages "heartbeat" monitoring for long-running processes (e.g., download/quantification progress).

2.  **`workflow.py` (Orchestration)**:
    *   Manages the execution graph of a multi-step analysis.
    *   Handles configuration loading (`AmalgkitWorkflowConfig`).
    *   Determines step dependencies and execution order (`plan_workflow`).
    *   Checks for completed steps to support resume/restart logic.
    *   Auto-configures environment (e.g., `vdb-config` for SRA tools).

3.  **`configs.py`**:
    *   Helpers for validating and applying default configurations.

4.  **`steps/` (Legacy/Specific Helpers)**:
    *   *Note: Historical native implementations in this directory are being deprecated in favor of direct CLI calls in `amalgkit.py`.*
    *   `process_samples.py` may remain if it provides unique parallelization logic *on top of* the CLI commands (e.g., concurrent process management).

## Data Flow

1.  **Input**: User provides a `config.yaml` or parameters via API.
2.  **Planning**: `workflow.py` translates config into a sequence of `amalgkit` commands.
3.  **Execution**: `amalgkit.py` executes each command as a subprocess.
    *   *Path Resolution*: `amalgkit.py` and `workflow.py` handle specific path logic (e.g., `getfastq` creates a subdir, `integrate` needs to point to it).
4.  **Monitoring**: Real-time log monitoring tracks progress (files downloaded, samples processed).
5.  **Output**: Standard Amalgkit directory structure (`metadata/`, `fastq/`, `quant/`, `merge/`).

## Testing Strategy

*   **No Mocks**: Tests must execute the real `amalgkit` CLI (or verify the command construction if the CLI isn't installed in the test env, but ideally real execution on test data).
*   **End-to-End**: `verify_workflow.sh` runs a complete pipeline on a small dataset.
*   **Unit Tests**: Verify that Python parameters are correctly translated to CLI arguments.

## Directory Structure

```
src/metainformant/rna/
├── __init__.py
├── amalgkit.py         # Main CLI wrapper implementation
├── workflow.py         # Orchestration logic
├── configs.py          # Configuration helpers
├── AGENTS.md           # AI collaboration log
├── README.md           # User guide
└── SPEC.md             # This file
```

## Amalgkit Version Compatibility
*   Target: Latest `amalgkit` release (via `uv pip install amalgkit`).
*   The wrapper should be agnostic to minor version changes by passing through valid CLI args.
