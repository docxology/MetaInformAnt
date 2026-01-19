# SPEC: Core Engine

The core engine is designed to be a domain-agnostic workflow orchestrator for biological data processing.

## Design Goals

1. **Robustness**: Handle network failures and tool-specific errors gracefully.
2. **Observability**: Provide real-time feedback via terminal UI.
3. **Concurrency**: Efficiently manage multi-threaded data acquisition and processing.

## State Management

Workflow stages are tracked as `SampleStage` objects with associated `SampleState` enum values:
- `PENDING`: Waiting for execution.
- `RUNNING`: Currently active.
- `COMPLETED`: Successfully finished.
- `FAILED`: Encountered an error.

## Flow Control

The engine manages a pool of workers and schedules tasks based on dependency completion. Currently, the flow is linear (Download -> Extract -> Quant), but the architecture supports DAG-based workflows.
