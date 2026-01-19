# AI Agents: Core Engine

This document provides function indices and technical details for the `core.engine` module, intended for use by AI agents.

## Function Index

### `metainformant.core.engine.workflow_manager`

- **`WorkflowManager(config_path: Path, max_threads: int = 5)`**: Main entry point for workflow orchestration.
- **`WorkflowManager.add_sample(sample_id: str, sra_url: str, dest_path: Path) -> None`**: Register a sample for processing.
- **`WorkflowManager.run() -> dict[str, bool]`**: Execute the planned workflow stages.

## Implementation Details

The `WorkflowManager` implements a stage-based execution model:
1. **Download Phase**: Concurrent retrieval of raw data.
2. **Extraction Phase**: Standardizing data formats (e.g., SRA to FASTQ).
3. **Quantification Phase**: Domain-specific analysis (e.g., RNA-seq quantification).

It uses `metainformant.core.progress` for UI updates and handles error recovery at each stage.
