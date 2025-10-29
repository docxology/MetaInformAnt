# AI Agents in RNA Step Development

This document outlines AI assistance across the modular RNA workflow steps.

## AI Contributions

### Step Implementations
**Code Assistant Agent** developed:
- `metadata.py`, `integration.py`, and related utilities for metadata acquisition and harmonization
- `getfastq.py`, `quant.py`, and `merge.py` for download, quantification, and aggregation stages
- `cstmm.py`, `csca.py`, and `sanity.py` for statistical validation and consistency checks

### Orchestration Patterns
**Code Assistant Agent** aligned:
- Step signatures with `AmalgkitWorkflowConfig` planning and execution
- Structured `StepResult` outputs for downstream chaining
- Error handling and logging messages consistent with `metainformant.core.logging`

### Documentation Support
**Documentation Agent** contributed to:
- The accompanying README descriptions and integration examples
- Cross-references to `docs/rna/steps.md` and workflow tutorials
- Maintenance notes for adding new steps without breaking orchestration

## Maintenance Guidelines
- Ensure any new step exposes `run_step()` with standardized inputs/outputs
- Update integration tests in `tests/rna` when modifying shared utilities
- Flag AI-assisted edits for human review before merging into mainline branches

