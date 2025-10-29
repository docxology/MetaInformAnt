# AI Agents in RNA Script Development

This document records AI involvement in the high-throughput RNA processing scripts housed here.

## AI Contributions

### Pipeline Automation
**Code Assistant Agent** implemented:
- `batch_ena.py` for resilient ENA downloads, concurrent quantification, and log management
- Kallisto execution orchestration with dynamic thread allocation
- Resume and restart hooks consumed by `restart_batch.sh`

### Operational Guidance
**Documentation Agent** authored:
- Usage instructions mirrored in `README.md`
- Environment and dependency reminders for wget, kallisto, and uv workflows
- Integration notes tying script output paths to `docs/rna/examples`

### Validation Support
**Code Assistant Agent** cross-checked:
- CLI argument handling with amalgkit metadata expectations
- File placement within `output/` per repository rules
- Error handling branches used during large cohort reprocessing

## Maintenance Practices
- Capture real run logs when updating retry semantics or concurrency parameters
- Reflect new script flags in both this AGENTS file and the accompanying README
- Keep tests in `tests/rna` aligned with the automation behavior documented here

