# SPEC: Quality Scripts

Orchestration of sequencing quality assessment and contamination detection.

## Workflows

- `run_qc_pipeline.py`: Comprehensive analysis of FASTQ files including adapter detection and base quality profiling.
- `verify_contamination.py`: Searches for non-target species sequences in genomic data.

## Standards

- **Reporting**: QC reports are generated in HTML and JSON formats for both human review and automated parsing.
- **Integration**: Leverage the `metainformant.quality` core module for analysis logic.
