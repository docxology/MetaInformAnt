# SPEC: RNA Scripts

Orchestration logic for RNA-seq analysis using Amalgkit and custom MetaInformAnt methods.

## Workflow Pipeline

1. **Environment Verification**: Check for R, Rscript, and Amalgkit dependencies.
2. **Metadata Acquisition**: Download and curate species-specific metadata.
3. **Download & Extraction**: Parallel retrieval and FASTQ extraction.
4. **Quantification**: Run Amalgkit quant with dynamic thread allocation.
5. **Post-Processing**: Merge, curate, and normalize results.

## Key Scripts

- `run_workflow.py`: The main entry point for end-to-end RNA-seq processing.
- `verify_rna.py`: Comprehensive diagnostic script for checking pipeline integrity.
- `run_workflow_tui.py`: TUI-enhanced version of the main workflow.

## Error Handling

- **Automatic Fallbacks**: If Amalgkit extraction fails, the script falls back to a custom SRA-to-FASTQ interceptor.
- **Heartbeat Monitoring**: Ensures long-running dump processes Haven't hung.
