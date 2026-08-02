# Agent Directives: docs/rna/amalgkit/steps

## Role
Step-by-step documentation for each amalgkit workflow stage.

## Contents
Numbered documentation files (01-11) corresponding to workflow order:
- `01_metadata.md` - SRA metadata retrieval
- `02_dataset.md` - Workspace and rule-set initialization
- `03_select.md` - Sample selection criteria
- `04_getfastq.md` - FASTQ download process
- `05_integrate.md` - Local-read integration
- `06_quant.md` - Kallisto quantification
- `07_merge.md` - Sample merging
- `09_wsfilter.md` - Within-species filtering
- `10_finalize.md` - Final analysis tables
- `10_csfilter.md` - Optional cross-species filtering
- `11_sanity.md` - Validation checks

`cstmm`, `busco`, and `rerun` are current auxiliary commands documented in
the command reference rather than default numbered stages.

## Usage
Read in order for complete workflow understanding. Each file documents:
- Input requirements
- Output files
- Configuration options
- Common issues and solutions
