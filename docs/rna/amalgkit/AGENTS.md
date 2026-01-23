# Agent Directives: docs/rna/amalgkit

## Role
Documentation for amalgkit RNA-seq workflow integration.

## Directory Structure
- `steps/` - Step-by-step workflow documentation (01-11)

## Key Files
- `amalgkit.md` - Complete amalgkit integration guide
- `FUNCTIONS.md` - API function reference
- `commands.md` - CLI command reference
- `genome_preparation.md` - Genome setup for quantification
- `monitoring.md` - Workflow monitoring
- `R_INSTALLATION.md` - R dependency setup

## Workflow Steps (in steps/)
1. metadata - Fetch sample metadata from SRA
2. config - Generate workflow configuration
3. select - Select samples for processing
4. getfastq - Download FASTQ files
5. integrate - Prepare genome indices
6. quant - Quantify expression with Salmon
7. merge - Merge sample quantifications
8. cstmm - Cross-species TMM normalization
9. curate - Curate and filter results
10. csca - Cross-species comparative analysis
11. sanity - Validation checks
