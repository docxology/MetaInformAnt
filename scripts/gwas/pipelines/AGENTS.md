# Agent Directives: scripts/gwas/pipelines

## Role
Thin CLI wrappers for complete GWAS pipeline execution.

## Contents
- `run_analysis.py` - General GWAS analysis pipeline (thin wrapper over `metainformant.gwas.workflow.pbarbatus_full_scale`)
- `run_genome_scale_gwas.py` - Genome-scale GWAS
- `run_pbarbatus_analysis.py` - P. barbatus specific pipeline (thin wrapper over `metainformant.gwas.workflow.pbarbatus_comprehensive`)
- `run_pbarbatus_gwas.py` - P. barbatus GWAS workflow (thin wrapper over `metainformant.gwas.workflow.pbarbatus_end_to_end`)

## Usage
These scripts run the complete GWAS workflow:
QC -> Structure -> Association -> Visualization

```bash
uv run python scripts/gwas/pipelines/run_analysis.py --config config/gwas/species.yaml
```
