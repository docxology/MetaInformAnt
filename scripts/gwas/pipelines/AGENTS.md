# Agent Directives: scripts/gwas/pipelines

## Role
Complete GWAS pipeline execution scripts.

## Contents
- `run_analysis.py` - General GWAS analysis pipeline
- `run_genome_scale_gwas.py` - Genome-scale GWAS
- `run_pbarbatus_analysis.py` - P. barbatus specific pipeline
- `run_pbarbatus_gwas.py` - P. barbatus GWAS workflow

## Usage
These scripts run the complete GWAS workflow:
QC -> Structure -> Association -> Visualization

```bash
uv run python scripts/gwas/pipelines/run_analysis.py --config config/gwas/species.yaml
```
