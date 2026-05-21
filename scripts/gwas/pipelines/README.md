# PIPELINES

## Overview
Command-line helpers for GWAS Pipelines workflows. Scripts should remain thin wrappers around `src/metainformant/` implementations and be run from the repository root with `uv`.

## Contents
- [run_analysis.py](run_analysis.py)
- [run_genome_scale_gwas.py](run_genome_scale_gwas.py)
- [run_pbarbatus_analysis.py](run_pbarbatus_analysis.py)
- [run_pbarbatus_gwas.py](run_pbarbatus_gwas.py)

## Usage
```bash
uv run python scripts/gwas/pipelines/run_analysis.py --help
```
