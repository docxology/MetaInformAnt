# GWAS

## Overview
Command-line helpers for GWAS workflows. Scripts should remain thin wrappers around `src/metainformant/` implementations and be run from the repository root with `uv`.

## Contents
- [association/](association/)
- [pipelines/](pipelines/)
- [preparation/](preparation/)
- [qc/](qc/)
- [structure/](structure/)
- [visualization/](visualization/)
- [run_amellifera_gwas.py](run_amellifera_gwas.py)
- [run_pbarbatus_gwas.py](run_pbarbatus_gwas.py)
- [validate_all_methods.py](validate_all_methods.py)

## Usage
```bash
uv run python scripts/gwas/run_amellifera_gwas.py --help
```
