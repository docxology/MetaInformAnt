# PREPARATION

## Overview
Command-line helpers for GWAS Preparation workflows. Scripts should remain thin wrappers around `src/metainformant/` implementations and be run from the repository root with `uv`.

## Contents
- [download_genome_scale_data.sh](download_genome_scale_data.sh)
- [download_honeybee_variants.py](download_honeybee_variants.py)
- [generate_phenotypes.py](generate_phenotypes.py)
- [generate_synthetic_variants.py](generate_synthetic_variants.py)
- [query_bioproject_metadata.py](query_bioproject_metadata.py)

## Usage
```bash
uv run python scripts/gwas/preparation/download_honeybee_variants.py --help
```
