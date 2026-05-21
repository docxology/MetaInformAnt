# GWAS

## Overview
Configuration assets for GWAS workflows. Keep templates and examples synchronized with the shared loaders in `metainformant.core.utils.config`.

## Contents
- [gwas_amellifera.yaml](gwas_amellifera.yaml)
- [gwas_pbarbatus.yaml](gwas_pbarbatus.yaml)
- [gwas_template.yaml](gwas_template.yaml)

## Usage
```python
from metainformant.core.utils.config import load_mapping_from_file

config = load_mapping_from_file("config/gwas/gwas_amellifera.yaml")
```
