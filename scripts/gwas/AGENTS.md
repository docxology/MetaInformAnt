# Agent Directives: scripts/gwas

## Role
GWAS pipeline workflow scripts organized by analysis stage.

## Directory Structure
- `association/` - Association testing scripts
- `pipelines/` - Complete GWAS pipeline scripts
- `preparation/` - Data preparation and download
- `qc/` - Quality control scripts
- `structure/` - Population structure analysis
- `visualization/` - GWAS visualization scripts

## Key Workflows
1. Preparation: Download data, generate phenotypes
2. QC: Run quality control filters
3. Structure: PCA, kinship matrix calculation
4. Association: Run association tests
5. Visualization: Generate Manhattan, QQ plots

## Usage
```bash
# Full pipeline
uv run python scripts/gwas/pipelines/run_analysis.py --config config/gwas/species.yaml

# Individual steps
uv run python scripts/gwas/qc/run_qc.py
uv run python scripts/gwas/association/run_association.py
```
