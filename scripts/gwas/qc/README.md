# GWAS Quality Control

Scripts for filtering and validating variant data.

## Scripts

- `run_qc.py`: Standalone CLI for applying QC filters to VCF files.

## Usage

```bash
python scripts/gwas/qc/run_qc.py --vcf input.vcf --output filtered.vcf --min-maf 0.05
```
