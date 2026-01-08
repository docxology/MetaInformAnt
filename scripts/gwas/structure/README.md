# Population Structure Analysis

Scripts for analyzing and accounting for population structure.

## Scripts

- `run_pca.py`: Calculate Principal Components Analysis (PCA).
- `run_kinship.py`: Calculate Kinship/Relatedness matrices.

## Usage

```bash
python scripts/gwas/structure/run_pca.py --vcf filtered.vcf --output pca_results.json
python scripts/gwas/structure/run_kinship.py --vcf filtered.vcf --output kinship.json
```
