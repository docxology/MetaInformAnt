# eQTL Configuration

Configuration files for expression Quantitative Trait Loci (eQTL) analysis.

## Files

| File | Description |
|------|-------------|
| `eqtl_amellifera.yaml` | *Apis mellifera* eQTL analysis config |

## Quick Start

```bash
# Verify the eQTL YAML loads with the shared config helper
uv run python - <<'PY'
from metainformant.core.utils.config import load_mapping_from_file

config = load_mapping_from_file("config/eqtl/eqtl_amellifera.yaml")
print(sorted(config))
PY
```

## Configuration Sections

- **expression**: RNA-seq data paths, normalization settings
- **variants**: VCF file, MAF filters
- **annotations**: Gene positions for cis-window
- **cis_eqtl**: Window size, FDR thresholds
- **trans_eqtl**: Trans analysis settings (optional)
- **colocalization**: GWAS integration
- **output**: Results and plots directories

## Data Sources

Requires:

1. **Expression data**: `output/amalgkit/apis_mellifera_all/work/quant/`
2. **Variants**: `output/gwas/amellifera/variants/`
3. **Annotations**: `output/gwas/amellifera/genome/genomic.gff`

## See Also

- [config/gwas/](../gwas/) - GWAS configuration
- [config/amalgkit/](../amalgkit/) - RNA-seq pipeline config
