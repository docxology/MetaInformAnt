# Configuration Directory

Configuration files for METAINFORMANT workflows including RNA-seq, GWAS, and multi-omic analysis.

## Directory Structure

```
config/
├── amalgkit/           # RNA-seq workflow configurations
│   ├── amalgkit_pbarbatus_all.yaml    # Pogonomyrmex barbatus (active)
│   ├── amalgkit_template.yaml         # Template for new species
│   └── amalgkit_test.yaml             # Test configuration
├── gwas/               # GWAS workflow configurations
│   ├── gwas_pbarbatus.yaml            # P. barbatus GWAS config
│   └── gwas_template.yaml             # Template for new GWAS
├── archive/            # Inactive/deprecated configurations
├── config_base/        # Base configuration templates
├── phenotype/          # Phenotype analysis configurations
├── ncbi.yaml           # NCBI API settings
├── life_events_template.yaml
├── multiomics_template.yaml
├── networks_template.yaml
└── singlecell_template.yaml
```

## Quick Start

### Run RNA-seq Workflow

```bash
# Using CLI
uv run metainformant rna run-config --config config/amalgkit/amalgkit_pbarbatus_all.yaml

# Using script
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus_all.yaml
```

### Run GWAS Workflow

```bash
uv run metainformant gwas run --config config/gwas/gwas_pbarbatus.yaml
```

## Configuration Format

### Amalgkit RNA-seq Configuration

```yaml
work_dir: output/amalgkit/species/work
log_dir: output/amalgkit/species/logs
threads: 8
species_list: ["Species_name"]

genome:
  accession: "GCF_XXXXXXXXX.X"
  assembly_name: "Assembly_Name"
  dest_dir: "output/amalgkit/species/genome"
  include: ["genome", "gff3", "rna", "cds", "protein"]

# Workflow steps (each can be enabled/disabled)
metadata:
  enabled: true
  max_samples: 50
integrate:
  enabled: true
getfastq:
  enabled: true
  prefetch: true
quant:
  enabled: true
  index: "kallisto"
merge:
  enabled: true
curate:
  enabled: true
```

### GWAS Configuration

```yaml
work_dir: output/gwas/species
threads: 8

genome:
  accession: "GCF_XXXXXXXXX.X"
  dest_dir: "output/gwas/species/genome"

variants:
  vcf_path: "data/variants/cohort.vcf.gz"

samples:
  phenotype_file: "data/phenotypes/traits.tsv"
  trait_column: "phenotype"

qc:
  min_maf: 0.01
  max_missing: 0.1
  min_hwe_p: 1e-6

association:
  model: "linear"
  covariates: ["pc1", "pc2", "pc3"]
```

## Environment Variable Overrides

All configuration parameters can be overridden via environment variables:

| Prefix | Module | Example |
|--------|--------|---------|
| `AK_` | RNA/Amalgkit | `AK_THREADS=16` |
| `GWAS_` | GWAS | `GWAS_WORK_DIR=/scratch/gwas` |
| `DNA_` | DNA | `DNA_WORK_DIR=output/dna` |
| `NCBI_` | NCBI API | `NCBI_EMAIL=user@example.com` |

## Loading Configurations

```python
from metainformant.core.config import load_mapping_from_file, apply_env_overrides

# Load with environment overrides
config = load_mapping_from_file("config/amalgkit/amalgkit_pbarbatus_all.yaml")
config = apply_env_overrides(config, prefix="AK")
```

## Active Species

| Species | Config | NCBI Assembly |
|---------|--------|---------------|
| *Pogonomyrmex barbatus* | `amalgkit_pbarbatus_all.yaml` | GCF_000187915.1 |

See `config/archive/README.md` for archived species configurations.

## NCBI Configuration

The `ncbi.yaml` file contains API settings:

```yaml
email: "your.email@example.com"  # Required for NCBI API
rate_limit_delay: 0.34           # Seconds between requests
max_retries: 3
```

Override with `NCBI_EMAIL` environment variable.

## Related Documentation

- [RNA Workflow Guide](../docs/rna/workflow.md)
- [GWAS Workflow Guide](../docs/gwas/workflow.md)
- [Core Config Module](../src/metainformant/core/config.py)
