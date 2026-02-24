# RNA Scripts

Scripts for RNA-seq workflow execution, monitoring, and recovery.

## 📦 Scripts by Category

### Workflow Execution

| Script | Purpose |
|--------|---------|
| [`run_all_species.sh`](run_all_species.sh) | **Primary pipeline** — runs all species sequentially with per-sample concurrency |
| [`run_workflow.py`](run_workflow.py) | Single-species workflow orchestrator (all steps including merge, curate, sanity) |
| [`run_workflow_tui.py`](run_workflow_tui.py) | Terminal UI for interactive workflow control |

### Monitoring

| Script | Purpose |
|--------|---------|
| [`check_pipeline_status.py`](check_pipeline_status.py) | **Pipeline progress dashboard** — per-species quant completion |
| [`verify_rna.py`](verify_rna.py) | Comprehensive RNA module verification |
| [`check_environment.py`](check_environment.py) | Validate dependencies (kallisto, prefetch, etc.) |

### Recovery & Utilities

| Script | Purpose |
|--------|---------|
| [`filter_valid_samples.py`](filter_valid_samples.py) | Filter metadata to valid RNA-seq samples |
| [`setup_genome.py`](setup_genome.py) | Download and prepare reference genome |
| [`sync_quant_results.py`](sync_quant_results.py) | **Sync quant results** to git-trackable `output/amalgkit_results/` |

### Discovery & Validation

| Script | Purpose |
|--------|---------|
| [`discover_species.py`](discover_species.py) | Query NCBI for species RNA-seq data availability |
| [`validate_all_species_workflow.py`](validate_all_species_workflow.py) | End-to-end workflow validation |
| [`validate_configs.py`](validate_configs.py) | Validate all species config YAML files |

### Setup

| Script | Purpose |
|--------|---------|
| [`_setup_utils.py`](_setup_utils.py) | Internal setup utilities |
| [`install_r_deps.R`](install_r_deps.R) | Install R dependencies for curation |
| [`install_r_packages.sh`](install_r_packages.sh) | Shell wrapper for R package installation |
| [`run_rna_tests.sh`](run_rna_tests.sh) | Run RNA test suite |

## 🚀 Quick Start

### Run Full Pipeline (Recommended)

The pipeline processes all 23 species sequentially. Within each species, samples are processed concurrently (chunk-size 6) with per-sample `getfastq → quant → cleanup`.

```bash
# Run all species (background)
nohup bash scripts/rna/run_all_species.sh > output/amalgkit/run_all_species_incremental.log 2>&1 &

# Monitor progress
.venv/bin/python scripts/package/generate_custom_summary.py

# Check active downloads
ps -fC amalgkit | grep SRR

# Check per-species completion status
.venv/bin/python scripts/rna/check_pipeline_status.py
```

### Run Single Species

```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_acromyrmex_echinatior.yaml --stream --chunk-size 6
```

### Check Environment

```bash
python3 scripts/rna/check_environment.py
```

### Discover Species Data

```bash
python3 scripts/rna/discover_species.py --species "Apis mellifera"
```

## 📊 Dependencies

- **Python**: 3.11+
- **External**: kallisto, fasterq-dump (SRA Toolkit)
- **Optional**: R (for merge/curate steps)

## 🔗 Related

- [config/amalgkit/](../../config/amalgkit/) - Workflow configurations
- [config/amalgkit/amalgkit_faq.md](../../config/amalgkit/amalgkit_faq.md) - FAQ
- [src/metainformant/rna/](../../src/metainformant/rna/) - Core RNA module
