# RNA Scripts

Scripts for RNA-seq workflow execution, monitoring, and recovery.

## 📦 Scripts by Category

### Workflow Execution

| Script | Purpose |
|--------|---------|
| [`run_workflow.py`](run_workflow.py) | Main workflow orchestrator - runs full amalgkit pipeline |
| [`run_workflow_tui.py`](run_workflow_tui.py) | Terminal UI for interactive workflow control |
| [`process_apis_mellifera.py`](process_apis_mellifera.py) | Sequential processor for A. mellifera samples |
| [`process_apis_mellifera_parallel.py`](process_apis_mellifera_parallel.py) | **Parallel processor** (4 workers) for A. mellifera |
| [`process_samples_sequential.py`](process_samples_sequential.py) | Generic sequential sample processor |

### Monitoring

| Script | Purpose |
|--------|---------|
| [`monitor_apis_mellifera.py`](monitor_apis_mellifera.py) | Real-time progress dashboard with `--watch` mode |
| [`verify_rna.py`](verify_rna.py) | Comprehensive RNA module verification |
| [`check_environment.py`](check_environment.py) | Validate dependencies (kallisto, prefetch, etc.) |

### Recovery & Utilities

| Script | Purpose |
|--------|---------|
| [`recover_missing_batch.py`](recover_missing_batch.py) | Retry failed samples via ENA download |
| [`recover_missing_parallel.py`](recover_missing_parallel.py) | Parallel retry with multiple workers |
| [`filter_valid_samples.py`](filter_valid_samples.py) | Filter metadata to valid RNA-seq samples |
| [`setup_genome.py`](setup_genome.py) | Download and prepare reference genome |

### Discovery & Validation

| Script | Purpose |
|--------|---------|
| [`discover_species.py`](discover_species.py) | Query NCBI for species RNA-seq data availability |
| [`validate_all_species_workflow.py`](validate_all_species_workflow.py) | End-to-end workflow validation |
| [`verify_fallback.py`](verify_fallback.py) | Test fallback download mechanisms |
| [`verify_tui.py`](verify_tui.py) | TUI component verification |

### Setup

| Script | Purpose |
|--------|---------|
| [`_setup_utils.py`](_setup_utils.py) | Internal setup utilities |
| [`install_r_deps.R`](install_r_deps.R) | Install R dependencies for curation |
| [`install_r_packages.sh`](install_r_packages.sh) | Shell wrapper for R package installation |
| [`run_rna_tests.sh`](run_rna_tests.sh) | Run RNA test suite |

## 🚀 Quick Start

### Run A. mellifera Parallel Processing

```bash
# Start parallel processor (4 workers)
nohup uv run python scripts/rna/process_apis_mellifera_parallel.py &

# Monitor progress
uv run python scripts/rna/monitor_apis_mellifera.py --watch
```

### Check Environment

```bash
uv run python scripts/rna/check_environment.py
```

### Discover Species Data

```bash
uv run python scripts/rna/discover_species.py --species "Apis mellifera"
```

## 📊 Dependencies

- **Python**: 3.12+
- **External**: kallisto, prefetch, fasterq-dump (SRA Toolkit)
- **Optional**: R (for merge/curate steps)

## 🔗 Related

- [config/amalgkit/](../../config/amalgkit/) - Workflow configurations
- [config/amalgkit/amalgkit_faq.md](../../config/amalgkit/amalgkit_faq.md) - FAQ
- [src/metainformant/rna/](../../src/metainformant/rna/) - Core RNA module
