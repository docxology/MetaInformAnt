# RNA Scripts

Scripts for RNA-seq workflow execution, monitoring, and recovery.

## 📦 Scripts by Category

### Workflow Execution

| Script | Purpose |
|--------|---------|
| [`run_workflow.py`](run_workflow.py) | Main workflow orchestrator - runs full amalgkit pipeline |
| [`run_workflow_tui.py`](run_workflow_tui.py) | Terminal UI for interactive workflow control |
| [`process_apis_mellifera.py`](process_apis_mellifera.py) | Sequential processor for A. mellifera samples |
| [`process_apis_mellifera_parallel.py`](process_apis_mellifera_parallel.py) | **NCBI parallel processor** (10 workers) with 4-hour timeout |
| [`process_apis_mellifera_ena.py`](process_apis_mellifera_ena.py) | **ENA parallel processor** (20 workers) - 10x faster, skips extraction |
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

### Run A. mellifera Processing (Dual-Pipeline Strategy)

**Recommended**: Use ENA pipeline for 10x faster processing, then NCBI fallback for failures.

```bash
# Step 1: Run ENA pipeline (20 workers, pre-extracted FASTQs)
nohup python scripts/rna/process_apis_mellifera_ena.py &

# Monitor ENA progress
tail -f output/amalgkit/apis_mellifera_all/work/ena_parallel_processing.log

# Step 2: After ENA completes, run NCBI fallback for failures
nohup python scripts/rna/process_apis_mellifera_parallel.py &
```

**Performance Comparison:**

| Pipeline | Workers | Time/Sample | ETA (5000 samples) |
|----------|---------|-------------|-------------------|
| ENA | 20 | 0.5-5 min | ~12-24 hours |
| NCBI | 10 | 45-60 min | ~7 days |

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
