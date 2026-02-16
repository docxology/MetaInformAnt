# RNA Scripts

Scripts for RNA-seq workflow execution, monitoring, and recovery.

## 📦 Scripts by Category

### Workflow Execution

| Script | Purpose |
|--------|---------|
| [`run_all_species.py`](run_all_species.py) | **Primary multi-species orchestrator** — streaming ENA-first pipeline |
| [`run_workflow.py`](run_workflow.py) | Single-species workflow (merge, curate, sanity steps) |
| [`run_workflow_tui.py`](run_workflow_tui.py) | Terminal UI for interactive workflow control |
| [`process_apis_mellifera.py`](process_apis_mellifera.py) | Sequential processor for A. mellifera samples |
| [`process_apis_mellifera_parallel.py`](process_apis_mellifera_parallel.py) | NCBI parallel processor (10 workers) with 4-hour timeout |
| [`process_apis_mellifera_ena.py`](process_apis_mellifera_ena.py) | ENA parallel processor (20 workers) |
| [`process_samples_sequential.py`](process_samples_sequential.py) | Generic sequential sample processor |

### Monitoring

| Script | Purpose |
|--------|---------|
| [`check_pipeline_status.py`](check_pipeline_status.py) | **Pipeline progress dashboard** — per-species quant completion |
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
| [`sync_quant_results.py`](sync_quant_results.py) | **Sync quant results** to git-trackable `output/amalgkit_results/` |

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

### Run Multi-Species Pipeline (Recommended)

The **streaming orchestrator** processes all species end-to-end: download → quant → merge → curate → sanity.

```bash
# Run all species (8 workers, ENA-first with NCBI fallback)
nohup uv run --python 3.11 python scripts/rna/run_all_species.py \
  --workers 8 --threads 8 --max-gb 100 \
  > output/amalgkit/pipeline.log 2>&1 &

# Monitor progress
tail -f output/amalgkit/pipeline.log

# Check per-species completion status
uv run python scripts/rna/check_pipeline_status.py

# Sync results to git-trackable location
uv run python scripts/rna/sync_quant_results.py
```

### A. mellifera (Large Dataset)

For the 7,000+ sample A. mellifera dataset, use the dedicated ENA pipeline:

```bash
nohup python scripts/rna/process_apis_mellifera_ena.py &
tail -f output/amalgkit/apis_mellifera_all/work/ena_parallel_processing.log
```

| Pipeline | Workers | Time/Sample | Notes |
|----------|---------|-------------|-------|
| Streaming (all species) | 8-16 | 2-10 min | ENA-first, auto-fallback |
| ENA (A. mellifera) | 20 | 0.5-5 min | Pre-extracted FASTQs |
| NCBI fallback | 10 | 45-60 min | SRA Toolkit required |

### Check Environment

```bash
uv run python scripts/rna/check_environment.py
```

### Discover Species Data

```bash
uv run python scripts/rna/discover_species.py --species "Apis mellifera"
```

## 📊 Dependencies

- **Python**: 3.11+
- **External**: kallisto, prefetch, fasterq-dump (SRA Toolkit)
- **Optional**: R (for merge/curate steps)

## 🔗 Related

- [config/amalgkit/](../../config/amalgkit/) - Workflow configurations
- [config/amalgkit/amalgkit_faq.md](../../config/amalgkit/amalgkit_faq.md) - FAQ
- [src/metainformant/rna/](../../src/metainformant/rna/) - Core RNA module
