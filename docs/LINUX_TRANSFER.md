# Resuming Amalgkit on Linux

> **Snapshot**: 2026-03-02 — Pipeline running on Linux with SQLite-backed progress tracking.

## Prerequisites

- Linux machine with Python 3.14+ and `uv` installed
- Git clone of MetaInformAnt

## Quick Start

```bash
# 1. Clone the repo
git clone https://github.com/docxology/MetaInformAnt.git
cd MetaInformAnt

# 2. Install dependencies
uv sync

# 3. Install the pinned Amalgkit release
uv sync --extra rna

# 4. Check pipeline status (instant, uses SQLite DB)
uv run python scripts/rna/check_pipeline_status.py -v
```

## Pipeline Progress Tracking

Progress is tracked via an **SQLite database** at `output/amalgkit/pipeline_progress.db`.
Status queries read this database instead of scanning the filesystem, which is too
slow for species like `amellifera` (7,300+ samples).

### Checking Status

```bash
# Quick overview
uv run python scripts/rna/check_pipeline_status.py

# Per-state breakdown (pending/downloading/quantifying/quantified/failed)
uv run python scripts/rna/check_pipeline_status.py -v

# Show failed samples with errors
uv run python scripts/rna/check_pipeline_status.py --failed

# Check a single species
uv run python scripts/rna/check_pipeline_status.py --species amellifera

# Generate visual dashboard (PDF + PNG)
uv run python scripts/rna/check_pipeline_status.py --dashboard
```

### Sample States

Each sample moves through the following states:

```
pending → downloading → downloaded → quantifying → quantified
            ↘ failed       ↘ failed        ↘ failed
```

### Key Modules

- `src/metainformant/rna/engine/progress_db.py` — `ProgressDB` class (SQLite, thread-safe)
- `src/metainformant/rna/engine/progress_dashboard.py` — Dashboard visualization (PDF/PNG)

## Resume Commands

### Run the full streaming pipeline

```bash
# Full run for the 27-species configuration (processes smallest first, streams download → quant → delete)
uv run python scripts/rna/run_all_species.py --max-gb 6.0 --workers 12 --threads 12
```

The orchestrator automatically:
- Initializes all samples in the SQLite DB as `pending`
- Reconciles existing quant results from the filesystem on startup
- Skips already-quantified samples
- Records state transitions atomically

### Resume for a single species

```bash
# The pipeline resumes where it left off — just re-run:
uv run python scripts/rna/run_all_species.py
```

## Config Files

All species configs are in `config/amalgkit/amalgkit_<species>.yaml`. Each config contains:

- `work_dir`: Path to working directory (e.g., `output/amalgkit/<species>/work`)
- `steps.quant.index_dir`: Path to Kallisto index
- `steps.quant.fasta_dir`: Path to genome FASTA
- `steps.getfastq.max_bp`: Max base pairs per sample (50M)

## Key Directories

```
output/amalgkit/
  pipeline_progress.db             # SQLite progress database (instant status queries)
  <species>/
    fastq/                         # Downloaded FASTQ files (deleted after quant)
    work/
      metadata/metadata.tsv        # Sample metadata (from SRA)
      index/                       # Kallisto index for this species
      quant/<SRR_ID>/              # Kallisto quantification results
      wsfilter/                    # Within-species filtered tables
      finalize/                    # Final expression tables
      sanity/                      # Sanity check outputs
```

## Known Issues

1. **Size filtering**: Samples exceeding `--max-gb` are skipped during processing
2. **FASTQ cleanup**: FASTQs are deleted after successful quantification to save disk space
3. **Failed samples**: Check with `--failed` flag; common causes are corrupt downloads or missing index files
