# Resuming Amalgkit on Linux

> **Snapshot**: 2026-03-02 — Pipeline running on Linux with SQLite-backed progress tracking.

## Prerequisites

- Linux machine with Python 3.11+ and `uv` installed
- Git clone of MetaInformAnt

## Quick Start

```bash
# 1. Clone the repo
git clone https://github.com/docxology/MetaInformAnt.git
cd MetaInformAnt

# 2. Install dependencies
uv sync

# 3. Install amalgkit
uv pip install amalgkit

# 4. Check pipeline status (instant, uses SQLite DB)
uv run python scripts/rna/check_pipeline_status.py -v
```

## Pipeline Progress Tracking

Progress is tracked via an **SQLite database** at `output/amalgkit/pipeline_progress.db`.
This replaces the old approach of scanning the filesystem for files, which was too slow
for species like `amellifera` (7,300+ samples).

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
- `src/metainformant/rna/engine/progress_tracker.py` — **Deprecated** (JSON-based, kept for backward compat)

## Resume Commands

### Run the full streaming pipeline

```bash
# Run all 28 species (processes smallest first, streams download → quant → delete)
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
      curate/                      # Curated expression tables
      sanity/                      # Sanity check outputs
```

## Known Issues

1. **Size filtering**: Samples exceeding `--max-gb` are skipped during processing
2. **FASTQ cleanup**: FASTQs are deleted after successful quantification to save disk space
3. **Failed samples**: Check with `--failed` flag; common causes are corrupt downloads or missing index files

## Cloud Deployment (GCP)

For large-scale processing, deploy the pipeline to a high-core GCP VM. We use a standardized approach that separates the data disk from the compute instance, allowing you to use STANDARD provisioning to prevent unexpected preemptions.

### 1. Create a Persistent 4TB Disk
Store all your pipeline data on a standalone disk that will not be deleted if the VM is destroyed:
```bash
gcloud compute disks create metainformant-pipeline-data \
    --size=4096GB \
    --type=pd-standard \
    --zone=us-central1-a \
    --project=YOUR_PROJECT_ID
```

### 2. Launch a STANDARD VM (n2-standard-16)
Attach the existing disk to a new STANDARD instance (`n2-standard-16` provides 16 vCPUs and 64GB RAM):
```bash
gcloud compute instances create metainformant-pipeline \
    --zone=us-central1-a \
    --project=YOUR_PROJECT_ID \
    --machine-type=n2-standard-16 \
    --provisioning-model=STANDARD \
    --maintenance-policy=MIGRATE \
    --restart-on-failure \
    --disk=name=metainformant-pipeline-data,device-name=persistent-disk-0,mode=rw,boot=yes,auto-delete=no \
    --scopes=default
```
*(Note: SPOT instances are significantly cheaper but will be terminated within 24 hours. STANDARD instances guarantee completion without data loss or pipeline interruption).*

### 3. Launch the Pipeline Container
SSH into the instance and run the Docker container. Tune `PIPELINE_WORKERS` to ~1.5x your vCPU count (e.g., 24 workers for 16 cores):
```bash
sudo docker run -d \
    --name metainformant-pipeline-fresh \
    --restart unless-stopped \
    -v /opt/MetaInformAnt/output:/app/output \
    -v /opt/MetaInformAnt/config/amalgkit:/app/config/amalgkit \
    -e PIPELINE_WORKERS=24 \
    -e PIPELINE_THREADS=2 \
    -e PIPELINE_MAX_GB=50.0 \
    docxology/metainformant-pipeline:latest
```

### 4. Monitor & Teardown
- **Monitor:** `sudo docker logs --tail 100 -f metainformant-pipeline-fresh`
- **Teardown Compute:** `gcloud compute instances delete metainformant-pipeline` (leaves data disk intact)
- **Teardown Disk:** `gcloud compute disks delete metainformant-pipeline-data` (once results are downloaded)
