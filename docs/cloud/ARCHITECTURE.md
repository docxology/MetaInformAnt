# Cloud Module Architecture

## Overview

This page describes the historical GCP execution path. The Hymenoptera cloud
environment was decommissioned on 2026-04-15; local data-root verification is
the current evidence path.

The cloud deployment system uses a **thin orchestrator pattern**: a local CLI (`deploy_gcp.py`) shells out to `gcloud` CLI for all GCP operations. No Python SDK required.

## Module Design

```
src/metainformant/cloud/
 cloud_config.py # CloudConfig dataclass
 gcp_deployer.py # GCPDeployer class
```

### CloudConfig

Immutable configuration for a GCP VM:
- `project`: GCP project ID  
- `machine_type`: e.g. `n2-standard-16`
- `spot`: Use preemptible pricing (default: false, pipeline requires stability)
- `workers`/`threads`: Pipeline parallelism
- `gcs_bucket`: Optional result sync target

### GCPDeployer

Manages VM lifecycle via `subprocess` → `gcloud`:

| Method | Purpose |
|---|---|
| `create_vm()` | Provisions VM with startup script |
| `delete_vm()` | Tears down VM + disks |
| `get_vm_status()` | Checks instance state |
| `get_pipeline_status()` | SSH check of pipeline progress |
| `tail_logs()` | Remote log viewing |
| `download_results()` | SCP results to local machine |
| `full_deploy()` | End-to-end: create → wait for SSH → status |

## Pipeline on VM

The VM runs the same pipeline as locally:

1. **`prep_genomes.py`** — Downloads reference transcriptomes from NCBI FTP and builds kallisto indices
2. **`run_all_species.py`** — Launches `StreamingPipelineOrchestrator` for the configured 27-species cohort; the exact sample inventory is data-root dependent.
3. **`streaming_orchestrator.py`** — Phase 1 (task discovery/metadata) -> Phase 2 (bounded download → quant) -> Phase 3 (merge → filter → finalize)

### Container Configuration Requirements
The core pipeline executes within a Docker container (`metainformant:patched`). To ensure data persistence and prevent crash loops, the 4TB persistent data disk **must** be explicitly bind-mounted:
`-v /opt/MetaInformAnt/output:/app/output`
Starting the container with named volumes (e.g. `pipeline_data:/app/output`) will mask the underlying 4TB host drive and trigger a catastrophic re-download loop from NCBI.

## Data Flow

```
NCBI SRA download (ENA FTP/NCBI) FASTQ files

NCBI FTP transcriptome FASTA kallisto index

 FASTQ + index kallisto quant abundance.tsv

 merge merged_abundance.tsv

 wsfilter → finalize expression tables
```

## Result Download

Results are transferred via `gcloud compute scp`:
- `output/amalgkit/*/work/quant/` — Per-sample abundances
- `output/amalgkit/*/merged/` — Merged expression data
- `output/amalgkit/pipeline_progress.db` — SQLite status database
