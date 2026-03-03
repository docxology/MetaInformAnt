# Cloud Module Architecture

## Overview

The cloud deployment system uses a **thin orchestrator pattern**: a local CLI (`deploy_gcp.py`) shells out to `gcloud` CLI for all GCP operations. No Python SDK required.

## Module Design

```
src/metainformant/cloud/
├── cloud_config.py    # CloudConfig dataclass
└── gcp_deployer.py    # GCPDeployer class
```

### CloudConfig

Immutable configuration for a GCP VM:
- `project`: GCP project ID  
- `machine_type`: e.g. `n2-highcpu-32`
- `spot`: Use preemptible pricing (default: true)
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
2. **`run_all_species.py`** — Launches `StreamingPipelineOrchestrator` for all species
3. **`streaming_orchestrator.py`** — Per-species: metadata → filter → download → quant → merge → curate

## Data Flow

```
NCBI SRA ──▸ download (ENA FTP/NCBI) ──▸ FASTQ files
                                              │
NCBI FTP ──▸ transcriptome FASTA ──▸ kallisto index
                                              │
                              FASTQ + index ──▸ kallisto quant ──▸ abundance.tsv
                                                                       │
                                              merge ──▸ merged_abundance.tsv
                                                                       │
                                              curate ──▸ expression tables
```

## Result Download

Results are transferred via `gcloud compute scp`:
- `output/amalgkit/*/work/quant/` — Per-sample abundances
- `output/amalgkit/*/merged/` — Merged expression data
- `output/amalgkit/pipeline_progress.db` — SQLite status database
