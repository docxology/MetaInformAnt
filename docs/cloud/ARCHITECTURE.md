# Cloud Module Architecture

## Overview

This page describes the historical GCP execution path. The Hymenoptera cloud
environment was decommissioned on 2026-04-15; local data-root verification is
the current evidence path.

The cloud deployment system uses a **thin orchestrator pattern**: a local CLI (`deploy_gcp.py`) shells out to the `gcloud` CLI for all GCP operations. No Python SDK is required. On the VM, the startup script uses `gsutil` for the optional bucket sync.

## Module Design

```text
src/metainformant/cloud/
  cloud_config.py   # CloudConfig dataclass
  gcp_deployer.py   # GCPDeployer class
```

### CloudConfig

A mutable dataclass describing one GCP VM and its pipeline settings (`deploy_gcp.py` assigns fields after construction):
- `project`: GCP project ID
- `machine_type`: e.g. `n2-standard-16`
- `spot`: Use Spot provisioning (default: true; `spot=False` or `--no-spot` selects on-demand pricing)
- `workers`/`threads`: Pipeline parallelism
- `gcs_bucket`: Optional result sync target

### GCPDeployer

Manages VM lifecycle via `subprocess` → `gcloud`:

| Method | Purpose |
|---|---|
| `create_vm()` | Provisions VM with startup script |
| `delete_vm()` | Tears down VM + disks |
| `stop_vm()` / `start_vm()` | Stops the VM (keeping its disk) / starts it again |
| `get_vm_status()` | Checks instance state |
| `get_pipeline_status()` | SSH check of pipeline progress |
| `tail_logs()` | `docker logs` of the pipeline container over SSH |
| `get_startup_log()` | Startup script log over SSH |
| `download_results()` | Runs `scripts/cloud/download_results.sh` (copy on the VM, then SCP) |
| `full_deploy()` | End-to-end: create → wait for SSH → status |

## Pipeline on VM

The startup script builds the repository `Dockerfile` and starts one container, `metainformant-pipeline`. The container runs the same producer as a local run:

1. **`process_species.py`** builds a `StreamingPipelineOrchestrator` and calls `run_all` for one species (`--species apis_mellifera` in the startup script).
2. **Phase 1, species discovery** runs in parallel: when the kallisto index is missing it first downloads the reference transcriptome and runs `kallisto index` (`prepare_reference_genome`), then it prepares the metadata and runs selection.
3. **Phase 2, bounded sample execution**: ENA-first FASTQ download with an NCBI SRA fallback, then kallisto quantification. Progress is recorded in `pipeline_progress.db` (SQLite) under the data root.

The producer does not run `merge → wsfilter → finalize → sanity`. After downloading results, run those locally with `projects/hymenoptera_amalgkit/scripts/run_all_finalization.sh --data-root "$AMALGKIT_DATA_ROOT"`.

Without Docker, `scripts/cloud/vm_setup.sh` installs the tools natively and starts `scripts/rna/run_all_species.py` for every species config. `scripts/cloud/prep_genomes.py` can build all kallisto indexes before a run.

### Container mounts

`cloud_startup.sh` bind-mounts `/opt/MetaInformAnt/projects` at `/app/projects` and `/opt/MetaInformAnt/config` at `/app/config`; the image is built from the repository `Dockerfile` and tagged with `docker_image` (default `metainformant-pipeline`). The Dockerfile also declares `VOLUME /app/output`, which the startup script does not bind-mount.

The pipeline resumes from the SQLite progress database and the quantification outputs under its data root, so keep the data root on a host bind mount. Data written to a Docker-managed volume or to the container's own filesystem is lost when that volume or container is removed.

## Data Flow

1. **Reference**: the species YAML names `genome.ftp_url` and `genome.files.transcriptome_fasta`; when no index exists yet, the orchestrator downloads that transcriptome FASTA and builds `<species>.idx` with `kallisto index`.
2. **Acquisition**: samples are chosen by `metadata` and `select` (and the `max_gb` size cap), then downloaded from ENA, falling back to NCBI SRA.
3. **Quantification**: Amalgkit quantifies each sample with kallisto and writes its abundance table.
4. **Downstream** (after the producer stops): `merge → wsfilter → finalize → sanity`.

The full per-species plan is `metadata → select → getfastq → integrate → quant → merge → wsfilter → finalize → sanity`.

## Result Download

`download_results()` runs `scripts/cloud/download_results.sh`, which copies `/opt/MetaInformAnt/projects` on the VM to `/tmp/quant_export`, drops `*.fastq.gz`, `*.sra` and `index` entries from that copy, transfers it with `gcloud compute scp --recurse` and removes the temporary copy. Everything else under `projects/` is downloaded.

If the script file is missing, the method falls back to `gcloud compute scp --recurse` of these paths under `/opt/MetaInformAnt/projects/hymenoptera_amalgkit/data/`:
- `*/work/quant`: per-sample quantification output
- `*/merged`: merged expression data
- `pipeline_progress.db`: SQLite status database
