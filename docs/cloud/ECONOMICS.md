# Cloud Economics and Cost Management

Running the amalgkit pipeline on Google Cloud Platform is billed per VM-hour and per disk. This document covers what the METAINFORMANT cloud module sets (provisioning model, disk type and size), what it does not provide, and the pipeline settings that bound run time.

## Compute Economics: Standard vs. Spot Provisioning

`CloudConfig.spot` defaults to `True`, and `deploy_gcp.py deploy` uses Spot unless `--no-spot` is passed. The two provisioning models differ in price and in what happens when Google needs the capacity back.

### Spot Provisioning

- **Cost**: about 75% cheaper than on demand according to the `CloudConfig` docstring. The comment on the default `n2-highcpu-96` machine gives about $3.40/hr on demand and about $0.80/hr as a Spot VM.
- **Behavior**: `create_vm` adds `--provisioning-model SPOT --instance-termination-action STOP`. Google can stop the VM when it needs the capacity; a stopped VM keeps its disk and `start_vm()` starts it again. The startup script runs the container with `--restart unless-stopped`.
- **Use Case**: short, fault-tolerant runs.
- **Experience**: on the large (8,000+ sample) RNA-seq runs, Spot VMs were repeatedly stopped during quantification and the pipeline stalled until the VM was restarted.

### Standard Provisioning (On-Demand)

- **Cost**: the on-demand list price (about $3.40/hr for the default `n2-highcpu-96`, per the comment in `cloud_config.py`).
- **Behavior**: runs until it is stopped or deleted.
- **Use Case**: long quantification runs that cannot tolerate interruption; select it with `spot=False` or `--no-spot`.
- **Recommendation**: use standard provisioning for production runs of the full cohort. The code does not enforce this; Spot stays the default.

## Storage Economics: Disks

- **Boot disk type**: `create_vm` always passes `--boot-disk-type pd-standard`, and `CloudConfig` has no field to change it. The comment on `_DEFAULT_DISK_GB` in `cloud_config.py` says SSD, but the command line does not request one.
- **Boot disk size**: `disk_size_gb` is 500 GB in `CloudConfig` and 1000 GB in `deploy_gcp.py deploy`; `validate()` rejects anything below 100 GB.
- **Local SSDs**: `local_ssd_count` attaches that many NVMe local SSDs (`--local-ssd interface=NVME`, 375 GB each per the CLI help). The startup script combines them as RAID 0 and mounts them on `/opt/MetaInformAnt`.
- **FASTQ cleanup**: by default (`AMALGKIT_RECLAIM_RAW_AFTER_QUANT=yes`) the producer removes a sample's raw FASTQ/SRA inputs once its quantification has verified provenance.
- **Deletion**: `delete_vm()` and `deploy_gcp.py destroy` pass `--delete-disks all`, so the boot disk is deleted with the VM; download the results first. `stop_vm()` keeps the disk, and a kept disk is still billed. The module creates no separate data disk and never sets `--auto-delete`.

## Pipeline Settings That Bound Cost

Run time, and therefore cost, depends on the work the pipeline accepts:

1. **Sample size cap**: `CloudConfig.max_gb` (CLI `--max-gb`, passed to the container as `PIPELINE_MAX_GB`) defaults to 20.0. Samples whose `total_bases` exceed `max_gb` x 10^9 are skipped. Outside the cloud module, `run_all_species.py` and `process_species.py` default to 50.0 when `PIPELINE_MAX_GB` is unset.
2. **Empty downloads**: a downloaded FASTQ file of 0 bytes marks that sample as failed before it reaches quantification.
3. **Concurrency**: `workers` and `threads` (CLI `--workers` and `--threads`) set the number of parallel samples and the total thread budget; the startup script also uses `threads` as the container's `--cpus` limit.

## Summary

- Spot provisioning is the default; use `--no-spot` for long production runs.
- Download results before `destroy`: the boot disk is deleted with the VM.
- The module has no cost estimator, budget alert or automatic shutdown; delete the VM when the run is done.
- VMs are created with the `cloud-platform` access scope, so processes on the VM can call any Google API its service account may use, and that usage is billed to the project.
