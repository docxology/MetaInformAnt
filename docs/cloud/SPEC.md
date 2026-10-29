# Cloud Module Technical Specification

## Module: `metainformant.cloud`

**Python:** 3.11+
**Dependencies:** standard library only. The module runs the `gcloud` CLI on the local machine; the VM startup script uses `apt-get`, Docker and `git`, plus `gsutil` when a bucket is configured.
**Exports:** `CloudConfig`, `GCPDeployer`

There is no file-based configuration: a `CloudConfig` is built in Python or from the flags of `scripts/cloud/deploy_gcp.py`.

---

## API Reference

### `CloudConfig`

A mutable dataclass (`@dataclass`, not frozen) in `cloud_config.py`:

```python-snippet
@dataclass
class CloudConfig:
    project: str = ""
    zone: str = "us-central1-a"
    instance_name: str = "metainformant-pipeline"
    machine_type: str = "n2-highcpu-96"
    disk_size_gb: int = 500
    local_ssd_count: int = 0
    spot: bool = True
    max_gb: float = 20.0
    workers: int = 80
    threads: int = 96
    gcs_bucket: str = ""
    repo_url: str = "https://github.com/docxology/MetaInformAnt.git"
    repo_branch: str = "main"
    docker_image: str = "metainformant-pipeline"
    service_account_email: str = ""
    config_dir: str = "config/amalgkit"   # declared, not read by GCPDeployer
    output_dir: str = "output/amalgkit"   # declared, not read by GCPDeployer
    image_family: str = "debian-12"
    image_project: str = "debian-cloud"

    @property
    def startup_script_path(self) -> Path: ...        # <repo>/scripts/cloud/cloud_startup.sh
    def validate(self) -> list[str]: ...              # [] when valid
    def to_metadata(self) -> dict[str, str]: ...      # instance metadata for the startup script
```

`validate()` reports these errors:

| Condition | Message |
|-----------|---------|
| `project` is empty | `GCP project ID is required (--project)` |
| `workers < 1` | `Workers must be >= 1` |
| `threads < 1` | `Threads must be >= 1` |
| `disk_size_gb < 100` | `Disk size should be >= 100 GB for pipeline data` |

`to_metadata()` returns these keys, which the startup script reads. `create_vm` passes only the non-empty values (`--metadata key=value,...`), and the script falls back to its own default for a missing key.

| Metadata key | Source field | Startup script default |
|--------------|--------------|------------------------|
| `pipeline-max-gb` | `max_gb` | `20.0` |
| `pipeline-workers` | `workers` | `80` |
| `pipeline-threads` | `threads` | `96` |
| `pipeline-repo-url` | `repo_url` | the MetaInformAnt repository URL |
| `pipeline-repo-branch` | `repo_branch` | `main` |
| `pipeline-gcs-bucket` | `gcs_bucket` | empty (no sync) |
| `pipeline-docker-image` | `docker_image` | `metainformant-pipeline` |

### `GCPDeployer`

```python-snippet
class GCPDeployer:
    def __init__(self, config: CloudConfig) -> None: ...   # the config is kept as self.cfg

    @staticmethod
    def gcloud_installed() -> bool: ...

    # VM lifecycle
    def create_vm(self, dry_run: bool = False) -> dict[str, Any]: ...
    def delete_vm(self) -> bool: ...
    def stop_vm(self) -> bool: ...
    def start_vm(self) -> bool: ...

    # Status and monitoring
    def get_vm_status(self) -> dict[str, Any]: ...
    def get_pipeline_status(self) -> str: ...
    def tail_logs(self, lines: int = 50) -> str: ...
    def get_startup_log(self) -> str: ...

    # Data transfer
    def download_results(self, local_dir: str = "output/amalgkit") -> bool: ...
    def sync_to_gcs(self) -> bool: ...

    # Convenience
    def wait_for_ssh(self, max_wait: int = 300) -> bool: ...
    def full_deploy(self) -> dict[str, Any]: ...
```

Every `gcloud` call is `gcloud <args> --project <project> --format json --quiet`. SSH commands are `gcloud compute ssh <instance_name> --zone <zone> --command <command>`.

### Return values and failures

| Method | Runs | Returns | Failure |
|--------|------|---------|---------|
| `create_vm` | validates the config, then `gcloud compute instances create` (300 s timeout) | the first instance object of gcloud's JSON output, or `{"status": "created"}` if gcloud printed nothing; `{"dry_run": True, "command": "<gcloud command>"}` with `dry_run=True`, and nothing is run | `ValueError("Invalid config: ...")`; `FileNotFoundError("Startup script not found: ...")`; `subprocess.CalledProcessError` or `subprocess.TimeoutExpired` from gcloud |
| `delete_vm` | `gcloud compute instances delete <name> --delete-disks all` (120 s) | `True` | `False` after logging the error when gcloud exits non-zero; `TimeoutExpired` propagates |
| `stop_vm`, `start_vm` | `gcloud compute instances stop` or `start` (120 s) | `True` | `False` when gcloud exits non-zero; `TimeoutExpired` propagates |
| `get_vm_status` | `gcloud compute instances describe` | the parsed JSON dict (`{}` if gcloud printed nothing) | `{"status": "NOT_FOUND"}` whenever gcloud exits non-zero; `TimeoutExpired` propagates |
| `get_pipeline_status` | over SSH: `cd /opt/MetaInformAnt && python3 scripts/rna/check_pipeline_status.py` (30 s), printing `Pipeline not ready` if it fails | stdout text, or `No output` | `SSH failed — VM may still be starting up` |
| `tail_logs` | over SSH: `sudo docker logs --tail <lines> metainformant-pipeline` (30 s) | log text; `No logs found yet` when docker fails; `No output` when nothing is printed | `SSH failed` |
| `get_startup_log` | over SSH: `journalctl -u google-startup-scripts.service -n 80`, falling back to `tail -80 /var/log/syslog` (30 s) | log text, or `No output` | `SSH failed` |
| `download_results` | creates `local_dir`, then `bash scripts/cloud/download_results.sh --output <local_dir>` (7200 s); direct `gcloud compute scp` calls only if the script file is missing | `True` | `False` when a timeout expires; `subprocess.CalledProcessError` when the script exits non-zero |
| `sync_to_gcs` | over SSH: `gsutil -m rsync -r /opt/MetaInformAnt/projects/ gs://<bucket>/amalgkit/` (600 s) | `True` | `False` without `gcs_bucket`, on a non-zero exit or on a timeout |
| `wait_for_ssh` | `echo ok` over SSH every 10 s | `True` once it succeeds | `False` after `max_wait` seconds |
| `full_deploy` | `create_vm()`, `wait_for_ssh(300)`, `get_vm_status()`; prints a summary | `{"vm": ..., "ssh_ready": bool, "status": str}` | the exceptions of `create_vm` propagate |

### The instance `create_vm` builds

```text
gcloud compute instances create <instance_name>
  --zone <zone> --machine-type <machine_type>
  --boot-disk-size <disk_size_gb>GB --boot-disk-type pd-standard
  --image-family <image_family> --image-project <image_project>
  --metadata-from-file startup-script=<startup_script_path>
  --metadata <key=value,...>
  --scopes cloud-platform --tags metainformant-pipeline
  [--local-ssd interface=NVME]            repeated local_ssd_count times
  [--provisioning-model SPOT --instance-termination-action STOP]   when spot
  [--service-account <service_account_email>]                      when set
  --project <project> --format json --quiet
```

The boot disk type is always `pd-standard`. No network, firewall or external-address option is set, so gcloud defaults apply.

---

## Startup script (`scripts/cloud/cloud_startup.sh`)

Runs under `set -euo pipefail` and logs to `/var/log/metainformant-startup.log`. Steps:

1. Read the `pipeline-*` metadata keys above from the instance metadata server.
2. Install `docker.io` with `apt-get` if `docker` is missing.
3. Install `mdadm`; if local NVMe SSDs are present, assemble them as RAID 0 (`/dev/md0`), format ext4 and mount them on `/opt/MetaInformAnt`; otherwise create that directory on the boot disk.
4. Clone `pipeline-repo-url` at `pipeline-repo-branch` with `--depth 1` into `/opt/MetaInformAnt` (or fetch, check out and pull when a clone exists), then `docker build -t <image> .`.
5. Start the container `metainformant-pipeline` (`docker run -d --restart unless-stopped --cpus=<threads> --memory=90g`, mounts `/opt/MetaInformAnt/projects` at `/app/projects` and `/opt/MetaInformAnt/config` at `/app/config`, environment `PIPELINE_MAX_GB`, `PIPELINE_WORKERS`, `PIPELINE_THREADS`) running `python3 scripts/rna/process_species.py --species apis_mellifera --config-dir projects/hymenoptera_amalgkit/config/amalgkit --data-root "$AMALGKIT_DATA_ROOT"`.
6. With a bucket configured, write `/etc/cron.d/metainformant-sync`, which runs `gsutil -m rsync -r /opt/MetaInformAnt/projects/ gs://<bucket>/projects/` every 30 minutes and logs to `/var/log/metainformant-sync.log`.

The script does not define `AMALGKIT_DATA_ROOT`, and `set -u` makes the expansion in step 5 fatal when the variable is unset. The shallow clone in step 4 leaves the `projects/hymenoptera_amalgkit` submodule empty.

## Command-line interface (`scripts/cloud/deploy_gcp.py`)

`gcloud` must be on `PATH` for every subcommand, `deploy --dry-run` included; otherwise the script prints an install hint and exits with status 1.

| Subcommand | Flags (default) |
|------------|-----------------|
| `deploy` | `--project` (required), `--zone` (`us-central1-a`), `--name` (`metainformant-pipeline`), `--machine-type` (`n2-standard-32`), `--disk-gb` (1000), `--local-ssd-count` (0), `--no-spot` (Spot is the default), `--max-gb` (20.0), `--workers` (28), `--threads` (32), `--gcs-bucket` (none), `--dry-run` |
| `status` | `--project`, `--zone`, `--name` |
| `logs` | `--project`, `--zone`, `--name`, `--lines` (50) |
| `startup-log` | `--project`, `--zone`, `--name` |
| `download` | `--project`, `--zone`, `--name`, `--output` (`output/amalgkit`) |
| `stop`, `start`, `destroy` | `--project`, `--zone`, `--name` |

`deploy` calls `full_deploy()`, or prints the `gcloud` command with `--dry-run`. `status` prints `VM Status: <status>` and, for a `RUNNING` VM, the pipeline status text. `destroy` asks `[y/N]` before calling `delete_vm()`. For commands other than `deploy` the flags `--project`, `--zone` and `--name` default to an empty project, `us-central1-a` and `metainformant-pipeline`.

## Download script (`scripts/cloud/download_results.sh`)

```text
bash scripts/cloud/download_results.sh [--output DIR] [--species SPECIES]
```

`--output` defaults to `output/amalgkit`; `--species` is parsed but not used. The project, zone and VM come from the environment: `GCP_PROJECT` (default `cryptoptera`), `GCP_ZONE` (default `us-central1-a`) and `GCP_VM` (default `metainformant-pipeline`). The script:

1. over SSH, copies `/opt/MetaInformAnt/projects` into `/tmp/quant_export` on the VM (copy errors are ignored), deletes `*.fastq.gz`, `*.sra` and `index` entries from that copy, and fixes permissions;
2. downloads the copy with `gcloud compute scp --recurse` into `--output`;
3. removes `/tmp/quant_export` on the VM.

## Companion scripts (`scripts/cloud/`)

| Script | Purpose |
|--------|---------|
| `install_gcloud.sh` | Install the Google Cloud CLI on Linux: from the apt repository when `apt-get` exists, otherwise from the standalone tarball. Exits early if `gcloud` is already installed. |
| `prep_genomes.py` | For each `amalgkit_*.yaml` in `--config-dir` (default `config/amalgkit`; names containing `template` or `test` are skipped), download `genome.files.transcriptome_fasta` from `genome.ftp_url` with `wget` and build `<species>.idx` with `kallisto index -i`. `--threads` is accepted but not used. |
| `vm_setup.sh` | Native setup without Docker: installs fastp 0.24.0 and the SRA Toolkit 3.2.1 when missing, runs `pip3 install --break-system-packages -e .` plus an unpinned amalgkit from `git+https://github.com/kfuku52/amalgkit` (the project pins `v0.16.60` in the `rna` extra), prints tool versions, and starts `scripts/rna/run_all_species.py --max-gb 20.0 --workers 24 --threads 32` in the background. kallisto must already be installed. |
| `sync_cloud_down.sh` | Loop that runs `rsync` every 60 s from the SSH host alias `metainformant-pipeline.us-central1-a.cryptoptera` (created by `gcloud compute config-ssh`), skipping FASTQ and SRA files. |
| `categorize_errors.py`, `print_table.py` | Query `/app/output/amalgkit/pipeline_progress.db` inside the `metainformant-pipeline` container over `gcloud compute ssh` (VM `metainformant-pipeline`, zone `us-central1-a`, project `cryptoptera`, all hard-coded) and print grouped failure messages or a per-species state table. |

---

## Security Model

1. **Authentication:** the account `gcloud` is logged in with (`gcloud auth login`); the module stores no credentials.
2. **VM identity:** `--scopes cloud-platform` lets the VM's service account (the default one unless `service_account_email` is set) call any Google API its IAM roles allow.
3. **Network:** no firewall rule is created; SSH runs through `gcloud compute ssh` and `scp`.
4. **Data in transit:** SSH for command execution and file transfer; HTTPS for the repository clone and for `gsutil`.

---

## Testing

The tests in `tests/cloud/` use real files and a stub `gcloud` executable placed first on `PATH`, so they need no GCP credentials or network access. `pyproject.toml` registers no marker for runs against a live GCP project.

---

## Related Modules

- **[rna/](../rna/)** — the amalgkit RNA-seq pipeline the VM runs
- **[gwas/](../gwas/)** — no cloud runner; its reference-genome helpers in `metainformant.gwas.data.download` can run on a VM
