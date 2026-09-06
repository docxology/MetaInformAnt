# Cloud Module Technical Specification

## Module: `metainformant.cloud`

**Status:** Production-ready (used by the live RNA-seq producer)
**Python:** 3.11+
**Dependencies:** stdlib only (`json`, `logging`, `shutil`, `subprocess`, `time`, `pathlib`); requires the `gcloud` CLI on `PATH` at runtime.

---

## API Reference

### `CloudConfig`

```python
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
    config_dir: str = "config/amalgkit"
    output_dir: str = "output/amalgkit"
    image_family: str = "debian-12"
    image_project: str = "debian-cloud"

    @property
    def startup_script_path(self) -> Path: ...
    def validate(self) -> list[str]: ...
    def to_metadata(self) -> dict[str, str]: ...
```

`validate()` returns a list of human-readable error strings:

- Missing `project` ("GCP project ID is required (--project)")
- `workers < 1`, `threads < 1`
- `disk_size_gb < 100` (warning-level: pipeline data needs space)

`to_metadata()` emits the `pipeline-*` metadata keys consumed by
`scripts/cloud/cloud_startup.sh`: `pipeline-max-gb`, `pipeline-workers`,
`pipeline-threads`, `pipeline-repo-url`, `pipeline-repo-branch`,
`pipeline-gcs-bucket`, `pipeline-docker-image`.

### `GCPDeployer`

```python
class GCPDeployer:
    def __init__(self, config: CloudConfig) -> None: ...

    @staticmethod
    def gcloud_installed() -> bool: ...

    def create_vm(self, dry_run: bool = False) -> dict[str, Any]: ...
    def delete_vm(self) -> bool: ...
    def stop_vm(self) -> bool: ...
    def start_vm(self) -> bool: ...

    def get_vm_status(self) -> dict[str, Any]: ...
    def get_pipeline_status(self) -> str: ...
    def tail_logs(self, lines: int = 50) -> str: ...
    def get_startup_log(self) -> str: ...

    def download_results(self, local_dir: str = "output/amalgkit") -> bool: ...
    def sync_to_gcs(self) -> bool: ...

    def wait_for_ssh(self, max_wait: int = 300) -> bool: ...
    def full_deploy(self) -> dict[str, Any]: ...
```

**Behavioral contract:**

- All `gcloud` invocations go through `_run()` which appends
  `--project <cfg.project> --format json --quiet`.
- `create_vm(dry_run=True)` returns `{"dry_run": True, "command": "<gcloud …>"}`
  without executing. On success it parses the `--format json` instance list and
  returns the first instance dict (or `{"status": "created"}`).
- `create_vm` raises `ValueError` on invalid config and `FileNotFoundError`
  when `startup_script_path` does not exist.
- `delete_vm`/`stop_vm`/`start_vm` return `False` on
  `subprocess.CalledProcessError`; other subprocess failures propagate.
- `get_vm_status()` returns `{"status": "NOT_FOUND"}` when describe exits non-zero.
- SSH-based methods (`get_pipeline_status`, `tail_logs`, `get_startup_log`)
  return the remote stdout, `"No output"`, or an `"SSH failed…"` message —
  they never raise.
- `wait_for_ssh` polls `gcloud compute ssh … echo ok` every 10 s for up to
  `max_wait` seconds.
- `full_deploy` returns `{"vm": <create_vm result>, "ssh_ready": <bool>,
  "status": <describe status or "UNKNOWN">}`; `ssh_ready` reflects
  `wait_for_ssh`.

---

## CLI

`scripts/cloud/deploy_gcp.py` subcommands: `deploy`, `status`, `logs`,
`startup-log`, `download`, `stop`, `start`, `destroy`. `deploy` builds a
`CloudConfig` from `--project/--zone/--machine-type/--disk-gb/--spot/--max-gb/
--workers/--threads/--name/--gcs-bucket/--dry-run`; the remaining commands take
`--project/--zone/--name` and delegate to the matching `GCPDeployer` method.

`scripts/cloud/download_results.sh` accepts `--output DIR`;
`GCPDeployer.download_results` calls it when present and falls back to direct
`gcloud compute scp` of `*/work/quant`, `*/merged`, and `pipeline_progress.db`
from `/opt/MetaInformAnt/projects/hymenoptera_amalgkit/data/` otherwise.

---

## Testing Strategy

### Unit Tests (fast, no GCP)

Real subprocess paths exercised via a stub `gcloud` shell script placed first
on `PATH` (see `tests/cloud/test_gcp_deployer_depth.py`,
`tests/cloud/test_download_results_depth.py`); `CloudConfig` and static
methods tested directly (`tests/cloud/test_cloud.py`). No in-process
test doubles; the stub is a real external `gcloud` process.

---

## Related Modules

- **RNA documentation** (`docs/rna/`) — amalgkit pipeline (primary consumer)
- **Cloud scripts** (`scripts/cloud/`) — CLI, startup script, download/sync helpers
- **[core/io/download.py](../../../src/metainformant/core/io/download.py)** — shared download utilities
