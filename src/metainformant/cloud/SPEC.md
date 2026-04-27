# Cloud Module Technical Specification

## Module: `metainformant.cloud`

**Status:** Production-ready
**Version:** 0.3.0
**Python:** 3.11+
**Dependencies:** google-cloud-compute, google-cloud-storage (optional; falls back to `gcloud` CLI)

---

## API Reference

### `CloudConfig`

```python
@dataclass
class CloudConfig:
    project_id: str
    region: str = "us-central1"
    zone: str = "us-central1-a"
    machine_type: str = "n1-standard-8"
    disk_size_gb: int = 100
    boot_disk_type: str = "pd-ssd"
    docker_image: str | None = None
    preemptible: bool = False
    network: str = "default"
    service_account: str | None = None
    scopes: list[str] = field(default_factory=lambda: ["https://www.googleapis.com/auth/cloud-platform"])
    startup_script: str | Path = "scripts/cloud/cloud_startup.sh"
    shutdown_script: str | None = None
    metadata: dict[str, str] | None = None
    labels: dict[str, str] | None = None
    ssh_key_path: str | None = None
    max_retries: int = 3
    retry_delay_seconds: int = 30
```

### `GCPDeployer`

```python
class GCPDeployer:
    def __init__(self, config: CloudConfig, retry_config: RetryConfig | None = None,
                 enable_cloud_logging: bool = True) -> None: ...

    # Instance lifecycle
    def create_instance(self, name: str, labels: dict | None = None) -> ComputeInstance: ...
    def delete_instance(self, instance_name: str) -> None: ...
    def get_instance(self, instance_name: str) -> ComputeInstance | None: ...
    def list_instances(self) -> list[ComputeInstance]: ...
    def wait_for_ready(self, instance_name: str, timeout_seconds: int = 600) -> bool: ...

    # File transfer
    def upload_file(self, local_path: str | Path, remote_path: str, instance_name: str) -> None: ...
    def upload_directory(self, local_path: str | Path, remote_path: str, instance_name: str) -> None: ...
    def download_file(self, remote_path: str, local_path: str | Path, instance_name: str) -> None: ...
    def download_directory(self, remote_path: str, local_path: str | Path, instance_name: str) -> None: ...

    # Command execution
    def execute_command(self, command: str, instance_name: str, capture_output: bool = True,
                        timeout: int | None = None) -> subprocess.CompletedProcess: ...

    # Docker operations
    def build_docker_image(self, local_context: str | Path, tag: str, instance_name: str) -> str: ...
    def push_docker_image(self, image: str, registry: str, instance_name: str) -> str: ...
    def run_docker_container(self, image: str, command: list[str] | None = None,
                              volumes: dict | None = None, instance_name: str) -> Container: ...

    # Static utilities
    @staticmethod
    def list_available_machine_types(zone: str) -> list[str]: ...
    @staticmethod
    def get_current_cost(project_id: str) -> float: ...
    @staticmethod
    def create_budget_alert(project_id: str, amount_usd: float, email: str) -> None: ...
```

### `genome_prep.py`

```python
def prepare_genome_index(
    accession: str,
    output_dir: str | Path,
    tools: list[str] | None = None,
    overwrite: bool = False,
    bucket_name: str | None = None,
    instance_name: str | None = None
) -> dict[str, Path]: ...

def download_genome_from_ncbi(
    accession: str,
    output_dir: str | Path,
    include: list[str] | None = None,
    overwrite: bool = False
) -> Path: ...

def upload_to_cache(genome_dir: str | Path, bucket_name: str, accession: str) -> list[str]: ...
def download_from_cache(accession: str, bucket_name: str, output_dir: str | Path) -> Path | None: ...
```

---

## Type Signatures

### Data Structures

```python
from dataclasses import dataclass
from typing import TypedDict

@dataclass
class ComputeInstance:
    name: str
    zone: str
    status: str  # "RUNNING", "STOPPED", "TERMINATED"
    internal_ip: str
    external_ip: str | None
    creation_timestamp: datetime
    machine_type: str

class CloudLoggingEntry(TypedDict):
    timestamp: str
    severity: str  # "INFO", "WARNING", "ERROR"
    message: str
    labels: dict[str, str]

class InstanceStatus(TypedDict):
    name: str
    status: str
    cpu_utilization: float
    memory_usage_gb: float
    disk_io_read: int
    disk_io_write: int
```

---

## Configuration Schema (YAML)

```yaml
# Required
project_id: str
# Optional with defaults
region: str = "us-central1"
zone: str = "us-central1-a"
machine_type: str = "n1-standard-8"
disk_size_gb: int = 100
# Flags
preemptible: bool = False
# Paths
startup_script: str = "scripts/cloud/cloud_startup.sh"
# Containerization
docker_image: str | None = None
docker_registry: str | None = None
# IAM
service_account: str | None = None
scopes: list[str] = ["https://www.googleapis.com/auth/cloud-platform"]
# Metadata
metadata: dict[str, str] | None = None
labels: dict[str, str] | None = None
# Retry policy
max_retries: int = 3
retry_delay_seconds: int = 30
retryable_errors: list[str] = ["RESOURCE_EXHAUSTED", "UNAVAILABLE"]
# Monitoring
enable_stackdriver: bool = True
log_level: str = "INFO"
```

---

## Error Codes & Messages

| Error | Condition | Recovery |
|-------|-----------|----------|
| `CLOUD-AUTH-001` | GCP credentials not found | Run `gcloud auth application-default login` |
| `CLOUD-QUOTA-001` | CPU quota exceeded in zone | Request quota increase or use different zone |
| `CLOUD-INSTANCE-001` | VM creation failed | Check startup script syntax at `/var/log/startupscript.log` |
| `CLOUD-COMMAND-001` | Remote command non-zero exit | Inspect `stderr` from `execute_command()` |
| `CLOUD-TRANSFER-001` | File transfer timeout | Increase `ssh_timeout_seconds` or check network |
| `CLOUD-COST-001` | Estimated cost exceeds budget (optional guard) | Reduce machine type or runtime |

---

## Performance Characteristics

| Operation | Latency | Throughput |
|-----------|---------|------------|
| VM creation | 60–120s | N/A |
| Docker build & push | 2–10 min | N/A |
| Genome download (100MB from NCBI) | 30–120s | ~1 MB/s |
| Genome from GCS cache | 5–15s | ~50 MB/s |
| FASTQ upload (10 GB) | 5–15 min | ~10-20 MB/s (home upload) |
| RNA-seq pipeline | 2–8 hours | N/A (CPU-bound, ~16-32 cores) |
| Results download (1 GB) | 1–3 min | ~5–10 MB/s |

**Parallelism:** Multiple VMs can be launched simultaneously (quota-limited).

---

## Security Model

1. **Authentication**: Application Default Credentials or explicit service account key
2. **Authorization**: IAM roles: `compute.instanceAdmin.v1`, `storage.objectAdmin`
3. **Network**: VPC firewall rules (SSH port 22, internal IP only by default)
4. **Data at rest**: Encrypted with Google-managed keys (or CMEK if configured)
5. **Data in transit**: SSH for file transfer, HTTPS for GCS/Container Registry

---

## Testing Strategy

### Unit Tests (Fast, No GCP)

Mocked `gcloud` CLI calls via `subprocess.run` patch (allowed because mocking external CLI is not the same as mocking internal functions).

### Integration Tests (Real GCP)

Marked `@pytest.mark.cloud_integration` — requires:
- `GCP_PROJECT_ID` environment variable
- Credentials configured (`gcloud auth application-default login`)
- Billing enabled on project

Cost: ~$0.50–2.00 per full integration test run.

---

## Related Modules

- **[rna/](rna/)** — Amalgkit RNA-seq pipeline (primary cloud consumer)
- **[gwas/](gwas/)** — GWAS analysis (can be deployed to cloud)
- **[core/io/download.py](core/io/download.py)** — Download utilities (used by genome_prep)

---

## Changelog

**v0.3.0** (2025-04-15)
- Initial production release
- GCP Compute Engine integration
- Docker container orchestration
- Genome caching system
- Preemptible VM support

**v0.2.5** (2025-01-10)
- Added RetryConfig with exponential backoff
- Cloud Logging integration
- Cost tracker with budget alerts

**v0.2.0** (2024-11-01)
- First public preview (alpha)
- Basic VM creation/teardown
- File transfer via SCP
