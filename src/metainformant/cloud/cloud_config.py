"""Cloud deployment configuration for GCP VMs."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# Sensible defaults for RNA-seq pipeline workloads
_DEFAULT_ZONE = "us-central1-a"
_DEFAULT_MACHINE = "n2-highcpu-96"  # 96 vCPUs, 96 GB RAM — $3.40/hr (spot: ~$0.80/hr)
_DEFAULT_DISK_GB = 500              # SSD for fast I/O
_DEFAULT_IMAGE_FAMILY = "debian-12"
_DEFAULT_IMAGE_PROJECT = "debian-cloud"


@dataclass
class CloudConfig:
    """Configuration for a GCP compute instance running the pipeline.

    Attributes:
        project: GCP project ID.
        zone: GCP zone (e.g. us-central1-a).
        instance_name: VM instance name.
        machine_type: GCP machine type.
        disk_size_gb: Boot disk size in GB.
        spot: Use spot/preemptible pricing (~75% cheaper).
        max_gb: Max sample size in GB for the pipeline.
        workers: Number of parallel download/quant workers.
        threads: Total CPU threads for pipeline.
        gcs_bucket: Optional GCS bucket for result sync.
        repo_url: Git URL to clone on the VM.
        repo_branch: Git branch to checkout.
        docker_image: Name for the Docker image built on-VM.
        service_account_email: Optional SA for the VM.
    """

    project: str = ""
    zone: str = _DEFAULT_ZONE
    instance_name: str = "metainformant-pipeline"
    machine_type: str = _DEFAULT_MACHINE
    disk_size_gb: int = _DEFAULT_DISK_GB
    spot: bool = True
    max_gb: float = 20.0
    workers: int = 80
    threads: int = 96
    gcs_bucket: str = ""
    repo_url: str = "https://github.com/docxology/MetaInformAnt.git"
    repo_branch: str = "main"
    docker_image: str = "metainformant-pipeline"
    service_account_email: str = ""

    # Local paths (relative to project root)
    config_dir: str = "config/amalgkit"
    output_dir: str = "output/amalgkit"

    # Image
    image_family: str = _DEFAULT_IMAGE_FAMILY
    image_project: str = _DEFAULT_IMAGE_PROJECT

    @property
    def startup_script_path(self) -> Path:
        """Path to the cloud startup script."""
        return Path(__file__).resolve().parent.parent.parent.parent / "scripts" / "cloud" / "cloud_startup.sh"

    def validate(self) -> list[str]:
        """Return list of validation errors (empty = valid)."""
        errors = []
        if not self.project:
            errors.append("GCP project ID is required (--project)")
        if self.workers < 1:
            errors.append("Workers must be >= 1")
        if self.threads < 1:
            errors.append("Threads must be >= 1")
        if self.disk_size_gb < 100:
            errors.append("Disk size should be >= 100 GB for pipeline data")
        return errors

    def to_metadata(self) -> dict[str, str]:
        """Convert pipeline params to GCP instance metadata key-value pairs."""
        return {
            "pipeline-max-gb": str(self.max_gb),
            "pipeline-workers": str(self.workers),
            "pipeline-threads": str(self.threads),
            "pipeline-repo-url": self.repo_url,
            "pipeline-repo-branch": self.repo_branch,
            "pipeline-gcs-bucket": self.gcs_bucket,
            "pipeline-docker-image": self.docker_image,
        }
