"""GCP VM lifecycle management for MetaInformAnt pipelines.

Manages creating, monitoring, and tearing down GCP Compute Engine VMs
to run the amalgkit RNA-seq pipeline at scale.  Shells out to `gcloud`
CLI so no Python SDK is needed.
"""
from __future__ import annotations

import json
import logging
import shutil
import subprocess
import sys
import textwrap
import time
from pathlib import Path
from typing import Any, Optional

from metainformant.cloud.cloud_config import CloudConfig

logger = logging.getLogger(__name__)


class GCPDeployer:
    """Create and manage a GCP VM running the pipeline."""

    def __init__(self, config: CloudConfig) -> None:
        self.cfg = config

    # ── helpers ───────────────────────────────────────────────────────────

    @staticmethod
    def gcloud_installed() -> bool:
        """Check if gcloud CLI is available."""
        return shutil.which("gcloud") is not None

    def _run(self, args: list[str], *, check: bool = True,
             capture: bool = True, timeout: int = 120) -> subprocess.CompletedProcess:
        """Run a gcloud command."""
        cmd = ["gcloud"] + args + [
            "--project", self.cfg.project,
            "--format", "json",
            "--quiet",
        ]
        logger.debug("Running: %s", " ".join(cmd))
        return subprocess.run(
            cmd,
            capture_output=capture,
            text=True,
            check=check,
            timeout=timeout,
        )

    def _ssh_run(self, command: str, *, timeout: int = 60) -> subprocess.CompletedProcess:
        """Run a command on the VM via SSH."""
        return self._run([
            "compute", "ssh", self.cfg.instance_name,
            "--zone", self.cfg.zone,
            "--command", command,
        ], timeout=timeout)

    # ── VM lifecycle ─────────────────────────────────────────────────────

    def create_vm(self, dry_run: bool = False) -> dict[str, Any]:
        """Create a GCP VM with the pipeline startup script.

        Returns:
            Dict with instance details on success.
        """
        errors = self.cfg.validate()
        if errors:
            raise ValueError(f"Invalid config: {'; '.join(errors)}")

        startup_script = self.cfg.startup_script_path
        if not startup_script.exists():
            raise FileNotFoundError(f"Startup script not found: {startup_script}")

        # Build metadata string
        metadata_items = self.cfg.to_metadata()
        metadata_str = ",".join(f"{k}={v}" for k, v in metadata_items.items() if v)

        cmd = [
            "compute", "instances", "create", self.cfg.instance_name,
            "--zone", self.cfg.zone,
            "--machine-type", self.cfg.machine_type,
            "--boot-disk-size", f"{self.cfg.disk_size_gb}GB",
            "--boot-disk-type", "pd-standard",
            "--image-family", self.cfg.image_family,
            "--image-project", self.cfg.image_project,
            "--metadata-from-file", f"startup-script={startup_script}",
            "--metadata", metadata_str,
            "--scopes", "cloud-platform",
            "--tags", "metainformant-pipeline",
        ]

        if self.cfg.spot:
            cmd.extend(["--provisioning-model", "SPOT",
                         "--instance-termination-action", "STOP"])

        if self.cfg.service_account_email:
            cmd.extend(["--service-account", self.cfg.service_account_email])

        if dry_run:
            full_cmd = ["gcloud"] + cmd + ["--project", self.cfg.project]
            return {"dry_run": True, "command": " ".join(full_cmd)}

        result = self._run(cmd, timeout=300)
        instances = json.loads(result.stdout) if result.stdout else []
        return instances[0] if instances else {"status": "created"}

    def delete_vm(self) -> bool:
        """Delete the VM and its disk."""
        try:
            self._run([
                "compute", "instances", "delete", self.cfg.instance_name,
                "--zone", self.cfg.zone,
                "--delete-disks", "all",
            ], timeout=120)
            return True
        except subprocess.CalledProcessError as e:
            logger.error("Failed to delete VM: %s", e.stderr)
            return False

    def stop_vm(self) -> bool:
        """Stop (but keep) the VM."""
        try:
            self._run([
                "compute", "instances", "stop", self.cfg.instance_name,
                "--zone", self.cfg.zone,
            ], timeout=120)
            return True
        except subprocess.CalledProcessError:
            return False

    def start_vm(self) -> bool:
        """Start a stopped VM."""
        try:
            self._run([
                "compute", "instances", "start", self.cfg.instance_name,
                "--zone", self.cfg.zone,
            ], timeout=120)
            return True
        except subprocess.CalledProcessError:
            return False

    # ── status & monitoring ──────────────────────────────────────────────

    def get_vm_status(self) -> dict[str, Any]:
        """Get VM instance status."""
        try:
            result = self._run([
                "compute", "instances", "describe", self.cfg.instance_name,
                "--zone", self.cfg.zone,
            ])
            return json.loads(result.stdout) if result.stdout else {}
        except subprocess.CalledProcessError:
            return {"status": "NOT_FOUND"}

    def get_pipeline_status(self) -> str:
        """Check remote pipeline progress via SSH."""
        try:
            result = self._ssh_run(
                "cd /opt/MetaInformAnt && python3 scripts/rna/check_pipeline_status.py 2>/dev/null || echo 'Pipeline not ready'",
                timeout=30,
            )
            return result.stdout.strip() if result.stdout else "No output"
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            return "SSH failed — VM may still be starting up"

    def tail_logs(self, lines: int = 50) -> str:
        """Tail the pipeline log on the remote VM."""
        try:
            result = self._ssh_run(
                f"tail -{lines} /opt/MetaInformAnt/output/amalgkit/pipeline_restart.log 2>/dev/null "
                f"|| tail -{lines} /opt/MetaInformAnt/pipeline.log 2>/dev/null "
                f"|| echo 'No logs found yet'",
                timeout=30,
            )
            return result.stdout.strip() if result.stdout else "No output"
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            return "SSH failed"

    def get_startup_log(self) -> str:
        """Check the VM startup script log."""
        try:
            result = self._ssh_run(
                "sudo journalctl -u google-startup-scripts.service --no-pager -n 80 2>/dev/null "
                "|| tail -80 /var/log/syslog 2>/dev/null",
                timeout=30,
            )
            return result.stdout.strip() if result.stdout else "No output"
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            return "SSH failed"

    # ── data transfer ────────────────────────────────────────────────────

    def download_results(self, local_dir: str = "output/amalgkit") -> bool:
        """Download pipeline results from the VM.

        Uses gcloud compute scp to recursively copy the output directory.
        """
        local_path = Path(local_dir)
        local_path.mkdir(parents=True, exist_ok=True)

        try:
            # Download quant results, merged data, and the progress DB
            for subdir in ["*/work/quant", "*/merged", "pipeline_progress.db"]:
                subprocess.run([
                    "gcloud", "compute", "scp", "--recurse",
                    "--zone", self.cfg.zone,
                    "--project", self.cfg.project,
                    f"{self.cfg.instance_name}:/opt/MetaInformAnt/output/amalgkit/{subdir}",
                    str(local_path),
                ], check=False, timeout=3600)
            return True
        except subprocess.TimeoutExpired:
            logger.error("Download timed out after 1 hour")
            return False

    def sync_to_gcs(self) -> bool:
        """Trigger a GCS sync on the remote VM (if bucket configured)."""
        if not self.cfg.gcs_bucket:
            logger.warning("No GCS bucket configured")
            return False
        try:
            self._ssh_run(
                f"gsutil -m rsync -r /opt/MetaInformAnt/output/amalgkit/ "
                f"gs://{self.cfg.gcs_bucket}/amalgkit/",
                timeout=600,
            )
            return True
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            return False

    # ── convenience ──────────────────────────────────────────────────────

    def wait_for_ssh(self, max_wait: int = 300) -> bool:
        """Wait until SSH is available on the VM."""
        start = time.time()
        while time.time() - start < max_wait:
            try:
                self._ssh_run("echo ok", timeout=10)
                return True
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
                time.sleep(10)
        return False

    def full_deploy(self) -> dict[str, Any]:
        """End-to-end: create VM, wait for SSH, return status.

        Returns:
            Dict with deployment result and initial status.
        """
        print("🚀 Creating GCP VM...")
        result = self.create_vm()
        print(f"   ✓ VM created: {self.cfg.instance_name} ({self.cfg.machine_type})")
        print(f"   💰 Pricing: {'SPOT' if self.cfg.spot else 'ON-DEMAND'}")
        print(f"   🗺️  Zone: {self.cfg.zone}")

        print("\n⏳ Waiting for SSH access...")
        if self.wait_for_ssh(max_wait=300):
            print("   ✓ SSH ready")
        else:
            print("   ⚠ SSH not ready after 5 min — VM may still be booting")

        print("\n📋 VM Status:")
        status = self.get_vm_status()
        vm_status = status.get("status", "UNKNOWN")
        print(f"   Status: {vm_status}")

        print(f"\n🔧 Pipeline will start automatically via startup script.")
        print(f"   Workers: {self.cfg.workers} | Threads: {self.cfg.threads} | Max GB: {self.cfg.max_gb}")
        print(f"\n   Monitor with: python scripts/cloud/deploy_gcp.py status")
        print(f"   View logs:    python scripts/cloud/deploy_gcp.py logs")

        return {"vm": result, "ssh_ready": True, "status": vm_status}
