"""Depth tests for GCPDeployer subprocess paths using a stub gcloud on PATH.

Real-implementation tests: every test exercises the actual subprocess code
path via a real executable shim written to a temporary bin directory placed
first on PATH (same pattern as tests/gwas/test_gwas_sra_environment.py).
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from metainformant.cloud.cloud_config import CloudConfig
from metainformant.cloud.gcp_deployer import GCPDeployer

STUB = """#!/bin/sh
# Record argv and emit a canned JSON response.
printf '%s\n' "$*" >> "$GCLOUD_CALL_LOG"
if [ -n "$GCLOUD_STDOUT_FILE" ]; then cat "$GCLOUD_STDOUT_FILE"; fi
exit 0
"""


@pytest.fixture()
def stubbed_gcloud(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    """Put a real gcloud stub first on PATH and capture its argv log."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    tool = bin_dir / "gcloud"
    tool.write_text(STUB, encoding="utf-8")
    tool.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    call_log = tmp_path / "calls.log"
    call_log.touch()
    stdout_file = tmp_path / "stdout.json"
    monkeypatch.setenv("GCLOUD_CALL_LOG", str(call_log))
    monkeypatch.setenv("GCLOUD_STDOUT_FILE", str(stdout_file))
    return {"log": call_log, "stdout": stdout_file}


def _deployer(tmp_path: Path, **overrides) -> GCPDeployer:
    cfg = CloudConfig(project="proj", **overrides)
    return GCPDeployer(cfg)


def test_create_vm_dry_run_builds_full_command(tmp_path: Path) -> None:
    deployer = _deployer(tmp_path, local_ssd_count=2, spot=True)
    result = deployer.create_vm(dry_run=True)
    assert result["dry_run"] is True
    cmd = result["command"]
    assert "--provisioning-model SPOT" in cmd
    assert "--instance-termination-action STOP" in cmd
    assert cmd.count("--local-ssd interface=NVME") == 2
    assert "--boot-disk-size" in cmd and "--machine-type" in cmd
    assert "cloud_startup.sh" in cmd


def test_create_vm_without_spot_or_ssd(tmp_path: Path) -> None:
    deployer = _deployer(tmp_path, local_ssd_count=0, spot=False)
    result = deployer.create_vm(dry_run=True)
    cmd = result["command"]
    assert "--provisioning-model" not in cmd
    assert "--local-ssd" not in cmd


def test_create_vm_rejects_invalid_config(tmp_path: Path) -> None:
    deployer = GCPDeployer(CloudConfig(project=""))
    with pytest.raises(ValueError, match="Invalid config"):
        deployer.create_vm()


def test_create_vm_rejects_missing_startup_script(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    deployer = _deployer(tmp_path)
    missing = tmp_path / "nope.sh"
    monkeypatch.setattr(type(deployer.cfg), "startup_script_path", property(lambda self: missing))
    with pytest.raises(FileNotFoundError):
        deployer.create_vm()


def test_create_vm_parses_json_response(tmp_path: Path, stubbed_gcloud: dict[str, Path]) -> None:
    stubbed_gcloud["stdout"].write_text(json.dumps([{"name": "metainformant-pipeline", "status": "RUNNING"}]))
    deployer = _deployer(tmp_path)
    result = deployer.create_vm()
    assert result == {"name": "metainformant-pipeline", "status": "RUNNING"}
    calls = stubbed_gcloud["log"].read_text()
    assert "compute instances create metainformant-pipeline" in calls
    assert "--project proj" in calls


def test_vm_lifecycle_success_paths(tmp_path: Path, stubbed_gcloud: dict[str, Path]) -> None:
    deployer = _deployer(tmp_path)
    assert deployer.delete_vm() is True
    assert deployer.stop_vm() is True
    assert deployer.start_vm() is True
    calls = stubbed_gcloud["log"].read_text()
    assert "compute instances delete metainformant-pipeline --zone us-central1-a --delete-disks all" in calls
    assert "compute instances stop metainformant-pipeline" in calls
    assert "compute instances start metainformant-pipeline" in calls


def test_get_vm_status_returns_describe_json(tmp_path: Path, stubbed_gcloud: dict[str, Path]) -> None:
    stubbed_gcloud["stdout"].write_text(json.dumps({"status": "TERMINATED"}))
    deployer = _deployer(tmp_path)
    assert deployer.get_vm_status() == {"status": "TERMINATED"}


def test_ssh_paths_emit_compute_ssh_calls(tmp_path: Path, stubbed_gcloud: dict[str, Path]) -> None:
    stubbed_gcloud["stdout"].write_text("Pipeline: 42/100 quantified\n")
    deployer = _deployer(tmp_path)
    assert deployer.get_pipeline_status() == "Pipeline: 42/100 quantified"
    deployer.tail_logs(lines=10)
    deployer.get_startup_log()
    calls = stubbed_gcloud["log"].read_text()
    assert "compute ssh metainformant-pipeline --zone us-central1-a --command" in calls
    assert "--tail 10 metainformant-pipeline" in calls
    assert "google-startup-scripts.service" in calls


def test_sync_to_gcs_requires_bucket(tmp_path: Path, stubbed_gcloud: dict[str, Path]) -> None:
    deployer = _deployer(tmp_path, gcs_bucket="")
    assert deployer.sync_to_gcs() is False
    deployer = _deployer(tmp_path, gcs_bucket="my-bucket")
    assert deployer.sync_to_gcs() is True
    assert "gs://my-bucket/amalgkit/" in stubbed_gcloud["log"].read_text()
