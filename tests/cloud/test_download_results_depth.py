"""Depth tests for GCPDeployer.download_results.

Two real-path branches:
- Script branch: the deployer shells out to the REAL
  scripts/cloud/download_results.sh with --output; we verify the call
  contract (script accepts --output, deployer creates the local dir first).
- Fallback branch: with the script moved aside, the deployer falls back to
  direct `gcloud compute scp` calls, exercised via a real gcloud stub on
  PATH (same pattern as tests/gwas/test_gwas_sra_environment.py).
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from metainformant.cloud.cloud_config import CloudConfig
from metainformant.cloud.gcp_deployer import GCPDeployer

GCLOUD_STUB = """#!/bin/sh
printf '%s\\n' "$*" >> "$GCLOUD_CALL_LOG"
exit 0
"""


@pytest.fixture()
def stubbed_gcloud(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    tool = bin_dir / "gcloud"
    tool.write_text(GCLOUD_STUB, encoding="utf-8")
    tool.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    call_log = tmp_path / "calls.log"
    call_log.touch()
    monkeypatch.setenv("GCLOUD_CALL_LOG", str(call_log))
    return call_log


def _deployer() -> GCPDeployer:
    return GCPDeployer(CloudConfig(project="proj"))


def _repo_script() -> Path:
    return Path(__file__).resolve().parents[2] / "scripts" / "cloud" / "download_results.sh"


def test_download_results_script_branch_contract(tmp_path: Path) -> None:
    """The real repo script accepts --output DIR; the deployer creates the
    local directory before shelling out."""
    deployer = _deployer()
    assert _repo_script().exists()  # branch precondition

    local_dir = tmp_path / "results"
    proc = subprocess.run(
        ["bash", str(_repo_script()), "--output", str(local_dir)],
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
        env={**os.environ, "GCP_PROJECT": "proj", "GCP_ZONE": "us-central1-a"},
    )
    # The script accepts --output (an unknown-arg error would say so); it may
    # fail later without docker/gcloud access in CI, which is fine here.
    assert "Unknown argument" not in proc.stderr
    # Deployer contract: the local directory is created before shelling out.
    # The script itself needs docker/gcloud access, so the deployer call may
    # either succeed (tooling present) or surface the failure via
    # CalledProcessError (check=True) - both preserve the directory contract.
    contract_dir = tmp_path / "contract_dir"
    try:
        result = deployer.download_results(str(contract_dir))
        assert isinstance(result, bool)
    except subprocess.CalledProcessError:
        pass
    assert contract_dir.exists()


def test_download_results_fallback_branch_uses_scp(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, stubbed_gcloud: Path
) -> None:
    """With the repo script absent, deployer falls back to direct gcloud scp."""
    deployer = _deployer()
    real_script = _repo_script()
    hidden = real_script.with_suffix(".sh.hidden-by-test")
    shutil.move(str(real_script), str(hidden))
    try:
        local_dir = tmp_path / "fallback_results"
        result = deployer.download_results(str(local_dir))
        assert result is True
        assert local_dir.exists()
    finally:
        shutil.move(str(hidden), str(real_script))
    assert real_script.exists()  # restored

    calls = stubbed_gcloud.read_text()
    assert "compute scp --recurse" in calls
    assert "--zone us-central1-a" in calls
    assert "--project proj" in calls
    assert "/opt/MetaInformAnt/projects/hymenoptera_amalgkit/data/*/work/quant" in calls
    assert "pipeline_progress.db" in calls
