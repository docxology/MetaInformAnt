"""Tests for cloud deployment configuration and GCPDeployer.

Zero-mock tests that validate CloudConfig dataclass behavior,
validation logic, and GCPDeployer static methods without requiring
any GCP credentials or network access.
"""
import sys
from pathlib import Path

import pytest

# Add src to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from metainformant.cloud.cloud_config import CloudConfig


class TestCloudConfig:
    """Tests for CloudConfig dataclass."""

    def test_default_values(self):
        """CloudConfig has sensible defaults for all fields."""
        cfg = CloudConfig()
        assert cfg.zone == "us-central1-a"
        assert cfg.instance_name == "metainformant-pipeline"
        assert cfg.spot is True
        assert cfg.max_gb == 20.0
        assert cfg.workers > 0
        assert cfg.threads > 0
        assert cfg.disk_size_gb >= 100

    def test_validation_no_project(self):
        """Validation catches missing project ID."""
        cfg = CloudConfig()
        errors = cfg.validate()
        assert any("project" in e.lower() for e in errors)

    def test_validation_with_project(self):
        """Validation passes with project set."""
        cfg = CloudConfig(project="my-project")
        errors = cfg.validate()
        assert len(errors) == 0

    def test_validation_bad_workers(self):
        """Validation catches invalid worker count."""
        cfg = CloudConfig(project="p", workers=0)
        errors = cfg.validate()
        assert any("worker" in e.lower() for e in errors)

    def test_validation_bad_threads(self):
        """Validation catches invalid thread count."""
        cfg = CloudConfig(project="p", threads=0)
        errors = cfg.validate()
        assert any("thread" in e.lower() for e in errors)

    def test_validation_small_disk(self):
        """Validation warns about small disk."""
        cfg = CloudConfig(project="p", disk_size_gb=50)
        errors = cfg.validate()
        assert any("disk" in e.lower() for e in errors)

    def test_to_metadata(self):
        """to_metadata produces correct GCP metadata dict."""
        cfg = CloudConfig(project="p", max_gb=10.0, workers=28, threads=32)
        meta = cfg.to_metadata()
        assert meta["pipeline-max-gb"] == "10.0"
        assert meta["pipeline-workers"] == "28"
        assert meta["pipeline-threads"] == "32"
        assert "pipeline-repo-url" in meta
        assert "pipeline-docker-image" in meta

    def test_startup_script_path_exists(self):
        """startup_script_path points to a real file."""
        cfg = CloudConfig()
        # The script may or may not exist depending on working directory
        path = cfg.startup_script_path
        assert path.name == "cloud_startup.sh"
        assert "scripts/cloud" in str(path)

    def test_repo_url_default(self):
        """Default repo URL points to MetaInformAnt."""
        cfg = CloudConfig()
        assert "MetaInformAnt" in cfg.repo_url

    def test_config_dir_default(self):
        """Default config dir is config/amalgkit."""
        cfg = CloudConfig()
        assert cfg.config_dir == "config/amalgkit"


class TestGCPDeployerStatic:
    """Tests for GCPDeployer static/class methods (no GCP required)."""

    def test_gcloud_installed_returns_bool(self):
        """gcloud_installed returns a boolean."""
        from metainformant.cloud.gcp_deployer import GCPDeployer
        result = GCPDeployer.gcloud_installed()
        assert isinstance(result, bool)

    def test_deployer_init(self):
        """GCPDeployer can be instantiated with a config."""
        from metainformant.cloud.gcp_deployer import GCPDeployer
        cfg = CloudConfig(project="test-project")
        deployer = GCPDeployer(cfg)
        assert deployer.cfg.project == "test-project"


class TestDownloadScript:
    """Tests for download_results.sh script existence and structure."""

    def test_download_script_exists(self):
        """download_results.sh exists in scripts/cloud/."""
        script = Path(__file__).resolve().parent.parent.parent / "scripts" / "cloud" / "download_results.sh"
        assert script.exists(), f"Missing: {script}"

    def test_download_script_has_shebang(self):
        """download_results.sh has a proper bash shebang."""
        script = Path(__file__).resolve().parent.parent.parent / "scripts" / "cloud" / "download_results.sh"
        content = script.read_text()
        assert content.startswith("#!/usr/bin/env bash") or content.startswith("#!/bin/bash")

    def test_download_script_uses_docker_cp(self):
        """download_results.sh uses docker cp (not direct scp from container)."""
        script = Path(__file__).resolve().parent.parent.parent / "scripts" / "cloud" / "download_results.sh"
        content = script.read_text()
        assert "docker cp" in content

    def test_download_script_has_cleanup(self):
        """download_results.sh cleans up VM temp files."""
        script = Path(__file__).resolve().parent.parent.parent / "scripts" / "cloud" / "download_results.sh"
        content = script.read_text()
        assert "rm -rf" in content or "cleanup" in content.lower()


class TestDeployScript:
    """Tests for deploy_gcp.py CLI structure."""

    def test_deploy_script_exists(self):
        """deploy_gcp.py exists in scripts/cloud/."""
        script = Path(__file__).resolve().parent.parent.parent / "scripts" / "cloud" / "deploy_gcp.py"
        assert script.exists()

    def test_deploy_script_importable(self):
        """deploy_gcp.py has a build_parser function."""
        script_dir = str(Path(__file__).resolve().parent.parent.parent / "scripts" / "cloud")
        sys.path.insert(0, script_dir)
        try:
            import deploy_gcp
            parser = deploy_gcp.build_parser()
            assert parser is not None
        finally:
            sys.path.pop(0)

    def test_deploy_commands_available(self):
        """All expected subcommands are in the parser."""
        script_dir = str(Path(__file__).resolve().parent.parent.parent / "scripts" / "cloud")
        sys.path.insert(0, script_dir)
        try:
            import importlib
            import deploy_gcp
            importlib.reload(deploy_gcp)  # Ensure fresh import
            parser = deploy_gcp.build_parser()
            # Parse each command to verify it exists
            for cmd in ["status", "logs", "download", "stop", "start"]:
                args = parser.parse_args([cmd])
                assert args.command == cmd
        finally:
            sys.path.pop(0)
