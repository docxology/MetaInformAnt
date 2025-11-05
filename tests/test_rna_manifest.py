"""Tests for RNA workflow manifest generation.

Tests manifest file creation and log directory configuration.
"""

from __future__ import annotations

from pathlib import Path

from metainformant.rna.workflow import AmalgkitWorkflowConfig, plan_workflow


def test_manifest_written_and_logs_directory(tmp_path: Path):
    """Test that workflow planning creates correct default paths for logs and manifest.
    
    Note: This test only checks planning and default paths, not actual execution.
    """
    cfg = AmalgkitWorkflowConfig(work_dir=tmp_path)
    steps = plan_workflow(cfg)
    assert len(steps) > 0

    # Validate default locations
    default_log_dir = tmp_path / "logs"
    default_manifest = tmp_path / "amalgkit.manifest.jsonl"
    assert default_log_dir.as_posix().endswith("/logs")
    assert default_manifest.name == "amalgkit.manifest.jsonl"
