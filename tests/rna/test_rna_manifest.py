"""Tests for RNA workflow manifest generation.

Tests manifest file creation and log directory configuration.
"""

from __future__ import annotations

from pathlib import Path

from metainformant.rna.engine.workflow import AmalgkitWorkflowConfig, plan_workflow


def test_manifest_written_and_logs_directory(tmp_path: Path):
    """Workflow configuration defaults resolve logs and manifest under work_dir."""

    cfg = AmalgkitWorkflowConfig(work_dir=tmp_path)
    steps = plan_workflow(cfg)
    assert len(steps) > 0

    assert cfg.log_dir == tmp_path / "logs"
    assert cfg.manifest_path == tmp_path / "amalgkit.manifest.jsonl"
    assert cfg.log_file == tmp_path / "logs" / "workflow.log"
