"""Tests for preflight manifest generation when amalgkit is missing.

Tests that workflows create preflight manifests documenting missing dependencies.
"""

from __future__ import annotations

from pathlib import Path


def test_preflight_manifest_when_amalgkit_missing(tmp_path: Path):
    """Test that workflow creates preflight manifest when amalgkit CLI is not available."""
    from metainformant.rna.amalgkit import check_cli_available
    from metainformant.rna.workflow import AmalgkitWorkflowConfig, execute_workflow

    ok, _ = check_cli_available()
    if ok:
        import pytest

        pytest.skip("amalgkit present; this test targets missing CLI scenario")

    cfg = AmalgkitWorkflowConfig(work_dir=tmp_path / "work", threads=1)
    codes = execute_workflow(cfg, check=False)
    assert codes == [127]

    manifest = cfg.work_dir / "amalgkit.manifest.jsonl"
    assert manifest.exists()
    text = manifest.read_text()
    assert "amalgkit -h" in text
    assert "preflight" in text
