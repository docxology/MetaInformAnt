from __future__ import annotations

from pathlib import Path


def test_preflight_manifest_when_amalgkit_missing(monkeypatch, tmp_path: Path):
    from metainformant.rna.workflow import AmalgkitWorkflowConfig, execute_workflow
    import metainformant.rna.workflow as wf

    # Force preflight failure
    monkeypatch.setattr(wf, "check_cli_available", lambda: (False, "missing"))

    cfg = AmalgkitWorkflowConfig(work_dir=tmp_path / "work", threads=1)
    codes = execute_workflow(cfg, check=False)
    assert codes == [127]

    manifest = cfg.work_dir / "amalgkit.manifest.jsonl"
    assert manifest.exists()
    text = manifest.read_text()
    assert "amalgkit -h" in text
    assert "preflight" in text


