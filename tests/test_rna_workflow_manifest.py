from pathlib import Path
import json

from metainformant.rna.workflow import AmalgkitWorkflowConfig, execute_workflow
from metainformant.rna.amalgkit import check_cli_available


def test_manifest_written_with_records(tmp_path: Path):
    ok, _ = check_cli_available()
    if not ok:
        import pytest
        pytest.skip("amalgkit not available on PATH")

    cfg = AmalgkitWorkflowConfig(work_dir=tmp_path / "work", threads=1)
    # Provide minimal metadata to allow downstream steps to run if applicable
    (cfg.work_dir / "metadata").mkdir(parents=True, exist_ok=True)
    (cfg.work_dir / "metadata" / "metadata.tsv").write_text("tissue\nbrain\n", encoding="utf-8")

    codes = execute_workflow(cfg, check=False)
    assert len(codes) >= 1

    manifest = cfg.work_dir / "amalgkit.manifest.jsonl"
    assert manifest.exists()
    lines = manifest.read_text().strip().splitlines()
    assert len(lines) >= 1
    rec = json.loads(lines[0])
    assert "step" in rec and "return_code" in rec and "duration_seconds" in rec


