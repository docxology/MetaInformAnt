from __future__ import annotations

from pathlib import Path


def test_workflow_skips_steps_when_missing_deps(tmp_path: Path):
    from metainformant.rna.amalgkit import check_cli_available
    from metainformant.rna.deps import check_step_dependencies
    from metainformant.rna.workflow import AmalgkitWorkflowConfig, execute_workflow

    ok, _ = check_cli_available()
    if not ok:
        import pytest

        pytest.skip("amalgkit not available on PATH")

    cfg = AmalgkitWorkflowConfig(work_dir=tmp_path / "work", threads=1)
    (cfg.work_dir / "metadata").mkdir(parents=True, exist_ok=True)
    (cfg.work_dir / "metadata" / "metadata.tsv").write_text("tissue\nbrain\n", encoding="utf-8")
    # Execute; tests rely on real dependency checks to naturally skip where needed
    codes = execute_workflow(cfg, check=False)
    # If any step lacks deps, 126 will appear; otherwise we at least ran metadata
    assert any(c in (0, 126, 204) for c in codes)
