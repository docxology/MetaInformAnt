from pathlib import Path
import json

from metainformant.rna.workflow import AmalgkitWorkflowConfig, execute_workflow


class DummyResult:
    def __init__(self, rc: int, out: str = "ok", err: str = ""):
        self.returncode = rc
        self.stdout = out
        self.stderr = err


def test_manifest_written_with_records(monkeypatch, tmp_path: Path):
    # monkeypatch step runners to avoid calling external amalgkit
    from metainformant.rna import steps as step_mod

    def fake_runner(params, **kwargs):
        return DummyResult(0, out="x", err="")

    monkeypatch.setattr(step_mod, "STEP_RUNNERS", {
        "metadata": fake_runner,
        "integrate": fake_runner,
        "config": fake_runner,
        "select": fake_runner,
        "getfastq": fake_runner,
        "quant": fake_runner,
        "merge": fake_runner,
        "cstmm": fake_runner,
        "curate": fake_runner,
        "csca": fake_runner,
        "sanity": fake_runner,
    })

    # also bypass preflight check
    import metainformant.rna.workflow as wf
    monkeypatch.setattr(wf, "check_cli_available", lambda: (True, "ok"))

    cfg = AmalgkitWorkflowConfig(work_dir=tmp_path / "work", threads=1)
    codes = execute_workflow(cfg, check=True)
    assert all(c == 0 for c in codes)

    manifest = cfg.work_dir / "amalgkit.manifest.jsonl"
    assert manifest.exists()
    lines = manifest.read_text().strip().splitlines()
    assert len(lines) >= 3
    rec = json.loads(lines[0])
    assert "step" in rec and "return_code" in rec and "duration_seconds" in rec


