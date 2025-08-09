from __future__ import annotations

import subprocess
from pathlib import Path


class DummyCompleted:
    def __init__(self, returncode: int = 0, stdout: str = "ok", stderr: str = "warn"):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


def test_run_amalgkit_writes_logs(monkeypatch, tmp_path: Path):
    from metainformant.rna.amalgkit import run_amalgkit

    # monkeypatch subprocess.run to avoid external call
    def fake_run(cmd, cwd=None, env=None, capture_output=True, text=True, check=False):
        assert isinstance(cmd, list) and cmd[:2] == ["amalgkit", "metadata"]
        return DummyCompleted(0, stdout="STDOUT", stderr="STDERR")

    monkeypatch.setattr(subprocess, "run", fake_run)

    logs = tmp_path / "logs"
    res = run_amalgkit(
        "metadata",
        {"threads": 1},
        work_dir=tmp_path / "work",
        log_dir=logs,
        step_name="metadata",
        check=False,
    )
    assert res.returncode == 0
    # Two log files written
    written = list(logs.glob("*.metadata.*.log"))
    assert any(p.name.endswith("stdout.log") for p in written)
    assert any(p.name.endswith("stderr.log") for p in written)
    # Content matches
    for p in written:
        content = p.read_text()
        assert content in ("STDOUT", "STDERR")


