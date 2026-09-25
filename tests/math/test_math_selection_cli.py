from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_selection_replay_removed(tmp_path: Path) -> None:
    """The fabricated ``math selection replay`` figure generator is removed.

    Invoking it must fail as an invalid CLI choice and write no figures;
    this guards against re-adding a command that fabricated PNG outputs
    without a real implementation behind them.
    """
    dest = tmp_path / "selection_experiments"
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "math",
        "selection",
        "replay",
        "--dest",
        str(dest),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0, result.stdout
    assert "invalid choice" in (result.stderr + result.stdout)
    assert not dest.exists()
