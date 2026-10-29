"""Invalid CLI options must refuse before importing the scientific stack."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from tests.rna.test_ehd_runtime_paths import require_operator_files

pytestmark = pytest.mark.external_tool


@pytest.mark.parametrize("value", ["0", "9", "invalid"])
def test_invalid_options_refuse_before_scientific_imports(tmp_path: Path, value: str) -> None:
    source = require_operator_files("continuous.py", "campaign.py", "checkpoints.py")
    private = tmp_path / "operations"
    private.mkdir()
    for name in ("continuous.py", "campaign.py", "checkpoints.py"):
        shutil.copy2(source / name, private / name)
    command = [
        sys.executable,
        "-v",
        str(private / "continuous.py"),
        "--quant-slots",
        value,
        "--data-root",
        str(tmp_path / "data"),
        "--config-dir",
        str(tmp_path / "configs"),
    ]
    environment = dict(os.environ)
    result = subprocess.run(command, env=environment, capture_output=True, text=True, timeout=20, check=False)
    (tmp_path / "import-trace.log").write_text(result.stderr)
    assert result.returncode == 2 and "quant" in result.stderr
    assert "import 'metainformant" not in result.stderr
    assert not (private / ".run.lock").exists()
    assert not (private / "active_runtime.json").exists()
