"""Tests for the RNA surface of the metainformant CLI."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def _cli_env() -> dict[str, str]:
    """Environment with the repository src on PYTHONPATH."""

    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).parent.parent / "src")
    return env


def test_rna_info_lists_subpackages():
    """'rna info' prints the module summary and exits 0."""

    result = subprocess.run(
        [sys.executable, "-m", "metainformant", "rna", "info"],
        capture_output=True,
        text=True,
        env=_cli_env(),
    )
    assert result.returncode == 0, result.stderr
    assert "RNA-seq Analysis Module" in result.stdout
    assert "amalgkit" in result.stdout
    assert "engine" in result.stdout


def test_rna_plan_is_not_a_supported_subcommand(tmp_path: Path):
    """No 'rna plan' subcommand exists; the CLI must say so, not crash."""

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "metainformant",
            "rna",
            "plan",
            "--work-dir",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
        env=_cli_env(),
    )
    assert result.returncode != 0
    assert "ModuleNotFoundError" not in result.stderr
