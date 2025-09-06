from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_rna_plan_cli_lists_expected_steps(tmp_path: Path):
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "rna",
        "plan",
        "--work-dir",
        str(tmp_path),
        "--threads",
        "2",
        "--species",
        "Apis_mellifera",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    out = result.stdout
    assert "metadata" in out
    assert "sanity" in out


def test_rna_plan_species_cli_includes_species_and_tissue(tmp_path: Path):
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "rna",
        "plan-species",
        "--work-dir",
        str(tmp_path),
        "--threads",
        "2",
        "--taxon-id",
        "7460",
        "--tissue",
        "brain",
        "--tissue",
        "thorax",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    out = result.stdout
    assert "metadata" in out and "select" in out
    assert "taxon-id" in out
    assert "brain" in out and "thorax" in out
