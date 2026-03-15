from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_protein_cli_taxon_ids(tmp_path: Path):
    path = tmp_path / "taxon.txt"
    path.write_text("# comment\n9606\n10090\nfoo\n7227\n")

    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "protein",
        "taxon-ids",
        "--file",
        str(path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    out = result.stdout.strip().split()
    assert set(out) >= {"9606", "10090", "7227"}
