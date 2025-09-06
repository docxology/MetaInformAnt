from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_protein_cli_rmsd_ca(tmp_path: Path):
    pdb1 = tmp_path / "a.pdb"
    pdb2 = tmp_path / "b.pdb"
    pdb1.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      2  CA  GLY A   2       1.000   0.000   0.000  1.00 20.00           C\n"
    )
    pdb2.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      2  CA  GLY A   2       0.000   1.000   0.000  1.00 20.00           C\n"
    )

    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "protein",
        "rmsd-ca",
        "--pdb-a",
        str(pdb1),
        "--pdb-b",
        str(pdb2),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0
    val = float(r.stdout.strip())
    assert 0.0 <= val < 2.0
