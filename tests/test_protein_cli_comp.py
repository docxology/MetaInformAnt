from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_protein_cli_comp(tmp_path: Path):
    fasta = ">a\nMKKLL\n>b\nGGGG\n"
    p = tmp_path / "x.faa"
    p.write_text(fasta)

    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "protein",
        "comp",
        "--fasta",
        str(p),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0
    out = r.stdout.strip().splitlines()
    assert out and out[0].startswith("a\t") and "," in out[0]
    assert any(line.startswith("b\tG:") for line in out)


