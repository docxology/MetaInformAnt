from __future__ import annotations

from pathlib import Path

import numpy as np


def test_read_pdb_ca_and_rmsd(tmp_path: Path):
    from metainformant.protein.structure.general import compute_rmsd_kabsch
    from metainformant.protein.structure.io import read_pdb_ca_coordinates

    pdb_text = (
        "ATOM      1  N   ALA A   1      11.104  13.207   9.479  1.00 20.00           N\n"
        "ATOM      2  CA  ALA A   1      12.560  13.414   9.567  1.00 20.00           C\n"
        "ATOM      3  C   ALA A   1      13.026  13.969  10.935  1.00 20.00           C\n"
        "ATOM      4  N   GLY A   2      14.104  14.207  11.479  1.00 20.00           N\n"
        "ATOM      5  CA  GLY A   2      15.560  14.414  11.567  1.00 20.00           C\n"
    )
    p = tmp_path / "toy.pdb"
    p.write_text(pdb_text)
    ca = read_pdb_ca_coordinates(p)
    assert len(ca) == 2

    A = np.array(ca)
    # Translate and rotate B
    theta = np.pi / 6
    Rz = np.array([[np.cos(theta), -np.sin(theta), 0.0], [np.sin(theta), np.cos(theta), 0.0], [0.0, 0.0, 1.0]])
    B = (A @ Rz) + np.array([2.0, -1.0, 0.5])
    rmsd = compute_rmsd_kabsch(A, B)
    assert rmsd < 1e-6
