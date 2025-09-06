from __future__ import annotations

from pathlib import Path
from typing import List


def read_pdb_ca_coordinates(pdb_path: Path) -> List[tuple[float, float, float]]:
    """Parse a PDB file and return list of CA atom coordinates as (x,y,z).

    Minimal PDB parser sufficient for unit tests; ignores altloc/occupancy.
    """
    coords: List[tuple[float, float, float]] = []
    for line in Path(pdb_path).read_text().splitlines():
        if not line.startswith("ATOM") and not line.startswith("HETATM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except Exception:
            continue
        coords.append((x, y, z))
    return coords
