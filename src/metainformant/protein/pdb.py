from __future__ import annotations

from pathlib import Path

import requests


def fetch_pdb_structure(pdb_id: str, out_dir: Path, *, fmt: str = "pdb") -> Path:
    """Download a PDB structure by ID in the chosen format (pdb|cif).

    Returns the path to the downloaded file.
    """
    pdb_id_clean = pdb_id.strip().lower()
    out_dir.mkdir(parents=True, exist_ok=True)
    suffix = ".pdb" if fmt == "pdb" else ".cif"
    out_path = out_dir / f"{pdb_id_clean}{suffix}"

    if fmt == "pdb":
        url = f"https://files.rcsb.org/download/{pdb_id_clean}.pdb"
    else:
        url = f"https://files.rcsb.org/download/{pdb_id_clean}.cif"

    resp = requests.get(url, stream=True, timeout=60)
    resp.raise_for_status()
    with open(out_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=8192):
            if chunk:
                fh.write(chunk)
    return out_path


