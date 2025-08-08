from __future__ import annotations

from pathlib import Path

import requests


def build_alphafold_url(uniprot_acc: str, *, version: int = 4, fmt: str = "pdb") -> str:
    acc = uniprot_acc.strip().upper()
    if fmt == "pdb":
        return f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{version}.pdb"
    elif fmt == "cif":
        return f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{version}.cif"
    else:
        raise ValueError("fmt must be 'pdb' or 'cif'")


def fetch_alphafold_model(uniprot_acc: str, out_dir: Path, *, version: int = 4, fmt: str = "pdb") -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    url = build_alphafold_url(uniprot_acc, version=version, fmt=fmt)
    suffix = ".pdb" if fmt == "pdb" else ".cif"
    out_path = out_dir / f"AF-{uniprot_acc.upper()}-F1-model_v{version}{suffix}"
    resp = requests.get(url, stream=True, timeout=60)
    resp.raise_for_status()
    with open(out_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=8192):
            if chunk:
                fh.write(chunk)
    return out_path


