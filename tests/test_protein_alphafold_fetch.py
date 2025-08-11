from __future__ import annotations

from pathlib import Path


def _check_online(url: str) -> bool:
    try:
        import requests
        requests.get(url, timeout=5)
        return True
    except Exception:
        return False


def test_alphafold_fetch_builds_and_writes(tmp_path: Path):
    from metainformant.protein.alphafold import build_alphafold_url, fetch_alphafold_model

    url = build_alphafold_url("P69905", version=4, fmt="pdb")
    assert url.endswith("AF-P69905-F1-model_v4.pdb")

    if not _check_online("https://alphafold.ebi.ac.uk"):  # pragma: no cover - skip offline
        import pytest
        pytest.skip("No network access for AlphaFold")

    out = fetch_alphafold_model("P69905", tmp_path, version=4, fmt="pdb")
    assert out.exists() and out.suffix == ".pdb"
    assert out.stat().st_size > 0


