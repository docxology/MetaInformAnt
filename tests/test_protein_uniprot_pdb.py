from __future__ import annotations

from pathlib import Path


def _check_online(url: str) -> bool:
    try:
        import requests
        requests.get(url, timeout=5)
        return True
    except Exception:
        return False


def test_uniprot_id_mapping_real_network():
    if not _check_online("https://rest.uniprot.org"):  # pragma: no cover - skip offline
        import pytest
        pytest.skip("No network access for UniProt")
    from metainformant.protein.uniprot import map_ids_uniprot

    result = map_ids_uniprot(["P69905"])  # hemoglobin alpha
    assert isinstance(result, dict)
    assert "P69905" in result


def test_pdb_download_builds_paths_and_writes(tmp_path: Path):
    if not _check_online("https://files.rcsb.org"):  # pragma: no cover - skip offline
        import pytest
        pytest.skip("No network access for PDB")
    from metainformant.protein.pdb import fetch_pdb_structure

    out = fetch_pdb_structure("1CRN", tmp_path, fmt="pdb")
    assert out.exists() and out.suffix == ".pdb"
    assert out.stat().st_size > 0


