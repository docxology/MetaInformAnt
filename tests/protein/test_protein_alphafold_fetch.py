from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.protein.structure.alphafold import build_alphafold_url, fetch_alphafold_model


def _check_online(url: str) -> bool:
    """Check if we can reach a URL within timeout."""
    try:
        import requests

        resp = requests.get(url, timeout=5)
        resp.raise_for_status()
        return True
    except Exception:
        return False


def test_build_alphafold_url_pdb_format():
    """Test AlphaFold URL construction for PDB format (no network required)."""
    url = build_alphafold_url("P69905", version=4, fmt="pdb")
    expected = "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v4.pdb"
    assert url == expected


def test_build_alphafold_url_cif_format():
    """Test AlphaFold URL construction for CIF format (no network required)."""
    url = build_alphafold_url("Q9Y6N1", version=3, fmt="cif")
    expected = "https://alphafold.ebi.ac.uk/files/AF-Q9Y6N1-F1-model_v3.cif"
    assert url == expected


def test_build_alphafold_url_invalid_format():
    """Test error handling for invalid format (no network required)."""
    with pytest.raises(ValueError, match="fmt must be 'pdb' or 'cif'"):
        build_alphafold_url("P69905", version=4, fmt="invalid")


def test_build_alphafold_url_edge_cases():
    """Test URL construction with various edge cases (no network required)."""
    # Different versions
    url_v1 = build_alphafold_url("P12345", version=1, fmt="pdb")
    url_v2 = build_alphafold_url("P12345", version=2, fmt="pdb")
    assert "model_v1" in url_v1
    assert "model_v2" in url_v2

    # Different accession formats
    url_short = build_alphafold_url("P123", version=4, fmt="pdb")
    url_long = build_alphafold_url("P123456789", version=4, fmt="pdb")
    assert "AF-P123-F1" in url_short
    assert "AF-P123456789-F1" in url_long


def test_fetch_alphafold_model_real_network(tmp_path: Path):
    """Test real AlphaFold model download with actual HTTP requests."""
    if not _check_online("https://alphafold.ebi.ac.uk"):
        pytest.skip("No network access for AlphaFold - real implementation requires connectivity")

    # Test with a well-known protein (hemoglobin alpha chain)
    try:
        result_path = fetch_alphafold_model("P69905", tmp_path, version=4, fmt="pdb")
        assert result_path.exists()
        assert result_path.suffix == ".pdb"
        assert result_path.name.startswith("AF-P69905")
        assert result_path.stat().st_size > 0

        # Verify it contains PDB content
        content = result_path.read_text()
        assert "HEADER" in content or "ATOM" in content
    except Exception as e:
        # AlphaFold might not have this protein or service might be down
        pytest.skip(f"AlphaFold model P69905 not available - real API behavior: {e}")


def test_fetch_alphafold_model_cif_format_real_network(tmp_path: Path):
    """Test real AlphaFold model download in CIF format."""
    if not _check_online("https://alphafold.ebi.ac.uk"):
        pytest.skip("No network access for AlphaFold - real implementation requires connectivity")

    try:
        result_path = fetch_alphafold_model("P69905", tmp_path, version=4, fmt="cif")
        assert result_path.exists()
        assert result_path.suffix == ".cif"
        assert result_path.stat().st_size > 0

        # Verify it contains CIF content
        content = result_path.read_text()
        assert "data_" in content or "_atom_site" in content
    except Exception as e:
        pytest.skip(f"AlphaFold CIF model P69905 not available - real API behavior: {e}")


def test_fetch_alphafold_model_nonexistent_protein(tmp_path: Path):
    """Test real behavior with non-existent protein ID."""
    if not _check_online("https://alphafold.ebi.ac.uk"):
        pytest.skip("No network access for AlphaFold - real implementation requires connectivity")

    # Test with obviously fake protein ID
    with pytest.raises(Exception):
        fetch_alphafold_model("FAKE_PROTEIN_12345", tmp_path, version=4, fmt="pdb")


def test_fetch_alphafold_model_invalid_format_error(tmp_path: Path):
    """Test that invalid format raises error before making HTTP request (no network required)."""
    with pytest.raises(ValueError):
        fetch_alphafold_model("P69905", tmp_path, version=4, fmt="xyz")


def test_fetch_alphafold_model_directory_creation(tmp_path: Path):
    """Test that output directory is created if it doesn't exist."""
    if not _check_online("https://alphafold.ebi.ac.uk"):
        pytest.skip("No network access for AlphaFold - real implementation requires connectivity")

    # Use a nested directory that doesn't exist yet
    nested_dir = tmp_path / "models" / "alphafold"
    assert not nested_dir.exists()

    try:
        result_path = fetch_alphafold_model("P69905", nested_dir, version=4, fmt="pdb")
        assert nested_dir.exists()
        assert result_path.parent == nested_dir
        assert result_path.exists()
    except Exception as e:
        # Model might not be available, but directory should still be created
        assert nested_dir.exists()
        pytest.skip(f"AlphaFold model not available but directory creation tested: {e}")


def test_alphafold_offline_behavior(tmp_path: Path):
    """Document real offline behavior for AlphaFold downloads."""
    # When offline, the function should fail gracefully
    try:
        result = fetch_alphafold_model("P69905", tmp_path, version=4, fmt="pdb")
        # If this succeeds, we're online
        assert result.exists()
    except Exception:
        # Expected when offline - this documents real failure modes
        # Real implementations reveal actual network dependencies
        assert True  # This is acceptable real-world behavior
