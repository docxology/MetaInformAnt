"""Tests for DNA restriction site finding."""
from metainformant.dna import restriction


def test_find_restriction_sites_simple() -> None:
    """Test finding HindIII restriction sites."""
    seq = "AAGCTTAAGCTT"  # HindIII site AAGCTT occurs twice
    sites = restriction.find_restriction_sites(seq, {"HindIII": "AAGCTT"})
    assert sites["HindIII"] == [0, 6]


def test_find_restriction_sites_with_enzyme_dict() -> None:
    """Test finding restriction sites returns dict."""
    seq = "GAATTCGAATTC"  # EcoRI site GAATTC occurs twice
    sites = restriction.find_restriction_sites(seq, {"EcoRI": "GAATTC"})
    assert isinstance(sites, dict)
    assert "EcoRI" in sites
    assert len(sites["EcoRI"]) >= 2


def test_find_restriction_sites_no_match() -> None:
    """Test when no restriction sites are found."""
    seq = "AAAAAAAAAAA"  # No HindIII sites
    sites = restriction.find_restriction_sites(seq, {"HindIII": "AAGCTT"})
    assert sites["HindIII"] == []
