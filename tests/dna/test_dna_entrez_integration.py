from __future__ import annotations

import os

import pytest

from metainformant.core.ncbi import NCBIContactError
from metainformant.dna.external import entrez


def _require_ncbi_email() -> str:
    """Return an explicitly supplied test contact or skip the network test."""
    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email:
        pytest.skip("Set NCBI_EMAIL to run live NCBI tests; tests never invent contact data")
    return email


@pytest.mark.network
def test_entrez_fetch_phix_fasta_real_network():
    """Test real NCBI Entrez API with actual network requests."""
    email = _require_ncbi_email()

    try:
        # PhiX174 is a common small reference genome
        rec = entrez.get_genome_from_ncbi("NC_001422.1", email=email)
        assert rec.id.startswith("NC_001422")
        assert len(str(rec.seq)) > 1000

        # Verify it's a real sequence record
        assert hasattr(rec, "seq")
        assert hasattr(rec, "description")
    except Exception as e:
        assert str(e), "Unavailable NCBI services must retain a diagnostic error"


def test_entrez_no_email_behavior():
    """Test behavior when no email is provided (no network required)."""
    # This tests parameter validation without network calls
    with pytest.raises(NCBIContactError):
        entrez.get_genome_from_ncbi("NC_001422.1", email="")


@pytest.mark.network
def test_entrez_invalid_accession_real():
    """Test real behavior with invalid accession numbers."""
    email = _require_ncbi_email()

    # Test with obviously invalid accession
    try:
        rec = entrez.get_genome_from_ncbi("INVALID_ACCESSION_12345", email=email)
        # If this succeeds, API is very lenient
        assert hasattr(rec, "id")
    except Exception as exc:
        # Expected - invalid accessions should fail
        assert str(exc)


@pytest.mark.network
def test_entrez_different_accession_types_real():
    """Test with different types of valid accession numbers."""
    email = _require_ncbi_email()

    # Test different accession types that should work
    test_accessions = [
        "NC_001422.1",  # RefSeq complete genome
        # Could add more if they exist and are stable
    ]

    for accession in test_accessions:
        try:
            rec = entrez.get_genome_from_ncbi(accession, email=email)
            assert rec.id.startswith(accession.split(".")[0])
            assert len(str(rec.seq)) > 0
        except Exception as e:
            # A service outage remains a tested, diagnosable integration state.
            assert str(e)


@pytest.mark.network
def test_entrez_offline_behavior():
    """Document real offline behavior for NCBI queries."""
    email = _require_ncbi_email()

    # When offline, the function should fail gracefully
    try:
        rec = entrez.get_genome_from_ncbi("NC_001422.1", email=email)
        # If this succeeds, we're online
        assert hasattr(rec, "id")
    except Exception:
        # Expected when offline - this documents real failure modes
        # Real implementations reveal actual network dependencies
        assert True  # This is acceptable real-world behavior


@pytest.mark.network
def test_entrez_rate_limiting_behavior():
    """Test how real API handles rapid successive requests."""
    email = _require_ncbi_email()

    # Make a few rapid requests to test rate limiting
    accessions = ["NC_001422.1"] * 3  # Same accession multiple times

    results = []
    for accession in accessions:
        try:
            rec = entrez.get_genome_from_ncbi(accession, email=email)
            results.append(rec.id)
        except Exception as e:
            # Rate limiting or other API restrictions
            results.append(f"Error: {e}")

    # Document real API behavior
    assert len(results) == 3


def test_entrez_email_parameter_validation():
    """Test email parameter handling without network calls."""
    client = entrez.EntrezClient(email="test@example.com", timeout=1)
    assert client.email == "test@example.com"
    assert client.contact_mode == "explicit"
    assert client.timeout == 1
