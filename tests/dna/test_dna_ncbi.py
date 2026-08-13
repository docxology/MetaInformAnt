"""Tests for NCBI DNA data functions."""

from __future__ import annotations

import importlib.util

import pytest

from metainformant.core.ncbi import NCBIContactError
from metainformant.dna.external import ncbi


def _check_ncbi_datasets_installed() -> bool:
    """Check if ncbi-datasets-pylib is installed."""
    try:
        return importlib.util.find_spec("ncbi.datasets") is not None
    except ModuleNotFoundError:
        return False


@pytest.mark.skipif(
    _check_ncbi_datasets_installed(), reason="ncbi-datasets installed; test only for missing dependency case"
)
def test_ncbi_datasets_optional_dependency_errors():
    """Test that clear errors are raised when ncbi-datasets is not installed."""
    # This test only runs when ncbi-datasets is NOT installed
    # Expect a clear error if ncbi-datasets-pylib is not installed
    try:
        ncbi.get_accession_by_tax_id("9606")
        # If we get here, either: function works, or different error type
        pytest.skip("ncbi function available, skipping missing dependency test")
    except RuntimeError:
        pass  # Expected
    except (ImportError, AttributeError):
        # Also acceptable - module may not have function
        pass


def test_ncbi_module_importable():
    """Test that the ncbi module can be imported."""
    assert ncbi is not None


def test_ncbi_client_requires_explicit_contact(monkeypatch: pytest.MonkeyPatch):
    """DNA NCBI clients enforce the shared contact policy before networking."""
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    with pytest.raises(NCBIContactError):
        ncbi.NCBIClient()


def test_ncbi_client_supports_explicit_anonymous_opt_in(monkeypatch: pytest.MonkeyPatch):
    """Anonymous mode is explicit and leaves the request identity unset."""
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    with pytest.warns(RuntimeWarning, match="anonymous NCBI"):
        client = ncbi.NCBIClient(allow_anonymous=True)
    assert client.email is None
    assert client.contact_mode == "anonymous"


def test_ncbi_client_rejects_non_positive_timeout(monkeypatch: pytest.MonkeyPatch):
    """External requests cannot be configured with an unbounded timeout."""
    monkeypatch.setenv("NCBI_EMAIL", "test@example.org")
    with pytest.raises(ValueError, match="timeout"):
        ncbi.NCBIClient(timeout=0)
