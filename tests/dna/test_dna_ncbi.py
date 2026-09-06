"""Tests for NCBI DNA data functions."""

from __future__ import annotations

import pytest

from metainformant.core.ncbi import NCBIContactError
from metainformant.dna.external import ncbi


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
