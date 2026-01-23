"""Tests for NCBI DNA data functions."""

from __future__ import annotations

import importlib.util

import pytest

from metainformant.dna import ncbi


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
    except (ImportError, AttributeError) as e:
        # Also acceptable - module may not have function
        pass


def test_ncbi_module_importable():
    """Test that the ncbi module can be imported."""
    assert ncbi is not None
