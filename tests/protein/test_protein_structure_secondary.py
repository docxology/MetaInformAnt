"""Tests for protein secondary structure prediction."""

from __future__ import annotations

import pytest


def test_secondary_structure_module_exists():
    """Test that the secondary structure module is importable."""
    # Import the main structure module
    from metainformant.protein import structure

    assert structure is not None


def test_general_structure_functions():
    """Test basic structure functions are available."""
    import numpy as np

    from metainformant.protein.structure.general import compute_rmsd_kabsch

    # Simple test
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    rmsd = compute_rmsd_kabsch(coords, coords)
    assert rmsd < 1e-9
