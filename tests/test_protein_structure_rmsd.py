"""Tests for protein structure RMSD calculation."""
from __future__ import annotations

import numpy as np

from metainformant.protein.structure.general import compute_rmsd_kabsch


def test_kabsch_rmsd_zero_for_identical():
    """Test RMSD is zero for identical coordinates."""
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-1.0, 0.5, 2.5]])
    rmsd = compute_rmsd_kabsch(coords, coords.copy())
    assert abs(rmsd - 0.0) < 1e-9


def test_kabsch_rmsd_invariant_to_rotation_translation():
    """Test RMSD is invariant to rotation and translation."""
    A = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
    # Rotate around z by 90 degrees and translate
    theta = np.pi / 2
    Rz = np.array([[np.cos(theta), -np.sin(theta), 0.0], [np.sin(theta), np.cos(theta), 0.0], [0.0, 0.0, 1.0]])
    B = (A @ Rz) + np.array([10.0, -5.0, 2.0])
    rmsd = compute_rmsd_kabsch(A, B)
    assert rmsd < 1e-6
