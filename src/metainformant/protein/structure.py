from __future__ import annotations

import numpy as np


def compute_rmsd_kabsch(coords_ref: np.ndarray, coords_mobile: np.ndarray) -> float:
    """Compute RMSD after optimal superposition using the Kabsch algorithm.

    Shapes must be (N, 3) for both arrays. Returns RMSD as float.
    """
    a = np.asarray(coords_ref, dtype=float)
    b = np.asarray(coords_mobile, dtype=float)
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 3:
        raise ValueError("Coordinate arrays must be of shape (N, 3) and equal")
    # Center the coordinates
    a_center = a.mean(axis=0)
    b_center = b.mean(axis=0)
    A = a - a_center
    B = b - b_center
    # Covariance and SVD
    H = A.T @ B
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1.0, 1.0, np.sign(d)])
    R = Vt.T @ D @ U.T
    B_rot = B @ R
    diff = A - B_rot
    rmsd = np.sqrt((diff * diff).sum() / a.shape[0])
    return float(rmsd)
