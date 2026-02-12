"""Protein structure analysis utilities.

This module provides tools for analyzing protein 3D structures,
including RMSD calculations, structural alignments, and geometric analysis.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def compute_rmsd_kabsch(coords_ref: np.ndarray, coords_mobile: np.ndarray) -> float:
    """Calculate RMSD between two protein structures using Kabsch algorithm.

    The Kabsch algorithm finds the optimal rotation and translation to
    minimize RMSD between two sets of coordinates.

    Args:
        coords_ref: Reference coordinates (n_atoms, 3)
        coords_mobile: Mobile coordinates to align (n_atoms, 3)

    Returns:
        RMSD value in Angstroms

    Raises:
        ValueError: If coordinate arrays have incompatible shapes

    Example:
        >>> ref_coords = np.random.rand(100, 3)
        >>> mobile_coords = np.random.rand(100, 3)
        >>> rmsd = compute_rmsd_kabsch(ref_coords, mobile_coords)
        >>> rmsd >= 0
        True
    """
    if coords_ref.shape != coords_mobile.shape:
        raise ValueError("Coordinate arrays must have the same shape")

    if coords_ref.shape[1] != 3:
        raise ValueError("Coordinates must be 3D")

    # Center both structures
    ref_centroid = np.mean(coords_ref, axis=0)
    mobile_centroid = np.mean(coords_mobile, axis=0)

    ref_centered = coords_ref - ref_centroid
    mobile_centered = coords_mobile - mobile_centroid

    # Compute covariance matrix
    covariance = np.dot(mobile_centered.T, ref_centered)

    # SVD
    U, s, Vt = np.linalg.svd(covariance)

    # Check for reflection
    if np.linalg.det(np.dot(U, Vt)) < 0:
        U[:, -1] *= -1

    # Optimal rotation matrix R = U @ Vt (minimizes ||ref - mobile @ R||)
    rotation = np.dot(U, Vt)

    # Apply rotation: aligned = mobile @ R
    aligned_mobile = np.dot(mobile_centered, rotation)

    # Calculate RMSD
    diff = ref_centered - aligned_mobile
    rmsd = np.sqrt(np.sum(diff**2) / len(diff))

    return rmsd


def compute_rmsd_simple(coords_ref: np.ndarray, coords_mobile: np.ndarray) -> float:
    """Calculate simple RMSD without structural alignment.

    Args:
        coords_ref: Reference coordinates
        coords_mobile: Mobile coordinates

    Returns:
        RMSD value

    Example:
        >>> ref_coords = np.random.rand(50, 3)
        >>> mobile_coords = np.random.rand(50, 3)
        >>> rmsd = compute_rmsd_simple(ref_coords, mobile_coords)
        >>> rmsd >= 0
        True
    """
    if coords_ref.shape != coords_mobile.shape:
        raise ValueError("Coordinate arrays must have the same shape")

    diff = coords_ref - coords_mobile
    rmsd = np.sqrt(np.sum(diff**2) / len(diff))

    return rmsd


def align_structures_kabsch(coords_ref: np.ndarray, coords_mobile: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """Align two protein structures using Kabsch algorithm.

    Args:
        coords_ref: Reference coordinates
        coords_mobile: Mobile coordinates to align

    Returns:
        Tuple of (aligned_coords, rotation_matrix, rmsd)

    Example:
        >>> ref = np.random.rand(100, 3)
        >>> mobile = np.random.rand(100, 3)
        >>> aligned, rotation, rmsd = align_structures_kabsch(ref, mobile)
        >>> aligned.shape == ref.shape
        True
    """
    if coords_ref.shape != coords_mobile.shape:
        raise ValueError("Coordinate arrays must have the same shape")

    # Center both structures
    ref_centroid = np.mean(coords_ref, axis=0)
    mobile_centroid = np.mean(coords_mobile, axis=0)

    ref_centered = coords_ref - ref_centroid
    mobile_centered = coords_mobile - mobile_centroid

    # Compute covariance matrix
    covariance = np.dot(mobile_centered.T, ref_centered)

    # SVD
    U, s, Vt = np.linalg.svd(covariance)

    # Check for reflection
    if np.linalg.det(np.dot(U, Vt)) < 0:
        U[:, -1] *= -1

    # Optimal rotation matrix R = U @ Vt
    rotation = np.dot(U, Vt)

    # Apply rotation: aligned = mobile @ R
    aligned_mobile = np.dot(mobile_centered, rotation)

    # Translate back to reference centroid
    aligned_mobile += ref_centroid

    # Calculate RMSD
    diff = ref_centered - (aligned_mobile - ref_centroid)
    rmsd = np.sqrt(np.sum(diff**2) / len(diff))

    return aligned_mobile, rotation, rmsd


def calculate_radius_of_gyration(coords: np.ndarray) -> float:
    """Calculate radius of gyration for a protein structure.

    The radius of gyration measures the compactness of the structure.

    Args:
        coords: Atomic coordinates (n_atoms, 3)

    Returns:
        Radius of gyration in Angstroms

    Example:
        >>> coords = np.random.rand(100, 3) * 50
        >>> rg = calculate_radius_of_gyration(coords)
        >>> rg >= 0
        True
    """
    centroid = np.mean(coords, axis=0)
    distances = np.linalg.norm(coords - centroid, axis=1)
    rg = np.sqrt(np.mean(distances**2))

    return rg


def calculate_center_of_mass(coords: np.ndarray, masses: Optional[np.ndarray] = None) -> np.ndarray:
    """Calculate center of mass for a protein structure.

    Args:
        coords: Atomic coordinates (n_atoms, 3)
        masses: Atomic masses (optional, defaults to equal masses)

    Returns:
        Center of mass coordinates (3,)

    Example:
        >>> coords = np.random.rand(50, 3)
        >>> com = calculate_center_of_mass(coords)
        >>> com.shape
        (3,)
    """
    if masses is None:
        masses = np.ones(len(coords))

    if len(masses) != len(coords):
        raise ValueError("Masses array must match coordinates")

    total_mass = np.sum(masses)
    if total_mass == 0:
        return np.mean(coords, axis=0)

    weighted_coords = coords * masses[:, np.newaxis]
    com = np.sum(weighted_coords, axis=0) / total_mass

    return com


def calculate_inertia_tensor(coords: np.ndarray, masses: Optional[np.ndarray] = None) -> np.ndarray:
    """Calculate inertia tensor for a protein structure.

    Args:
        coords: Atomic coordinates (n_atoms, 3)
        masses: Atomic masses (optional)

    Returns:
        3x3 inertia tensor matrix

    Example:
        >>> coords = np.random.rand(50, 3)
        >>> tensor = calculate_inertia_tensor(coords)
        >>> tensor.shape
        (3, 3)
    """
    if masses is None:
        masses = np.ones(len(coords))

    if len(masses) != len(coords):
        raise ValueError("Masses array must match coordinates")

    # Center coordinates
    com = calculate_center_of_mass(coords, masses)
    coords_centered = coords - com

    # Calculate inertia tensor (vectorized)
    x = coords_centered[:, 0]
    y = coords_centered[:, 1]
    z = coords_centered[:, 2]

    I = np.zeros((3, 3))
    I[0, 0] = np.sum(masses * (y**2 + z**2))
    I[1, 1] = np.sum(masses * (x**2 + z**2))
    I[2, 2] = np.sum(masses * (x**2 + y**2))
    I[0, 1] = I[1, 0] = -np.sum(masses * x * y)
    I[0, 2] = I[2, 0] = -np.sum(masses * x * z)
    I[1, 2] = I[2, 1] = -np.sum(masses * y * z)

    return I


def find_principal_axes(coords: np.ndarray, masses: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Find principal axes of a protein structure.

    Args:
        coords: Atomic coordinates
        masses: Atomic masses (optional)

    Returns:
        Tuple of (eigenvalues, eigenvectors)

    Example:
        >>> coords = np.random.rand(50, 3)
        >>> eigenvals, eigenvecs = find_principal_axes(coords)
        >>> len(eigenvals) == 3
        True
    """
    I = calculate_inertia_tensor(coords, masses)

    # Diagonalize
    eigenvals, eigenvecs = np.linalg.eigh(I)

    # Sort by eigenvalue (largest first)
    idx = np.argsort(eigenvals)[::-1]
    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]

    return eigenvals, eigenvecs


def calculate_structural_statistics(coords: np.ndarray) -> Dict[str, float]:
    """Calculate comprehensive structural statistics.

    Args:
        coords: Atomic coordinates

    Returns:
        Dictionary with structural statistics

    Example:
        >>> coords = np.random.rand(100, 3) * 50
        >>> stats = calculate_structural_statistics(coords)
        >>> "radius_of_gyration" in stats
        True
    """
    stats = {}

    # Basic properties
    stats["n_atoms"] = len(coords)
    stats["center_of_mass"] = calculate_center_of_mass(coords).tolist()

    # Size measures
    stats["radius_of_gyration"] = calculate_radius_of_gyration(coords)

    # Bounding box
    min_coords = np.min(coords, axis=0)
    max_coords = np.max(coords, axis=0)
    stats["bounding_box"] = {
        "min": min_coords.tolist(),
        "max": max_coords.tolist(),
        "dimensions": (max_coords - min_coords).tolist(),
    }

    # Principal axes
    eigenvals, eigenvecs = find_principal_axes(coords)
    stats["principal_axes"] = {"eigenvalues": eigenvals.tolist(), "eigenvectors": eigenvecs.tolist()}

    return stats


def identify_secondary_structure_elements(coords: np.ndarray, backbone_atoms: List[int]) -> List[Dict[str, Any]]:
    """Identify secondary structure elements from coordinates.

    Args:
        coords: All atomic coordinates
        backbone_atoms: Indices of backbone atoms

    Returns:
        List of secondary structure elements

    Example:
        >>> coords = np.random.rand(200, 3)
        >>> backbone = list(range(0, 200, 4))  # Every 4th atom
        >>> elements = identify_secondary_structure_elements(coords, backbone)
        >>> isinstance(elements, list)
        True
    """
    # This is a simplified implementation
    # Real implementation would use DSSP or similar algorithms

    backbone_coords = coords[backbone_atoms]

    elements = []

    # Simple helix detection based on distance patterns
    for i in range(len(backbone_coords) - 4):
        # Check for alpha helix pattern (3.8 Å rise, 100° rotation)
        segment = backbone_coords[i : i + 5]

        if len(segment) >= 4:
            # Calculate distances between every 4th residue
            dist1 = np.linalg.norm(segment[0] - segment[3])
            dist2 = np.linalg.norm(segment[1] - segment[4]) if len(segment) > 4 else 0

            if 5.0 < dist1 < 7.0:  # Approximate helix distance
                elements.append({"type": "helix", "start_residue": i, "end_residue": i + 4, "confidence": 0.8})

    return elements


def calculate_solvent_accessible_surface_area(coords: np.ndarray, probe_radius: float = 1.4) -> float:
    """Calculate solvent accessible surface area.

    Args:
        coords: Atomic coordinates
        probe_radius: Radius of solvent probe (default: 1.4 Å for water)

    Returns:
        Surface area in Å²

    Example:
        >>> coords = np.random.rand(50, 3) * 20
        >>> area = calculate_solvent_accessible_surface_area(coords)
        >>> area >= 0
        True
    """
    # This is a simplified implementation
    # Real SASA calculation requires more complex algorithms

    # Simple approximation using convex hull
    try:
        from scipy.spatial import ConvexHull

        hull = ConvexHull(coords)
        # Approximate surface area from convex hull
        surface_area = hull.area
    except ImportError:
        # Fallback: approximate as sphere
        radius = np.max(np.linalg.norm(coords - np.mean(coords, axis=0), axis=1))
        surface_area = 4 * np.pi * (radius + probe_radius) ** 2

    return surface_area
