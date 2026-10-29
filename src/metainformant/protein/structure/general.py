"""Protein structure analysis utilities.

This module provides tools for analyzing protein 3D structures: Kabsch superposition and RMSD, the TM-score,
Shrake-Rupley solvent accessible surface areas, P-SEA secondary structure assignment from a C-alpha trace, and
simple geometric descriptors (radius of gyration, center of mass, inertia tensor).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.spatial import cKDTree

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def _kabsch_rotation(ref_centered: np.ndarray, mobile_centered: np.ndarray) -> np.ndarray:
    """Proper rotation ``R`` (determinant +1) that minimises ``||ref - mobile @ R||`` for centred (n, 3) arrays.

    The rotation comes from the singular value decomposition of the covariance matrix. When ``U @ Vt`` is a
    reflection (determinant -1, which happens for mirror-image point sets) the singular vector of the smallest
    singular value is negated, so the result is always a rotation and never a reflection.
    """
    covariance = np.dot(mobile_centered.T, ref_centered)
    u, _, vt = np.linalg.svd(covariance)
    if np.linalg.det(np.dot(u, vt)) < 0:
        u[:, -1] *= -1
    return np.asarray(np.dot(u, vt))


def _kabsch_transform(coords_ref: np.ndarray, coords_mobile: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Rotation ``R`` and translation ``t`` with ``coords_mobile @ R + t`` superposed on ``coords_ref``."""
    ref_centroid = coords_ref.mean(axis=0)
    mobile_centroid = coords_mobile.mean(axis=0)
    rotation = _kabsch_rotation(coords_ref - ref_centroid, coords_mobile - mobile_centroid)
    return rotation, ref_centroid - np.dot(mobile_centroid, rotation)


def _as_points(coords: Any, name: str = "coords") -> np.ndarray:
    """``coords`` as a finite float array of shape (n, 3)."""
    points = np.asarray(coords, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n, 3), got {points.shape}")
    if not np.isfinite(points).all():
        raise ValueError(f"{name} must contain only finite coordinates")
    return points


def _dihedral_degrees(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> np.ndarray:
    """Dihedral angle in degrees (IUPAC sign: positive for a right-handed twist) of rows of four points."""
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - np.sum(b0 * b1, axis=-1, keepdims=True) * b1
    w = b2 - np.sum(b2 * b1, axis=-1, keepdims=True) * b1
    x = np.sum(v * w, axis=-1)
    y = np.sum(np.cross(b1, v) * w, axis=-1)
    return np.asarray(np.degrees(np.arctan2(y, x)))


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
        >>> rng = np.random.default_rng(0)
        >>> ref_coords = rng.random((100, 3))
        >>> mobile_coords = rng.random((100, 3))
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

    # Optimal rotation R (minimizes ||ref - mobile @ R||); the reflection check is part of the helper
    rotation = _kabsch_rotation(ref_centered, mobile_centered)

    # Apply rotation: aligned = mobile @ R
    aligned_mobile = np.dot(mobile_centered, rotation)

    # Calculate RMSD
    diff = ref_centered - aligned_mobile
    rmsd = np.sqrt(np.sum(diff**2) / len(diff))

    return float(rmsd)


def compute_rmsd_simple(coords_ref: np.ndarray, coords_mobile: np.ndarray) -> float:
    """Calculate simple RMSD without structural alignment.

    Args:
        coords_ref: Reference coordinates
        coords_mobile: Mobile coordinates

    Returns:
        RMSD value

    Example:
        >>> rng = np.random.default_rng(0)
        >>> ref_coords = rng.random((50, 3))
        >>> mobile_coords = rng.random((50, 3))
        >>> rmsd = compute_rmsd_simple(ref_coords, mobile_coords)
        >>> rmsd >= 0
        True
    """
    if coords_ref.shape != coords_mobile.shape:
        raise ValueError("Coordinate arrays must have the same shape")

    diff = coords_ref - coords_mobile
    rmsd = np.sqrt(np.sum(diff**2) / len(diff))

    return float(rmsd)


def align_structures_kabsch(coords_ref: np.ndarray, coords_mobile: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """Align two protein structures using Kabsch algorithm.

    Args:
        coords_ref: Reference coordinates
        coords_mobile: Mobile coordinates to align

    Returns:
        Tuple of (aligned_coords, rotation_matrix, rmsd)

    Example:
        >>> rng = np.random.default_rng(0)
        >>> ref = rng.random((100, 3))
        >>> mobile = rng.random((100, 3))
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

    # Optimal rotation R (SVD of the covariance matrix, with the reflection check)
    rotation = _kabsch_rotation(ref_centered, mobile_centered)

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
        >>> rng = np.random.default_rng(0)
        >>> coords = rng.random((100, 3)) * 50
        >>> rg = calculate_radius_of_gyration(coords)
        >>> rg >= 0
        True
    """
    centroid = np.mean(coords, axis=0)
    distances = np.linalg.norm(coords - centroid, axis=1)
    rg = np.sqrt(np.mean(distances**2))

    return float(rg)


def calculate_center_of_mass(coords: np.ndarray, masses: Optional[np.ndarray] = None) -> np.ndarray:
    """Calculate center of mass for a protein structure.

    Args:
        coords: Atomic coordinates (n_atoms, 3)
        masses: Atomic masses (optional, defaults to equal masses)

    Returns:
        Center of mass coordinates (3,)

    Example:
        >>> rng = np.random.default_rng(0)
        >>> coords = rng.random((50, 3))
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
        >>> rng = np.random.default_rng(0)
        >>> coords = rng.random((50, 3))
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

    inertia = np.zeros((3, 3))
    inertia[0, 0] = np.sum(masses * (y**2 + z**2))
    inertia[1, 1] = np.sum(masses * (x**2 + z**2))
    inertia[2, 2] = np.sum(masses * (x**2 + y**2))
    inertia[0, 1] = inertia[1, 0] = -np.sum(masses * x * y)
    inertia[0, 2] = inertia[2, 0] = -np.sum(masses * x * z)
    inertia[1, 2] = inertia[2, 1] = -np.sum(masses * y * z)

    return inertia


def find_principal_axes(coords: np.ndarray, masses: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Find principal axes of a protein structure.

    Args:
        coords: Atomic coordinates
        masses: Atomic masses (optional)

    Returns:
        Tuple of (eigenvalues, eigenvectors)

    Example:
        >>> rng = np.random.default_rng(0)
        >>> coords = rng.random((50, 3))
        >>> eigenvals, eigenvecs = find_principal_axes(coords)
        >>> len(eigenvals) == 3
        True
    """
    inertia = calculate_inertia_tensor(coords, masses)

    # Diagonalize
    eigenvals, eigenvecs = np.linalg.eigh(inertia)

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
        >>> rng = np.random.default_rng(0)
        >>> coords = rng.random((100, 3)) * 50
        >>> stats = calculate_structural_statistics(coords)
        >>> "radius_of_gyration" in stats
        True
    """
    stats: Dict[str, Any] = {}

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


# ---------------------------------------------------------------------------------------------------------------------
# P-SEA secondary structure assignment from the C-alpha trace
# ---------------------------------------------------------------------------------------------------------------------

#: Consecutive C-alpha atoms closer than the first or farther than the second bound (angstrom) are not linked by a
#: virtual bond, so they mark a chain break. The ideal virtual bond is 3.8 A; a cis peptide shortens it to about 2.9 A.
_CA_BOND_RANGE: Tuple[float, float] = (2.9, 4.6)

#: P-SEA windows of an alpha-helix (Labesse et al. 1997), each as (centre, half-width). ``d2``, ``d3`` and ``d4`` are
#: the distances CA(i)-CA(i+2), CA(i)-CA(i+3) and CA(i)-CA(i+4) in angstrom, ``angle`` the virtual bond angle
#: CA(i-1)-CA(i)-CA(i+1) and ``torsion`` the virtual dihedral CA(i-1)-CA(i)-CA(i+1)-CA(i+2), both in degrees.
#: The numbers were transcribed from the P-SEA parameter set without access to the paper in this repository's
#: environment. The tests check that the geometry of an ideal helix and strand falls inside them; confirm them against
#: the table of the paper before relying on them for anything beyond that.
PSEA_HELIX_WINDOWS: Dict[str, Tuple[float, float]] = {
    "d2": (5.5, 0.5),
    "d3": (5.3, 0.5),
    "d4": (6.4, 0.6),
    "angle": (89.0, 12.0),
    "torsion": (50.0, 20.0),
}

#: P-SEA windows of a beta-strand (Labesse et al. 1997), in the same layout as :data:`PSEA_HELIX_WINDOWS`.
PSEA_STRAND_WINDOWS: Dict[str, Tuple[float, float]] = {
    "d2": (6.7, 0.6),
    "d3": (9.9, 0.9),
    "d4": (12.4, 1.1),
    "angle": (124.0, 14.0),
    "torsion": (-170.0, 45.0),
}

_PSEA_WINDOW = 5  # residues i..i+4 span the longest P-SEA distance, CA(i)-CA(i+4); also the shortest element


def _within_window(values: np.ndarray, window: Tuple[float, float], *, circular: bool = False) -> np.ndarray:
    """Whether ``values`` lie within ``centre +/- half_width`` (angles in degrees wrap around when ``circular``)."""
    centre, half_width = window
    delta = values - centre
    if circular:
        delta = (delta + 180.0) % 360.0 - 180.0
    return np.asarray(np.abs(delta) <= half_width)


def _covered_residues(window_ok: np.ndarray, width: int, n_residues: int) -> np.ndarray:
    """Residues in at least one satisfied window; ``window_ok[s]`` says the window of ``width`` residues at s holds."""
    starts = np.zeros(n_residues, dtype=int)
    starts[: len(window_ok)] = window_ok
    return np.asarray(np.convolve(starts, np.ones(width, dtype=int))[:n_residues] > 0)


def _true_runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    """Maximal runs of True in ``mask`` as (first index, last index)."""
    padded = np.concatenate([[False], np.asarray(mask, dtype=bool), [False]]).astype(int)
    edges = np.flatnonzero(np.diff(padded))
    return [(int(first), int(stop) - 1) for first, stop in zip(edges[::2], edges[1::2])]


def _chain_segments(ca: np.ndarray) -> List[Tuple[int, int]]:
    """(first, stop) index ranges of the runs of C-alpha atoms that are linked by virtual bonds of normal length."""
    if len(ca) == 0:
        return []
    bond = np.linalg.norm(ca[1:] - ca[:-1], axis=1)
    breaks = (np.flatnonzero((bond < _CA_BOND_RANGE[0]) | (bond > _CA_BOND_RANGE[1])) + 1).tolist()
    return list(zip([0, *breaks], [*breaks, len(ca)]))


def _psea_candidates(
    windows: Dict[str, Tuple[float, float]],
    distances: Tuple[np.ndarray, np.ndarray, np.ndarray],
    angle: np.ndarray,
    torsion: np.ndarray,
    n_residues: int,
) -> np.ndarray:
    """Residues of an unbroken segment (at least five residues) that lie in a window fitting one P-SEA conformation.

    A window of five residues i..i+4 fits a conformation when either its three distances d(i,i+2), d(i,i+3) and
    d(i,i+4) all lie in the windows of the conformation, or the virtual angles at CA(i+1) and CA(i+2) and the virtual
    torsions about CA(i+1)-CA(i+2) and CA(i+2)-CA(i+3) all do. The residues of every fitting window are candidates.
    """
    d2, d3, d4 = distances
    n_windows = n_residues - _PSEA_WINDOW + 1  # windows i = 0 .. n-5
    distance_ok = (
        _within_window(d2[:n_windows], windows["d2"])
        & _within_window(d3[:n_windows], windows["d3"])
        & _within_window(d4[:n_windows], windows["d4"])
    )
    # entry j: the angle at CA(j+1) and the torsion CA(j)..CA(j+3) fit; a window needs the entries j = i and j = i+1
    pair_ok = _within_window(angle, windows["angle"]) & _within_window(torsion, windows["torsion"], circular=True)
    angle_ok = pair_ok[:-1] & pair_ok[1:]
    return _covered_residues(distance_ok | angle_ok, _PSEA_WINDOW, n_residues)


def _psea_segment_states(ca: np.ndarray) -> np.ndarray:
    """P-SEA state (``"H"``, ``"E"`` or ``"C"``) of every residue of one unbroken C-alpha segment."""
    n_residues = len(ca)
    states = np.full(n_residues, "C", dtype="U1")
    if n_residues < _PSEA_WINDOW:
        return states

    d2 = np.linalg.norm(ca[2:] - ca[:-2], axis=1)  # entry i: CA(i)-CA(i+2)
    d3 = np.linalg.norm(ca[3:] - ca[:-3], axis=1)  # entry i: CA(i)-CA(i+3)
    d4 = np.linalg.norm(ca[4:] - ca[:-4], axis=1)  # entry i: CA(i)-CA(i+4)
    toward_previous = ca[:-2] - ca[1:-1]
    toward_next = ca[2:] - ca[1:-1]
    cosine = np.sum(toward_previous * toward_next, axis=1) / (
        np.linalg.norm(toward_previous, axis=1) * np.linalg.norm(toward_next, axis=1)
    )
    angle = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))[: n_residues - 3]  # entry j: angle at CA(j+1)
    torsion = _dihedral_degrees(ca[:-3], ca[1:-2], ca[2:-1], ca[3:])  # entry j: torsion of CA(j)..CA(j+3)

    helix = _psea_candidates(PSEA_HELIX_WINDOWS, (d2, d3, d4), angle, torsion, n_residues)
    # A residue of a helix is not a strand residue; what is left of a strand must still span a whole window.
    strand = _psea_candidates(PSEA_STRAND_WINDOWS, (d2, d3, d4), angle, torsion, n_residues) & ~helix
    for first, last in _true_runs(strand):
        if last - first + 1 >= _PSEA_WINDOW:
            states[first : last + 1] = "E"
    states[helix] = "H"
    return states


def assign_secondary_structure_psea(ca_coords: np.ndarray) -> List[str]:
    """Assign a secondary structure state to every residue of a C-alpha trace with the P-SEA algorithm.

    P-SEA (Labesse et al. 1997) needs nothing but C-alpha coordinates. For the residues of the trace it measures the
    distances CA(i)-CA(i+2), CA(i)-CA(i+3) and CA(i)-CA(i+4), the virtual bond angle at CA(i) and the virtual
    torsion about CA(i)-CA(i+1), and compares them with the windows of an alpha-helix (:data:`PSEA_HELIX_WINDOWS`:
    5.5 +/- 0.5, 5.3 +/- 0.5 and 6.4 +/- 0.6 A; 89 +/- 12 and 50 +/- 20 degrees) and of a beta-strand
    (:data:`PSEA_STRAND_WINDOWS`: 6.7 +/- 0.6, 9.9 +/- 0.9 and 12.4 +/- 1.1 A; 124 +/- 14 and -170 +/- 45 degrees).

    How the windows become an assignment:

    - A window of five residues i..i+4 fits a conformation when its three distances all lie in the windows of the
      conformation, or when the virtual bond angles at CA(i+1) and CA(i+2) and the virtual torsions about
      CA(i+1)-CA(i+2) and CA(i+2)-CA(i+3) all do. Every residue of a fitting window is a candidate for that
      conformation.
    - Candidate runs of the helix windows are helices (``"H"``). Of the remaining strand candidates, the runs of at
      least five residues are strands (``"E"``). Everything else is coil (``"C"``). No element is shorter than five
      residues, because the longest distance, CA(i)-CA(i+4), spans five.

    The windows are the P-SEA parameters (see the note at :data:`PSEA_HELIX_WINDOWS`). The five-residue window rule
    above, the minimum length, the priority of helix over strand and the break handling are this implementation's
    reading of the method, not text of the paper, so assignments can differ from the original program, mostly at the
    ends of elements.

    The trace is split into segments wherever consecutive C-alpha atoms are not 2.9-4.6 A apart (a chain break), and
    each segment is assigned on its own, so no element reaches across a gap. The distance criterion does not depend
    on handedness: the mirror image of a helix passes it. The labels are those of
    :mod:`metainformant.protein.structure.secondary`; P-SEA's own letters are a, b and c for ``"H"``, ``"E"`` and
    ``"C"``.

    Args:
        ca_coords: C-alpha coordinates (n_residues, 3) in chain order, in angstrom

    Returns:
        One state per residue, ``"H"``, ``"E"`` or ``"C"``

    Raises:
        ValueError: If ``ca_coords`` is not an (n, 3) array of finite numbers

    Reference:
        Labesse G, Colloc'h N, Pothier J, Mornon JP. P-SEA: a new efficient assignment of secondary structure from
        C alpha trace of proteins. Comput Appl Biosci 1997;13(3):291-295.

    Example:
        >>> turn = np.radians(100.0) * np.arange(8)  # an ideal alpha-helix: 2.3 A radius, 1.5 A rise, 100 degrees
        >>> helix = np.column_stack([2.3 * np.cos(turn), 2.3 * np.sin(turn), 1.5 * np.arange(8)])
        >>> assign_secondary_structure_psea(helix)
        ['H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']
    """
    ca = _as_points(ca_coords, "ca_coords")
    states = np.full(len(ca), "C", dtype="U1")
    for first, stop in _chain_segments(ca):
        states[first:stop] = _psea_segment_states(ca[first:stop])
    return [str(state) for state in states]


def identify_secondary_structure_elements(
    coords: np.ndarray, ca_indices: Optional[Sequence[int]] = None
) -> List[Dict[str, Any]]:
    """Identify helices and strands in a C-alpha trace with the P-SEA algorithm.

    This groups the per-residue states of :func:`assign_secondary_structure_psea` into elements; see there for the
    criteria. P-SEA provides no confidence, so none is reported.

    Args:
        coords: Coordinates (n, 3). Without ``ca_indices`` they are the C-alpha atoms themselves, one per residue in
            chain order.
        ca_indices: Optional indices of the C-alpha atoms in ``coords``, one per residue in chain order (so ``coords``
            may hold all atoms). Other backbone atoms must not be listed: the criteria are defined on C-alpha atoms.

    Returns:
        Elements in chain order, each a dict with ``type`` (``"helix"`` or ``"strand"``), ``start_residue`` and
        ``end_residue`` (0-based positions in the C-alpha list, both included) and ``length`` (residues)

    Raises:
        ValueError: If the coordinates are not an (n, 3) array of finite numbers, or an index is out of range

    Example:
        >>> turn = np.radians(100.0) * np.arange(12)  # an ideal alpha-helix: 2.3 A radius, 1.5 A rise, 100 degrees
        >>> helix = np.column_stack([2.3 * np.cos(turn), 2.3 * np.sin(turn), 1.5 * np.arange(12)])
        >>> identify_secondary_structure_elements(helix)
        [{'type': 'helix', 'start_residue': 0, 'end_residue': 11, 'length': 12}]
    """
    points = _as_points(coords)
    if ca_indices is None:
        ca = points
    else:
        indices = np.asarray(ca_indices)
        if indices.ndim != 1 or (indices.size and not np.issubdtype(indices.dtype, np.integer)):
            raise ValueError("ca_indices must be a one-dimensional sequence of integers")
        if indices.size and (indices.min() < 0 or indices.max() >= len(points)):
            raise ValueError(f"ca_indices must lie in [0, {len(points)}), got {indices.min()}..{indices.max()}")
        ca = points[indices.astype(int)]

    states = np.array(assign_secondary_structure_psea(ca), dtype="U1")
    elements: List[Dict[str, Any]] = []
    for offset, stop in _chain_segments(ca):  # an element never reaches across a chain break
        for state, kind in (("H", "helix"), ("E", "strand")):
            for first, last in _true_runs(states[offset:stop] == state):
                elements.append(
                    {
                        "type": kind,
                        "start_residue": offset + first,
                        "end_residue": offset + last,
                        "length": last - first + 1,
                    }
                )
    elements.sort(key=lambda element: element["start_residue"])
    return elements


# ---------------------------------------------------------------------------------------------------------------------
# Solvent accessible surface area (Shrake-Rupley)
# ---------------------------------------------------------------------------------------------------------------------

#: Van der Waals radii in angstrom by element symbol (upper case): Bondi, J. Phys. Chem. 68, 441 (1964), with the
#: calcium radius of Mantina et al., J. Phys. Chem. A 113, 5806 (2009).
VAN_DER_WAALS_RADII: Dict[str, float] = {
    "H": 1.20,
    "HE": 1.40,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "NA": 2.27,
    "MG": 1.73,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "K": 2.75,
    "CA": 2.31,
    "NI": 1.63,
    "CU": 1.40,
    "ZN": 1.39,
    "SE": 1.90,
    "BR": 1.85,
    "CD": 1.58,
    "I": 1.98,
    "HG": 1.55,
}

#: Default water probe radius in angstrom
WATER_PROBE_RADIUS = 1.4

#: Default number of test points per atomic sphere in the Shrake-Rupley calculation
DEFAULT_SASA_POINTS = 960


def van_der_waals_radius(element: str) -> float:
    """Van der Waals radius in angstrom of an element given by its symbol (case-insensitive).

    Raises:
        ValueError: If :data:`VAN_DER_WAALS_RADII` has no radius for the element; nothing is guessed
    """
    symbol = str(element).strip().upper()
    if symbol not in VAN_DER_WAALS_RADII:
        raise ValueError(
            f"no van der Waals radius for element {element!r}; known elements: {', '.join(sorted(VAN_DER_WAALS_RADII))}"
        )
    return VAN_DER_WAALS_RADII[symbol]


def golden_spiral_points(n_points: int) -> np.ndarray:
    """``n_points`` nearly evenly spread unit vectors on the sphere (golden spiral / Fibonacci lattice).

    Point k has ``z = 1 - (2k + 1) / n_points`` and longitude ``k * pi * (3 - sqrt(5))``, so every point stands for the
    same area of the sphere. The set is deterministic.

    Raises:
        ValueError: If ``n_points`` is not a positive integer

    Example:
        >>> golden_spiral_points(4).shape
        (4, 3)
    """
    if isinstance(n_points, bool) or int(n_points) != n_points or n_points < 1:
        raise ValueError(f"n_points must be a positive integer, got {n_points!r}")
    k = np.arange(int(n_points), dtype=float)
    z = 1.0 - (2.0 * k + 1.0) / n_points
    ring_radius = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    longitude = k * math.pi * (3.0 - math.sqrt(5.0))
    return np.column_stack([ring_radius * np.cos(longitude), ring_radius * np.sin(longitude), z])


def _atomic_radii(radii: Union[None, float, Sequence[float], np.ndarray], n_atoms: int) -> np.ndarray:
    """Per-atom radii: the carbon radius for all when ``radii`` is None, one value for all, or one per atom."""
    if radii is None:
        values = np.full(n_atoms, VAN_DER_WAALS_RADII["C"])
    else:
        values = np.asarray(radii, dtype=float)
        if values.ndim == 0:
            values = np.full(n_atoms, float(values))
        elif values.shape != (n_atoms,):
            raise ValueError(f"radii must be a number or have one value per atom ({n_atoms}), got shape {values.shape}")
    if not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError("radii must be positive and finite")
    return values


def calculate_atomic_sasa(
    coords: np.ndarray,
    radii: Union[None, float, Sequence[float], np.ndarray] = None,
    probe_radius: float = WATER_PROBE_RADIUS,
    n_points: int = DEFAULT_SASA_POINTS,
) -> np.ndarray:
    """Solvent accessible surface area of every atom by the Shrake-Rupley algorithm.

    Each atom is a sphere of its van der Waals radius plus the probe radius. ``n_points`` test points are spread over
    the surface of each such sphere (:func:`golden_spiral_points`); a point is accessible when it lies inside no other
    atom's sphere, and the atom's area is the sphere area ``4 pi (r + probe)^2`` times the accessible fraction of the
    points. The result approximates the exact area, and the approximation improves as ``n_points`` grows; an isolated
    atom is exact and two overlapping spheres match the analytic spherical-cap formula to about ``1 / n_points``.

    Args:
        coords: Atom centres (n_atoms, 3) in angstrom
        radii: Van der Waals radii in angstrom: one per atom, one value for all, or None. Bare coordinates carry no
            element, so None gives every atom the carbon radius (1.70 A); use :func:`van_der_waals_radius` to build
            radii from elements, or :func:`metainformant.protein.structure.analysis.calculate_structure_sasa` for a
            structure dictionary.
        probe_radius: Probe radius in angstrom (1.4 A is a water molecule); 0 gives the van der Waals surface
        n_points: Test points per sphere

    Returns:
        Accessible area of each atom in square angstrom, shape (n_atoms,)

    Raises:
        ValueError: If the coordinates are not (n, 3) and finite, a radius is not positive, the probe radius is
            negative or ``n_points`` is not a positive integer

    Reference:
        Shrake A, Rupley JA. Environment and exposure to solvent of protein atoms. Lysozyme and insulin.
        J Mol Biol 1973;79:351-371.

    Example:
        >>> centres = np.array([[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]])  # two atoms too far apart to touch
        >>> bool(np.allclose(calculate_atomic_sasa(centres, radii=1.6), 4 * np.pi * (1.6 + 1.4) ** 2))
        True
    """
    points = _as_points(coords)
    n_atoms = len(points)
    if not np.isfinite(probe_radius) or probe_radius < 0:
        raise ValueError(f"probe_radius must be a non-negative number, got {probe_radius!r}")
    sphere = golden_spiral_points(n_points)
    atom_radii = _atomic_radii(radii, n_atoms)
    if n_atoms == 0:
        return np.zeros(0)

    expanded = atom_radii + probe_radius  # the sphere traced by the centre of the probe
    areas = np.empty(n_atoms)
    neighbours = cKDTree(points).query_ball_point(points, expanded + expanded.max())
    for i in range(n_atoms):
        others = np.array([j for j in neighbours[i] if j != i], dtype=int)
        if others.size:
            others = others[np.linalg.norm(points[others] - points[i], axis=1) < expanded[i] + expanded[others]]
        n_accessible = len(sphere)
        if others.size:
            surface = points[i] + expanded[i] * sphere
            squared = np.sum((surface[:, np.newaxis, :] - points[others][np.newaxis, :, :]) ** 2, axis=2)
            n_accessible -= int(np.count_nonzero(np.any(squared < expanded[others] ** 2, axis=1)))
        areas[i] = 4.0 * math.pi * expanded[i] ** 2 * n_accessible / len(sphere)
    return areas


def calculate_solvent_accessible_surface_area(
    coords: np.ndarray,
    probe_radius: float = WATER_PROBE_RADIUS,
    radii: Union[None, float, Sequence[float], np.ndarray] = None,
    n_points: int = DEFAULT_SASA_POINTS,
) -> float:
    """Total solvent accessible surface area of a set of atoms (Shrake-Rupley, see :func:`calculate_atomic_sasa`).

    The area follows from the atomic radii, the probe and every concavity of the atom set; the area of a convex hull
    of the centres would know none of them.

    Args:
        coords: Atom centres (n_atoms, 3) in angstrom
        probe_radius: Probe radius in angstrom (default: 1.4 A for water)
        radii: Van der Waals radii in angstrom (one per atom, one for all) or None for the carbon radius of every atom
        n_points: Test points per atomic sphere; more points give a closer approximation

    Returns:
        Surface area in square angstrom (0.0 for no atoms)

    Example:
        >>> area = calculate_solvent_accessible_surface_area(np.zeros((1, 3)), radii=1.7)  # one carbon atom
        >>> abs(area - 4 * math.pi * (1.7 + 1.4) ** 2) < 1e-9
        True
    """
    return math.fsum(calculate_atomic_sasa(coords, radii=radii, probe_radius=probe_radius, n_points=n_points))


# ---------------------------------------------------------------------------------------------------------------------
# TM-score
# ---------------------------------------------------------------------------------------------------------------------

#: Smallest d0 of the TM-score in angstrom; the published d0 formula is not used at or below 21 residues
_TM_D0_MIN = 0.5
_TM_MIN_ALIGNED = 3
_TM_MAX_ITERATIONS = 20


def tm_score_d0(n_reference: int) -> float:
    """Distance scale d0 of the TM-score for a reference structure of ``n_reference`` residues, in angstrom.

    ``d0 = 1.24 * (L - 15) ** (1 / 3) - 1.8`` (Zhang & Skolnick 2004). The formula falls below 0.5 A at 21 residues
    and is meaningless below 16, so, as in the published TM-score program, d0 is 0.5 A for ``L <= 21``.

    Raises:
        ValueError: If ``n_reference`` is below 1

    Example:
        >>> round(tm_score_d0(100), 4)
        3.6521
        >>> tm_score_d0(20)
        0.5
    """
    if n_reference < 1:
        raise ValueError(f"n_reference must be at least 1, got {n_reference}")
    if n_reference <= 21:
        return _TM_D0_MIN
    return float(1.24 * (n_reference - 15) ** (1.0 / 3.0) - 1.8)


def tm_score_from_distances(distances: Sequence[float], n_reference: int) -> float:
    """TM-score of aligned residue pairs at the given distances, without any superposition.

    ``TM = (1 / L) * sum_i 1 / (1 + (d_i / d0) ** 2)`` over the aligned pairs, with ``L = n_reference`` the number of
    residues of the reference structure (residues without a partner count as 0) and d0 from :func:`tm_score_d0`.

    Args:
        distances: Distance in angstrom of each aligned pair after superposition
        n_reference: Number of residues of the reference structure, at least the number of pairs

    Raises:
        ValueError: If a distance is negative or not finite, or there are more pairs than reference residues

    Example:
        >>> d0 = tm_score_d0(22)
        >>> tm_score_from_distances([0.0] * 11 + [d0] * 11, 22)  # 11 perfect pairs and 11 at d0: (11 + 11 / 2) / 22
        0.75
    """
    values = np.asarray(distances, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("distances must be a one-dimensional sequence of finite, non-negative numbers")
    if len(values) > n_reference:
        raise ValueError(f"{len(values)} aligned pairs cannot exceed the {n_reference} residues of the reference")
    d0 = tm_score_d0(n_reference)
    return float(np.sum(1.0 / (1.0 + (values / d0) ** 2)) / n_reference)


def _within_cutoff(distances: np.ndarray, cutoff: float) -> np.ndarray:
    """Pairs closer than ``cutoff``; the cutoff grows in steps of 0.5 A until at least three pairs qualify."""
    selected = distances < cutoff
    while np.count_nonzero(selected) < min(_TM_MIN_ALIGNED, len(distances)):
        cutoff += 0.5
        selected = distances < cutoff
    return np.asarray(selected)


def _tm_score_search(
    coords_ref: np.ndarray, coords_mobile: np.ndarray, n_reference: int
) -> Tuple[float, np.ndarray, np.ndarray]:
    """Superposition that maximises the TM-score, found as in the TM-score program of Zhang & Skolnick.

    The superposition is seeded with every fragment of length n, n/2, n/4, n/8, n/16 and 4 residues (halving stops
    at 4). Each seed is refined for up to 20 rounds: the pairs closer than a cutoff of max(min(d0, 8), 4.5) A
    (widened until three pairs qualify) are superposed with the Kabsch algorithm, and the TM-score of the result is
    kept if it is the best so far. The score is the exact TM-score of the returned superposition, and a lower bound
    of the maximum over all superpositions.
    """
    n_aligned = len(coords_ref)
    d0 = tm_score_d0(n_reference)
    search_cutoff = min(max(d0, 4.5), 8.0)
    smallest = min(4, n_aligned)

    lengths: List[int] = []
    for halvings in range(5):
        length = n_aligned // 2**halvings
        if length <= smallest:
            lengths.append(smallest)
            break
        lengths.append(length)
    else:
        lengths.append(smallest)

    def score_of(distances: np.ndarray) -> float:
        return float(np.sum(1.0 / (1.0 + (distances / d0) ** 2)) / n_reference)

    best_score = -1.0
    best_rotation = np.eye(3)
    best_translation = np.zeros(3)
    perfect = n_aligned / n_reference - 1e-12  # no superposition can score higher

    for length in lengths:
        for first in range(n_aligned - length + 1):
            fragment = slice(first, first + length)
            rotation, translation = _kabsch_transform(coords_ref[fragment], coords_mobile[fragment])
            distances = np.linalg.norm(np.dot(coords_mobile, rotation) + translation - coords_ref, axis=1)
            previous: Optional[np.ndarray] = None
            for step in range(_TM_MAX_ITERATIONS + 1):
                score = score_of(distances)
                if score > best_score:
                    best_score, best_rotation, best_translation = score, rotation, translation
                    if best_score >= perfect:
                        return best_score, best_rotation, best_translation
                if step == _TM_MAX_ITERATIONS:
                    break
                selected = _within_cutoff(distances, search_cutoff)
                if previous is not None and np.array_equal(selected, previous):
                    break  # the superposition would not change any more
                previous = selected
                rotation, translation = _kabsch_transform(coords_ref[selected], coords_mobile[selected])
                distances = np.linalg.norm(np.dot(coords_mobile, rotation) + translation - coords_ref, axis=1)
    return best_score, best_rotation, best_translation


def calculate_tm_score(
    coords_ref: np.ndarray, coords_mobile: np.ndarray, n_reference: Optional[int] = None
) -> Dict[str, Any]:
    """TM-score of two structures whose residues correspond one to one (Zhang & Skolnick 2004).

    The TM-score is the maximum over rigid-body superpositions of ``(1 / L) * sum_i 1 / (1 + (d_i / d0) ** 2)``, where
    ``d_i`` are the distances of the paired residues after superposition, ``L`` is the length of the reference and
    ``d0`` comes from :func:`tm_score_d0`. It is 1 for identical structures, independent of protein length, and a
    score below about 0.17 is what unrelated structures give. The maximising superposition is searched for as the
    TM-score program does (fragment seeds refined by iterated Kabsch superposition), so the value is the TM-score of an
    actual superposition and a lower bound of the true maximum, and it does not change when either structure is
    rotated or translated.

    Args:
        coords_ref: Reference positions (n, 3), one per residue (normally the C-alpha atoms), in angstrom
        coords_mobile: Positions of the corresponding residues of the other structure, same shape
        n_reference: Length ``L`` that normalises the score: the number of residues of the whole reference structure.
            Default: ``n``, which is right when every reference residue has a partner.

    Returns:
        Dictionary with ``tm_score``, ``d0`` (angstrom), ``n_aligned``, ``n_reference``, and the superposition that
        gives the score as ``rotation`` (3, 3) and ``translation`` (3,): ``coords_mobile @ rotation + translation`` lies
        on ``coords_ref``

    Raises:
        ValueError: If the arrays differ in shape, are not (n, 3) and finite, hold fewer than 3 pairs, or
            ``n_reference`` is smaller than the number of pairs

    Reference:
        Zhang Y, Skolnick J. Scoring function for automated assessment of protein structure template quality.
        Proteins 2004;57:702-710.

    Example:
        >>> line = np.column_stack([3.8 * np.arange(30), np.zeros(30), np.zeros(30)])
        >>> result = calculate_tm_score(line, line + np.array([5.0, -2.0, 1.0]))  # the same chain, moved
        >>> round(result["tm_score"], 6)
        1.0
    """
    ref = _as_points(coords_ref, "coords_ref")
    mobile = _as_points(coords_mobile, "coords_mobile")
    if ref.shape != mobile.shape:
        raise ValueError(f"coordinate arrays must have the same shape, got {ref.shape} and {mobile.shape}")
    n_aligned = len(ref)
    if n_aligned < _TM_MIN_ALIGNED:
        raise ValueError(f"the TM-score needs at least {_TM_MIN_ALIGNED} aligned residues, got {n_aligned}")
    length = n_aligned if n_reference is None else int(n_reference)
    if length < n_aligned:
        raise ValueError(f"n_reference ({length}) cannot be smaller than the number of aligned residues ({n_aligned})")

    score, rotation, translation = _tm_score_search(ref, mobile, length)
    return {
        "tm_score": score,
        "d0": tm_score_d0(length),
        "n_aligned": n_aligned,
        "n_reference": length,
        "rotation": rotation,
        "translation": translation,
    }
