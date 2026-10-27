"""Protein structure analysis utilities.

This module provides tools for analyzing protein 3D structures,
including domain identification, surface analysis, and structural motifs.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def calculate_contact_map(coords: np.ndarray, threshold: float = 8.0) -> np.ndarray:
    """Calculate residue contact map from atomic coordinates.

    Args:
        coords: Atomic coordinates (n_atoms, 3)
        threshold: Distance threshold for contacts (Å)

    Returns:
        Contact map matrix (n_residues, n_residues)

    Example:
        >>> rng = np.random.default_rng(0)
        >>> coords = rng.random((50, 3)) * 50
        >>> contact_map = calculate_contact_map(coords, threshold=10.0)
        >>> contact_map.shape[0] == contact_map.shape[1]
        True
    """
    if coords.shape[0] < 2:
        return np.array([[0]])

    # Calculate pairwise distances between all atoms
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff**2, axis=2))

    # Convert to contact map (1 if within threshold, 0 otherwise)
    contact_map = (distances <= threshold).astype(int)

    # Set diagonal to 0 (no self-contacts)
    np.fill_diagonal(contact_map, 0)

    return contact_map


def identify_domains(structure: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Identify structural domains in a protein structure.

    Args:
        structure: Structure dictionary from PDB loading

    Returns:
        List of identified domains

    Example:
        >>> # Assuming structure data exists
        >>> # domains = identify_domains(structure)
        >>> # isinstance(domains, list)
        >>> # True
    """
    # This is a simplified domain identification
    # Real implementation would use domain databases or clustering

    atoms = structure.get("atoms", [])
    if not atoms:
        return []

    # Group atoms by residue
    residues: Dict[Tuple[Any, Any], List[Any]] = {}
    for atom in atoms:
        key = (atom["chain_id"], atom["res_seq"])
        if key not in residues:
            residues[key] = []
        residues[key].append(atom)

    # Simple domain identification based on chain breaks
    domains = []
    current_domain = []
    prev_chain = None
    prev_res_seq: Any = None

    for (chain_id, res_seq), residue_atoms in sorted(residues.items()):
        if prev_chain is not None and (chain_id != prev_chain or res_seq > prev_res_seq + 10):
            # Chain break or large gap - start new domain
            if current_domain:
                domains.append(
                    {
                        "residues": current_domain,
                        "chain_id": prev_chain,
                        "start": min(current_domain),
                        "end": max(current_domain),
                        "size": len(current_domain),
                    }
                )
            current_domain = []

        current_domain.append(res_seq)
        prev_chain = chain_id
        prev_res_seq = res_seq

    # Add final domain
    if current_domain:
        domains.append(
            {
                "residues": current_domain,
                "chain_id": prev_chain,
                "start": min(current_domain),
                "end": max(current_domain),
                "size": len(current_domain),
            }
        )

    return domains


def calculate_surface_area(coords: np.ndarray, probe_radius: float = 1.4) -> float:
    """Calculate solvent accessible surface area.

    Args:
        coords: Atomic coordinates
        probe_radius: Solvent probe radius (Å)

    Returns:
        Surface area in Å²

    Example:
        >>> rng = np.random.default_rng(0)
        >>> coords = rng.random((100, 3)) * 30
        >>> area = calculate_surface_area(coords)
        >>> area > 0
        True
    """
    # Simplified SASA calculation
    # Real implementation would use more sophisticated algorithms

    if len(coords) < 4:
        # For small structures, approximate as sphere
        if len(coords) > 0:
            radius = np.max(np.linalg.norm(coords - np.mean(coords, axis=0), axis=1))
            return float(4 * np.pi * (radius + probe_radius) ** 2)
        return 0.0

    # Use convex hull as approximation
    try:
        from scipy.spatial import ConvexHull

        hull = ConvexHull(coords)

        # Approximate surface area (this is rough)
        # In practice, would need proper SASA calculation
        surface_area = hull.area

        # Add probe radius effect (simplified)
        surface_area *= 1 + probe_radius / 10.0

        return float(surface_area)

    except ImportError:
        # Fallback to bounding box approximation
        min_coords = np.min(coords, axis=0)
        max_coords = np.max(coords, axis=0)
        dimensions = max_coords - min_coords

        # Approximate as rectangular box surface area
        surface_area = 2 * (
            dimensions[0] * dimensions[1] + dimensions[1] * dimensions[2] + dimensions[2] * dimensions[0]
        )

        return float(surface_area)


#: Motif types :func:`analyze_structural_motifs` searches for by default.
SUPPORTED_MOTIF_TYPES: Tuple[str, ...] = ("helix", "extended", "zinc_site")

# C-alpha trace geometry. An ideal right-handed alpha-helix (radius 2.3 A, rise 1.5 A, 100 degrees per residue) has
# CA(i)-CA(i+3) = 5.05 A, CA(i)-CA(i+4) = 6.20 A and a virtual torsion of +50 degrees over four consecutive C-alpha
# atoms; an extended strand has CA(i)-CA(i+2) of about 6.6 A and CA(i)-CA(i+3) of about 10 A. The windows below
# are those ideal values with a margin that admits real, distorted elements (distances in A, torsions in degrees).
_CA_BOND_RANGE = (2.9, 4.6)  # consecutive C-alpha atoms outside it are not bonded: a chain break (ideal 3.8, cis 2.9)
_HELIX_D3_WINDOW = (4.5, 5.7)
_HELIX_D4_WINDOW = (5.5, 6.9)
_HELIX_TORSION_WINDOW = (20.0, 80.0)
_EXTENDED_D2_WINDOW = (6.0, 7.5)
_EXTENDED_D3_WINDOW = (8.8, 11.0)
_MIN_HELIX_RESIDUES = 6
_MIN_EXTENDED_RESIDUES = 4

# Zinc ligands: (residue name, atom name) -> one-letter code of the ligating residue type.
_ZINC_LIGAND_ATOMS = {
    ("CYS", "SG"): "C",
    ("HIS", "ND1"): "H",
    ("HIS", "NE2"): "H",
    ("ASP", "OD1"): "D",
    ("ASP", "OD2"): "D",
    ("GLU", "OE1"): "E",
    ("GLU", "OE2"): "E",
}
_ZINC_LIGAND_CUTOFF = 2.8  # A; Zn-S(Cys) is about 2.3 A, Zn-N(His) 2.0-2.1 A and Zn-O(Asp/Glu) 2.0-2.2 A


def _record_type(atom: Dict[str, Any]) -> str:
    """Record type (``ATOM`` or ``HETATM``) of an atom dict from either PDB loader of this package."""
    return str(atom.get("record_name") or atom.get("record_type") or "ATOM").strip().upper()


def _insertion_code(atom: Dict[str, Any]) -> str:
    return str(atom.get("icode") or atom.get("ins_code") or "").strip()


def _in_window(values: np.ndarray, window: Tuple[float, float]) -> np.ndarray:
    return (values >= window[0]) & (values <= window[1])


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


def _ca_segments(atoms: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """C-alpha atoms of the protein residues, per chain and split wherever consecutive C-alpha atoms are not bonded.

    Only ATOM records count, so a calcium ion (atom name ``CA`` in a HETATM record) is not a C-alpha, and only
    the first alternate location of a residue is used.
    """
    seen = set()
    chains: Dict[str, List[Tuple[int, np.ndarray]]] = {}
    for atom in atoms:
        if _record_type(atom) != "ATOM" or str(atom.get("name", "")).strip().upper() != "CA":
            continue
        residue = (str(atom.get("chain_id", "")), atom["res_seq"], _insertion_code(atom))
        if residue in seen:
            continue
        seen.add(residue)
        position = np.array([atom["x"], atom["y"], atom["z"]], dtype=float)
        chains.setdefault(residue[0], []).append((atom["res_seq"], position))

    segments: List[Dict[str, Any]] = []
    for chain_id, residues in chains.items():
        coords = np.array([position for _, position in residues])
        bond = np.linalg.norm(coords[1:] - coords[:-1], axis=1)
        starts = [0, *(np.flatnonzero((bond < _CA_BOND_RANGE[0]) | (bond > _CA_BOND_RANGE[1])) + 1).tolist()]
        for first, stop in zip(starts, [*starts[1:], len(residues)]):
            segments.append(
                {
                    "chain_id": chain_id,
                    "numbers": [number for number, _ in residues[first:stop]],
                    "coords": coords[first:stop],
                }
            )
    return segments


def _helix_residue_mask(coords: np.ndarray) -> np.ndarray:
    """Residues in a window of five whose C-alpha atoms form an alpha-helical turn (see the constants above)."""
    n = len(coords)
    mask = np.zeros(n, dtype=bool)
    if n < 5:
        return mask
    d3 = np.linalg.norm(coords[3:] - coords[:-3], axis=1)  # entry i: CA(i)-CA(i+3)
    d4 = np.linalg.norm(coords[4:] - coords[:-4], axis=1)  # entry i: CA(i)-CA(i+4)
    torsion = _dihedral_degrees(coords[:-3], coords[1:-2], coords[2:-1], coords[3:])  # entry i: CA(i..i+3)
    windows = n - 4
    ok = (
        _in_window(d3[:windows], _HELIX_D3_WINDOW)
        & _in_window(d3[1 : windows + 1], _HELIX_D3_WINDOW)
        & _in_window(d4, _HELIX_D4_WINDOW)
        & _in_window(torsion[:windows], _HELIX_TORSION_WINDOW)
        & _in_window(torsion[1 : windows + 1], _HELIX_TORSION_WINDOW)
    )
    for first in np.flatnonzero(ok):
        mask[first : first + 5] = True
    return mask


def _extended_residue_mask(coords: np.ndarray) -> np.ndarray:
    """Residues in a window of four whose C-alpha atoms are stretched out like a strand."""
    n = len(coords)
    mask = np.zeros(n, dtype=bool)
    if n < 4:
        return mask
    d2 = np.linalg.norm(coords[2:] - coords[:-2], axis=1)  # entry i: CA(i)-CA(i+2)
    d3 = np.linalg.norm(coords[3:] - coords[:-3], axis=1)  # entry i: CA(i)-CA(i+3)
    windows = n - 3
    ok = (
        _in_window(d2[:windows], _EXTENDED_D2_WINDOW)
        & _in_window(d2[1 : windows + 1], _EXTENDED_D2_WINDOW)
        & _in_window(d3, _EXTENDED_D3_WINDOW)
    )
    for first in np.flatnonzero(ok):
        mask[first : first + 4] = True
    return mask


def _runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    """Maximal runs of True in ``mask`` as (first index, last index)."""
    padded = np.concatenate([[False], mask, [False]]).astype(int)
    edges = np.flatnonzero(np.diff(padded))
    return [(int(first), int(stop) - 1) for first, stop in zip(edges[::2], edges[1::2])]


def _secondary_structure_elements(atoms: Sequence[Dict[str, Any]], wanted: Sequence[str]) -> List[Dict[str, Any]]:
    """Helices and extended segments of the C-alpha trace, in chain order."""
    detectors = (
        ("helix", _helix_residue_mask, _MIN_HELIX_RESIDUES),
        ("extended", _extended_residue_mask, _MIN_EXTENDED_RESIDUES),
    )
    found: List[Tuple[int, int, Dict[str, Any]]] = []
    for segment_index, segment in enumerate(_ca_segments(atoms)):
        for motif_type, residue_mask, minimum in detectors:
            if motif_type not in wanted:
                continue
            for first, last in _runs(residue_mask(segment["coords"])):
                if last - first + 1 < minimum:
                    continue
                element = {
                    "type": motif_type,
                    "chain_id": segment["chain_id"],
                    "start_residue": segment["numbers"][first],
                    "end_residue": segment["numbers"][last],
                    "length": last - first + 1,
                }
                found.append((segment_index, first, element))
    found.sort(key=lambda item: (item[0], item[1], item[2]["type"]))
    return [element for _, _, element in found]


def _is_zinc(atom: Dict[str, Any]) -> bool:
    element = str(atom.get("element") or "").strip().upper()
    if element:
        return element == "ZN"
    return _record_type(atom) == "HETATM" and str(atom.get("name", "")).strip().upper() == "ZN"


def _zinc_sites(atoms: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Zinc ions with the Cys, His, Asp and Glu side-chain atoms that coordinate them."""
    ligand_atoms = [
        (
            atom,
            _ZINC_LIGAND_ATOMS[
                (str(atom.get("res_name", "")).strip().upper(), str(atom.get("name", "")).strip().upper())
            ],
        )
        for atom in atoms
        if _record_type(atom) == "ATOM"
        and (str(atom.get("res_name", "")).strip().upper(), str(atom.get("name", "")).strip().upper())
        in _ZINC_LIGAND_ATOMS
    ]
    sites: List[Dict[str, Any]] = []
    for ion in atoms:
        if not _is_zinc(ion):
            continue
        ion_position = np.array([ion["x"], ion["y"], ion["z"]], dtype=float)
        ligands: List[Dict[str, Any]] = []
        residue_letters: Dict[Tuple[str, Any, str], str] = {}
        for atom, letter in ligand_atoms:
            distance = float(np.linalg.norm(np.array([atom["x"], atom["y"], atom["z"]], dtype=float) - ion_position))
            if distance <= _ZINC_LIGAND_CUTOFF:
                ligands.append(
                    {
                        "chain_id": str(atom.get("chain_id", "")),
                        "res_seq": atom["res_seq"],
                        "res_name": str(atom.get("res_name", "")).strip().upper(),
                        "atom": str(atom.get("name", "")).strip().upper(),
                        "distance": distance,
                    }
                )
                residue_letters[(str(atom.get("chain_id", "")), atom["res_seq"], _insertion_code(atom))] = letter
        if not ligands:
            continue
        counts = {letter: sum(1 for value in residue_letters.values() if value == letter) for letter in "CHDE"}
        pattern = "".join(f"{letter}{count if count > 1 else ''}" for letter, count in counts.items() if count)
        ligands.sort(key=lambda ligand: (ligand["distance"], ligand["chain_id"], ligand["res_seq"], ligand["atom"]))
        sites.append(
            {
                "type": "zinc_site",
                "chain_id": str(ion.get("chain_id", "")),
                "ion_res_seq": ion.get("res_seq"),
                "ion_serial": ion.get("serial"),
                "position": ion_position.tolist(),
                "n_ligand_residues": len(residue_letters),
                "coordination_pattern": pattern,
                # Four Cys/His ligands (C2H2, C3H, C4, ...) are the signature of a structural zinc finger site
                "zinc_finger_like": len(residue_letters) == 4 and set(residue_letters.values()) <= {"C", "H"},
                "ligands": ligands,
            }
        )
    return sites


def analyze_structural_motifs(
    structure: Dict[str, Any], motif_types: Optional[Sequence[str]] = None
) -> List[Dict[str, Any]]:
    """Find structural motifs from the atom coordinates of a structure.

    Three motif types are searched for, each by a geometric criterion on the atoms given:

    - ``"helix"``: right-handed alpha-helices, found on the C-alpha trace. A window of five residues belongs
      to a helix when CA(i)-CA(i+3) and CA(i+1)-CA(i+4) are 4.5-5.7 A, CA(i)-CA(i+4) is 5.5-6.9 A and the virtual
      torsions of the C-alpha atoms i..i+3 and i+1..i+4 are 20-80 degrees (an ideal helix has 5.05 A, 6.20 A and
      +50 degrees). Overlapping windows are merged and a helix has at least 6 residues. A left-handed helix has
      negative torsions and is not reported.
    - ``"extended"``: extended, strand-like segments: windows of four residues with CA(i)-CA(i+2) and
      CA(i+1)-CA(i+3) of 6.0-7.5 A and CA(i)-CA(i+3) of 8.8-11.0 A (an ideal strand has 6.6 A and 10 A), at
      least 4 residues. This is the stretched backbone geometry only: pairing of strands into a sheet is
      not tested, so an extended loop is reported as well.
    - ``"zinc_site"``: zinc ions (``ZN``) with the Cys SG, His ND1/NE2 and Asp/Glu carboxylate oxygen atoms
      within 2.8 A, as ligands. The result gives the coordination pattern (``"C2H2"``: two Cys and two His
      residues, a residue counted once however many of its atoms bind), and ``zinc_finger_like`` is true for
      four Cys/His ligands, the signature of a structural zinc finger site.

    Chains are split wherever consecutive C-alpha atoms are not bonded (3.8 A apart, 2.9-4.6 A allowed). Only
    ATOM records form the C-alpha trace and supply ligands.

    Helix-turn-helix and other supersecondary motifs are not detected: asking for ``"helix_turn_helix"`` raises
    ``NotImplementedError`` rather than returning guesses.

    Args:
        structure: Structure dictionary from PDB loading, with an ``atoms`` list (atom ``name``, ``res_name``,
            ``chain_id``, ``res_seq``, ``x``, ``y``, ``z`` and, for ions, ``element``)
        motif_types: Motif types to search for; defaults to every type of ``SUPPORTED_MOTIF_TYPES``

    Returns:
        List of motifs. A helix or extended segment has ``type``, ``chain_id``, ``start_residue`` and
        ``end_residue`` (residue numbers of the first and last residue) and ``length`` (residues); they come in
        chain order. A zinc site has ``type``, ``chain_id``, ``ion_res_seq``, ``ion_serial``, ``position``,
        ``n_ligand_residues``, ``coordination_pattern``, ``zinc_finger_like`` and ``ligands`` (chain, residue,
        atom and distance in A, nearest first), and follows the other motifs. Empty when nothing is found.

    Raises:
        NotImplementedError: If ``"helix_turn_helix"`` is requested
        ValueError: If a motif type is unknown, or the structure has coordinates but no atom records

    Example:
        >>> analyze_structural_motifs({"atoms": []})
        []
    """
    requested = (motif_types,) if isinstance(motif_types, str) else tuple(motif_types or SUPPORTED_MOTIF_TYPES)
    for name in requested:
        if name == "helix_turn_helix":
            raise NotImplementedError(
                "helix_turn_helix is not implemented: it is a supersecondary motif that needs more than "
                f"local geometry to identify. Supported motif types: {', '.join(SUPPORTED_MOTIF_TYPES)}."
            )
        if name not in SUPPORTED_MOTIF_TYPES:
            raise ValueError(f"Unknown motif type {name!r}; supported: {', '.join(SUPPORTED_MOTIF_TYPES)}")

    atoms = structure.get("atoms") or []
    if not atoms:
        coordinates = structure.get("coordinates")
        if coordinates is not None and len(coordinates) > 0:
            raise ValueError(
                "structure has coordinates but no 'atoms' records: atom and residue names are needed to locate "
                "the C-alpha trace and the metal ligands"
            )
        return []

    motifs = _secondary_structure_elements(atoms, requested)
    if "zinc_site" in requested:
        motifs.extend(_zinc_sites(atoms))
    return motifs


def calculate_structural_similarity(structure1: Dict[str, Any], structure2: Dict[str, Any]) -> Dict[str, float]:
    """Calculate structural similarity between two proteins.

    Args:
        structure1: First structure dictionary
        structure2: Second structure dictionary

    Returns:
        Similarity metrics

    Example:
        >>> # Assuming structure data exists
        >>> # similarity = calculate_structural_similarity(struct1, struct2)
        >>> # isinstance(similarity, dict)
        >>> # True
    """
    coords1 = structure1.get("coordinates")
    coords2 = structure2.get("coordinates")

    if coords1 is None or coords2 is None:
        return {"rmsd": float("inf"), "similarity": 0.0}

    # Use RMSD as similarity measure
    from .general import compute_rmsd_simple

    rmsd = compute_rmsd_simple(coords1, coords2)

    # Convert RMSD to similarity score (simplified)
    # Lower RMSD = higher similarity
    if rmsd == 0:
        similarity = 1.0
    else:
        similarity = 1.0 / (1.0 + rmsd / 10.0)  # Arbitrary scaling

    return {"rmsd": rmsd, "similarity": similarity, "n_atoms_1": len(coords1), "n_atoms_2": len(coords2)}


# Van der Waals radii (A, Bondi) by element; any other element gets the carbon radius.
_VDW_RADII = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "P": 1.80, "S": 1.80, "SE": 1.90}
_DEFAULT_VDW_RADIUS = 1.70
_WATER_RESIDUES = frozenset({"HOH", "WAT", "H2O", "DOD"})

# Scan directions of the buriedness test: the three axes and the four body diagonals of the grid.
_SCAN_DIRECTIONS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1))
_MAX_GRID_CELLS = 30_000_000


def _vdw_radius(atom: Dict[str, Any]) -> float:
    """Van der Waals radius of an atom dict, by its element (or, without one, the first letter of its name)."""
    symbol = str(atom.get("element") or "").strip().upper()
    if not symbol:
        symbol = str(atom.get("name") or "").strip().lstrip("0123456789")[:1].upper()
    return _VDW_RADII.get(symbol, _DEFAULT_VDW_RADIUS)


def _mark_spheres(
    shape: Tuple[int, int, int], origin: np.ndarray, spacing: float, centers: np.ndarray, radii: np.ndarray
) -> np.ndarray:
    """Boolean grid, True where a cell centre (``origin + spacing * index``) lies within ``radii`` of a centre."""
    marked = np.zeros(shape, dtype=bool)
    limit = np.array(shape)
    for center, radius in zip(centers, radii):
        low = np.maximum(np.floor((center - radius - origin) / spacing).astype(int), 0)
        high = np.minimum(np.ceil((center + radius - origin) / spacing).astype(int) + 1, limit)
        if np.any(high <= low):
            continue
        offsets = [(origin[axis] + spacing * np.arange(low[axis], high[axis]) - center[axis]) ** 2 for axis in range(3)]
        inside = offsets[0][:, None, None] + offsets[1][None, :, None] + offsets[2][None, None, :] <= radius**2
        marked[low[0] : high[0], low[1] : high[1], low[2] : high[2]] |= inside
    return marked


def _buriedness(occupied: np.ndarray, spacing: float, ray_length: float) -> np.ndarray:
    """For every grid cell, in how many of the seven scan directions it lies between occupied cells.

    A direction counts when the line through the cell meets an occupied cell within ``ray_length`` in both
    senses (a protein-solvent-protein event).
    """
    nx, ny, nz = occupied.shape
    axis_steps = max(1, int(ray_length / spacing))
    diagonal_steps = max(1, int(ray_length / (spacing * math.sqrt(3.0))))
    pad = axis_steps
    padded = np.pad(occupied, pad, mode="constant", constant_values=False)
    count = np.zeros(occupied.shape, dtype=np.int8)
    for dx, dy, dz in _SCAN_DIRECTIONS:
        steps = axis_steps if abs(dx) + abs(dy) + abs(dz) == 1 else diagonal_steps
        ahead = np.zeros(occupied.shape, dtype=bool)
        behind = np.zeros(occupied.shape, dtype=bool)
        for step in range(1, steps + 1):
            ahead |= padded[
                pad + step * dx : pad + step * dx + nx,
                pad + step * dy : pad + step * dy + ny,
                pad + step * dz : pad + step * dz + nz,
            ]
            behind |= padded[
                pad - step * dx : pad - step * dx + nx,
                pad - step * dy : pad - step * dy + ny,
                pad - step * dz : pad - step * dz + nz,
            ]
        count += ahead & behind
    return count


def identify_ligand_binding_sites(
    structure: Dict[str, Any],
    *,
    grid_spacing: float = 1.0,
    probe_radius: float = 1.4,
    min_buriedness: int = 5,
    ray_length: float = 8.0,
    min_volume: float = 30.0,
    lining_distance: float = 4.0,
) -> List[Dict[str, Any]]:
    """Find buried pockets that could bind a ligand, by grid-based buriedness (the LIGSITE approach).

    The protein is placed on a regular grid. A cell is *occupied* when it lies within the van der Waals radius of
    an atom, and *accessible* when a probe sphere of ``probe_radius`` centred on it touches no atom. Every cell
    is then scanned along the three axes and the four body diagonals of the grid: a direction counts as buried
    when occupied cells are met on both sides of the cell within ``ray_length`` (a protein-solvent-protein
    event). Accessible cells buried in at least ``min_buriedness`` of the seven directions are pocket cells;
    connected pocket cells (26-neighbourhood) form a pocket. Pockets are ranked by volume.

    Only ATOM records are used: waters, ions and ligands (HETATM) are left out, so a pocket occupied by a bound
    ligand is still found. A structure with no pocket, such as a compact solid or a flat surface, gives an empty
    list.

    Args:
        structure: Structure dictionary from PDB loading. Its ``atoms`` supply coordinates and elements; without
            atoms the ``coordinates`` array is used with carbon radii for every point.
        grid_spacing: Grid cell size in A. Volumes are cell counts times its cube, so they are estimates that
            improve as it shrinks.
        probe_radius: Radius in A of the probe that must fit in a pocket cell (1.4 A is a water molecule)
        min_buriedness: Directions (out of 7) in which a cell must be buried to count as pocket
        ray_length: Length in A of each scan in each sense; it must reach across the pockets of interest
        min_volume: Smallest pocket volume in A^3 that is reported
        lining_distance: Atoms within this distance in A of a pocket cell line the pocket

    Returns:
        List of pockets, the largest first, each a dictionary with ``rank``, ``center`` (centroid of the pocket
        cells, in A), ``volume`` (A^3), ``n_grid_points``, ``mean_buriedness`` (mean number of buried
        directions over its cells), ``lining_atoms`` (indices into ``structure["atoms"]``, or into
        ``coordinates`` without atoms) and ``lining_residues`` (``chain_id``, ``res_seq``, ``res_name`` of the
        lining residues; empty without atoms).

    Raises:
        ValueError: If a parameter is out of range, or the grid would exceed 30 million cells (increase
            ``grid_spacing``)
    """
    if grid_spacing <= 0 or ray_length <= 0:
        raise ValueError("grid_spacing and ray_length must be positive")
    if probe_radius < 0 or min_volume < 0 or lining_distance < 0:
        raise ValueError("probe_radius, min_volume and lining_distance must not be negative")
    if not 1 <= min_buriedness <= len(_SCAN_DIRECTIONS):
        raise ValueError(f"min_buriedness must be between 1 and {len(_SCAN_DIRECTIONS)}, got {min_buriedness}")

    atoms = structure.get("atoms") or []
    if atoms:
        indices = [
            i
            for i, atom in enumerate(atoms)
            if _record_type(atom) == "ATOM" and str(atom.get("res_name", "")).strip().upper() not in _WATER_RESIDUES
        ]
        coords = np.array([[atoms[i]["x"], atoms[i]["y"], atoms[i]["z"]] for i in indices], dtype=float)
        radii = np.array([_vdw_radius(atoms[i]) for i in indices], dtype=float)
    else:
        raw = structure.get("coordinates")
        coords = np.empty((0, 3)) if raw is None else np.asarray(raw, dtype=float).reshape(-1, 3)
        indices = list(range(len(coords)))
        radii = np.full(len(coords), _DEFAULT_VDW_RADIUS)
    if len(coords) == 0:
        return []

    padding = float(radii.max()) + probe_radius + grid_spacing
    origin = coords.min(axis=0) - padding
    shape_array = np.ceil((coords.max(axis=0) + padding - origin) / grid_spacing).astype(int) + 1
    if int(np.prod(shape_array.astype(float))) > _MAX_GRID_CELLS:
        raise ValueError(
            f"a grid of {int(np.prod(shape_array.astype(float))):,} cells exceeds the limit of {_MAX_GRID_CELLS:,}; "
            "increase grid_spacing"
        )
    shape = (int(shape_array[0]), int(shape_array[1]), int(shape_array[2]))

    occupied = _mark_spheres(shape, origin, grid_spacing, coords, radii)
    inaccessible = _mark_spheres(shape, origin, grid_spacing, coords, radii + probe_radius)
    buried = _buriedness(occupied, grid_spacing, ray_length)
    pocket_cells = ~inaccessible & (buried >= min_buriedness)
    if not pocket_cells.any():
        return []

    labels, _ = ndimage.label(pocket_cells, structure=np.ones((3, 3, 3), dtype=int))
    cell_index = np.argwhere(pocket_cells)  # same C order as the boolean indexing below
    cell_label = labels[pocket_cells]
    cell_buriedness = buried[pocket_cells]
    order = np.argsort(cell_label, kind="stable")
    groups = np.split(order, np.flatnonzero(np.diff(cell_label[order])) + 1)

    pockets: List[Dict[str, Any]] = []
    for group in groups:
        volume = len(group) * grid_spacing**3
        if volume < min_volume:
            continue
        points = origin + grid_spacing * cell_index[group]
        nearest, _ = cKDTree(points).query(coords, distance_upper_bound=lining_distance)
        lining = [indices[i] for i in np.flatnonzero(np.isfinite(nearest))]
        residues = sorted(
            {
                (str(atoms[i].get("chain_id", "")), atoms[i]["res_seq"], str(atoms[i].get("res_name", "")).strip())
                for i in lining
            }
            if atoms
            else set()
        )
        pockets.append(
            {
                "center": points.mean(axis=0).tolist(),
                "volume": float(volume),
                "n_grid_points": int(len(group)),
                "mean_buriedness": float(cell_buriedness[group].mean()),
                "lining_atoms": lining,
                "lining_residues": [
                    {"chain_id": chain_id, "res_seq": number, "res_name": name} for chain_id, number, name in residues
                ],
            }
        )

    pockets.sort(key=lambda pocket: (-pocket["volume"], pocket["center"]))
    for rank, pocket in enumerate(pockets, start=1):
        pocket["rank"] = rank
    return pockets


def analyze_protein_flexibility(structure: Dict[str, Any]) -> Dict[str, Any]:
    """Analyze protein structural flexibility.

    Args:
        structure: Structure dictionary

    Returns:
        Flexibility analysis results

    Example:
        >>> # Assuming structure data exists
        >>> # flexibility = analyze_protein_flexibility(structure)
        >>> # isinstance(flexibility, dict)
        >>> # True
    """
    atoms = structure.get("atoms", [])
    coords = structure.get("coordinates")

    if coords is None:
        return {"flexibility_score": 0.0, "rigid_regions": [], "flexible_regions": []}

    # Simple flexibility analysis based on B-factors
    b_factors = []
    for atom in atoms:
        b_factor = atom.get("temp_factor", 0.0)
        b_factors.append(b_factor)

    if not b_factors:
        return {"flexibility_score": 0.0, "rigid_regions": [], "flexible_regions": []}

    avg_b_factor = np.mean(b_factors)

    # Classify regions (simplified)
    rigid_regions = []
    flexible_regions = []

    # Group by residue
    residue_b_factors: Dict[Tuple[Any, Any], List[Any]] = {}
    for atom in atoms:
        key = (atom["chain_id"], atom["res_seq"])
        b_factor = atom.get("temp_factor", 0.0)
        if key not in residue_b_factors:
            residue_b_factors[key] = []
        residue_b_factors[key].append(b_factor)

    for (chain, res), b_vals in residue_b_factors.items():
        avg_res_b = np.mean(b_vals)
        if avg_res_b < avg_b_factor * 0.8:
            rigid_regions.append({"chain": chain, "residue": res, "b_factor": avg_res_b})
        elif avg_res_b > avg_b_factor * 1.2:
            flexible_regions.append({"chain": chain, "residue": res, "b_factor": avg_res_b})

    return {
        "flexibility_score": np.std(b_factors) / np.mean(b_factors) if np.mean(b_factors) > 0 else 0.0,
        "average_b_factor": avg_b_factor,
        "rigid_regions": rigid_regions,
        "flexible_regions": flexible_regions,
    }


def calculate_structural_alignment_quality(structure1: Dict[str, Any], structure2: Dict[str, Any]) -> Dict[str, float]:
    """Assess quality of structural alignment.

    Args:
        structure1: First structure
        structure2: Second structure

    Returns:
        Alignment quality metrics

    Example:
        >>> # Assuming structure data exists
        >>> # quality = calculate_structural_alignment_quality(struct1, struct2)
        >>> # isinstance(quality, dict)
        >>> # True
    """
    coords1 = structure1.get("coordinates")
    coords2 = structure2.get("coordinates")

    if coords1 is None or coords2 is None:
        return {"quality_score": 0.0, "aligned_atoms": 0}

    # Calculate RMSD
    from .general import compute_rmsd_simple

    rmsd = compute_rmsd_simple(coords1, coords2)

    # Calculate alignment quality (simplified)
    # Lower RMSD = higher quality
    quality_score = max(0.0, 1.0 - rmsd / 20.0)  # Arbitrary scaling

    return {
        "quality_score": quality_score,
        "rmsd": rmsd,
        "aligned_atoms": min(len(coords1), len(coords2)),
        "structure1_atoms": len(coords1),
        "structure2_atoms": len(coords2),
    }
