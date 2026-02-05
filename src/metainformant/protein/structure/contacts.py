"""Protein contact analysis utilities.

This module provides functions for analyzing protein contacts, including
residue-residue contacts, hydrogen bonds, salt bridges, and hydrophobic interactions.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def calculate_residue_contacts(
    coords: np.ndarray, residue_ranges: List[Tuple[int, int]], threshold: float = 8.0
) -> np.ndarray:
    """Calculate residue-residue contact map.

    Args:
        coords: Atomic coordinates (n_atoms, 3)
        residue_ranges: List of (start_atom, end_atom) for each residue
        threshold: Distance threshold for contacts (Å)

    Returns:
        Contact map matrix (n_residues, n_residues)

    Example:
        >>> coords = np.random.rand(50, 3) * 50
        >>> residue_ranges = [(0, 5), (5, 10), (10, 15)]  # 3 residues
        >>> contacts = calculate_residue_contacts(coords, residue_ranges)
        >>> contacts.shape
        (3, 3)
    """
    n_residues = len(residue_ranges)
    contact_map = np.zeros((n_residues, n_residues))

    # Calculate minimum distance between any atom of residue i and j
    for i in range(n_residues):
        start_i, end_i = residue_ranges[i]
        coords_i = coords[start_i:end_i]

        for j in range(i + 1, n_residues):
            start_j, end_j = residue_ranges[j]
            coords_j = coords[start_j:end_j]

            # Calculate all pairwise distances
            diff = coords_i[:, np.newaxis, :] - coords_j[np.newaxis, :, :]
            distances = np.sqrt(np.sum(diff**2, axis=2))

            # Check if any atom pair is within threshold
            if np.any(distances <= threshold):
                contact_map[i, j] = 1
                contact_map[j, i] = 1

    return contact_map


def identify_hydrogen_bonds(
    atoms: List[Dict[str, Any]], coords: np.ndarray, distance_threshold: float = 3.5, angle_threshold: float = 120.0
) -> List[Dict[str, Any]]:
    """Identify hydrogen bonds in protein structure.

    Args:
        atoms: List of atom dictionaries
        coords: Atomic coordinates
        distance_threshold: Maximum H-bond distance (Å)
        angle_threshold: Minimum H-bond angle (degrees)

    Returns:
        List of hydrogen bond dictionaries

    Example:
        >>> # Assuming atom data exists
        >>> # h_bonds = identify_hydrogen_bonds(atoms, coords)
        >>> # isinstance(h_bonds, list)
        >>> # True
    """
    h_bonds = []

    # Find potential donors and acceptors
    donors = []
    acceptors = []

    for i, atom in enumerate(atoms):
        element = atom.get("element", "").upper()
        name = atom.get("name", "").upper()

        # Nitrogen atoms (potential donors)
        if element == "N":
            donors.append(i)

        # Oxygen atoms (potential acceptors)
        elif element == "O":
            acceptors.append(i)

    # Check all donor-acceptor pairs
    for donor_idx in donors:
        donor_pos = coords[donor_idx]

        for acceptor_idx in acceptors:
            if donor_idx == acceptor_idx:
                continue

            acceptor_pos = coords[acceptor_idx]
            distance = np.linalg.norm(donor_pos - acceptor_pos)

            if distance <= distance_threshold:
                # Check angle (simplified - would need hydrogen position for accurate angle)
                # For now, just use distance-based criterion
                h_bond = {
                    "donor_atom": atoms[donor_idx],
                    "acceptor_atom": atoms[acceptor_idx],
                    "distance": distance,
                    "angle": None,  # Would calculate if H positions available
                    "strength": max(0, 1 - distance / distance_threshold),
                }
                h_bonds.append(h_bond)

    return h_bonds


def identify_salt_bridges(
    atoms: List[Dict[str, Any]], coords: np.ndarray, distance_threshold: float = 4.0
) -> List[Dict[str, Any]]:
    """Identify salt bridges between charged residues.

    Args:
        atoms: List of atom dictionaries
        coords: Atomic coordinates
        distance_threshold: Maximum salt bridge distance (Å)

    Returns:
        List of salt bridge dictionaries

    Example:
        >>> # Assuming atom data exists
        >>> # salt_bridges = identify_salt_bridges(atoms, coords)
        >>> # isinstance(salt_bridges, list)
        >>> # True
    """
    salt_bridges = []

    # Define charged residue types
    positive_residues = {"ARG", "LYS", "HIS"}
    negative_residues = {"ASP", "GLU"}

    # Group atoms by residue
    residue_atoms = {}
    for i, atom in enumerate(atoms):
        res_key = (atom["chain_id"], atom["res_seq"], atom["res_name"])
        if res_key not in residue_atoms:
            residue_atoms[res_key] = []
        residue_atoms[res_key].append((i, atom))

    # Find charged groups
    positive_groups = []
    negative_groups = []

    for res_key, atom_list in residue_atoms.items():
        chain_id, res_seq, res_name = res_key

        if res_name in positive_residues:
            # Find charged atoms (NZ for LYS, NH1/NH2 for ARG, ND1/NE2 for HIS)
            charged_atoms = []
            for atom_idx, atom in atom_list:
                if res_name == "LYS" and atom["name"] == "NZ":
                    charged_atoms.append(atom_idx)
                elif res_name == "ARG" and atom["name"] in ["NH1", "NH2"]:
                    charged_atoms.append(atom_idx)
                elif res_name == "HIS" and atom["name"] in ["ND1", "NE2"]:
                    charged_atoms.append(atom_idx)

            if charged_atoms:
                positive_groups.append(
                    {"residue": res_key, "atoms": charged_atoms, "center": np.mean(coords[charged_atoms], axis=0)}
                )

        elif res_name in negative_residues:
            # Find charged atoms (OD1/OD2 for ASP, OE1/OE2 for GLU)
            charged_atoms = []
            for atom_idx, atom in atom_list:
                if res_name == "ASP" and atom["name"] in ["OD1", "OD2"]:
                    charged_atoms.append(atom_idx)
                elif res_name == "GLU" and atom["name"] in ["OE1", "OE2"]:
                    charged_atoms.append(atom_idx)

            if charged_atoms:
                negative_groups.append(
                    {"residue": res_key, "atoms": charged_atoms, "center": np.mean(coords[charged_atoms], axis=0)}
                )

    # Check all positive-negative pairs
    for pos_group in positive_groups:
        for neg_group in negative_groups:
            distance = np.linalg.norm(pos_group["center"] - neg_group["center"])

            if distance <= distance_threshold:
                salt_bridge = {
                    "positive_residue": pos_group["residue"],
                    "negative_residue": neg_group["residue"],
                    "distance": distance,
                    "positive_atoms": pos_group["atoms"],
                    "negative_atoms": neg_group["atoms"],
                }
                salt_bridges.append(salt_bridge)

    return salt_bridges


def identify_hydrophobic_contacts(
    atoms: List[Dict[str, Any]], coords: np.ndarray, distance_threshold: float = 5.0
) -> List[Dict[str, Any]]:
    """Identify hydrophobic contacts between non-polar residues.

    Args:
        atoms: List of atom dictionaries
        coords: Atomic coordinates
        distance_threshold: Maximum hydrophobic contact distance (Å)

    Returns:
        List of hydrophobic contact dictionaries

    Example:
        >>> # Assuming atom data exists
        >>> # contacts = identify_hydrophobic_contacts(atoms, coords)
        >>> # isinstance(contacts, list)
        >>> # True
    """
    hydrophobic_contacts = []

    # Define hydrophobic residue types
    hydrophobic_residues = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}

    # Group atoms by residue
    residue_atoms = {}
    for i, atom in enumerate(atoms):
        res_key = (atom["chain_id"], atom["res_seq"], atom["res_name"])
        if res_key not in residue_atoms:
            residue_atoms[res_key] = []
        residue_atoms[res_key].append((i, atom))

    # Find hydrophobic groups (use CB atoms or CA for GLY)
    hydrophobic_groups = []

    for res_key, atom_list in residue_atoms.items():
        chain_id, res_seq, res_name = res_key

        if res_name in hydrophobic_residues:
            # Find hydrophobic atoms (CB for most, CA for GLY)
            hydrophobic_atoms = []
            for atom_idx, atom in atom_list:
                if atom["name"] in ["CB", "CG", "CG1", "CG2", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]:
                    hydrophobic_atoms.append(atom_idx)
                elif res_name == "GLY" and atom["name"] == "CA":
                    hydrophobic_atoms.append(atom_idx)

            if hydrophobic_atoms:
                hydrophobic_groups.append(
                    {
                        "residue": res_key,
                        "atoms": hydrophobic_atoms,
                        "center": np.mean(coords[hydrophobic_atoms], axis=0),
                    }
                )

    # Check all hydrophobic pairs
    for i, group1 in enumerate(hydrophobic_groups):
        for group2 in hydrophobic_groups[i + 1 :]:
            distance = np.linalg.norm(group1["center"] - group2["center"])

            if distance <= distance_threshold:
                contact = {
                    "residue1": group1["residue"],
                    "residue2": group2["residue"],
                    "distance": distance,
                    "atoms1": group1["atoms"],
                    "atoms2": group2["atoms"],
                }
                hydrophobic_contacts.append(contact)

    return hydrophobic_contacts


def analyze_contact_network(contact_map: np.ndarray) -> Dict[str, Any]:
    """Analyze the protein contact network.

    Args:
        contact_map: Residue contact map

    Returns:
        Network analysis results

    Example:
        >>> contact_map = np.random.randint(0, 2, (10, 10))
        >>> analysis = analyze_contact_network(contact_map)
        >>> isinstance(analysis, dict)
        True
    """
    n_residues = contact_map.shape[0]

    # Calculate network properties
    degree = np.sum(contact_map, axis=1)
    avg_degree = np.mean(degree)

    # Clustering coefficient (simplified)
    triangles = 0
    for i in range(n_residues):
        neighbors = np.where(contact_map[i])[0]
        if len(neighbors) < 2:
            continue

        # Count triangles
        for j in neighbors:
            for k in neighbors:
                if j < k and contact_map[j, k]:
                    triangles += 1

    clustering_coeff = triangles / max(1, np.sum(degree * (degree - 1) / 2))

    # Connected components
    visited = np.zeros(n_residues, dtype=bool)
    components = []

    for i in range(n_residues):
        if not visited[i]:
            component = []
            stack = [i]
            while stack:
                node = stack.pop()
                if not visited[node]:
                    visited[node] = True
                    component.append(node)
                    neighbors = np.where(contact_map[node])[0]
                    stack.extend(neighbors)
            components.append(component)

    return {
        "n_residues": n_residues,
        "n_contacts": int(np.sum(contact_map) / 2),  # Undirected
        "average_degree": avg_degree,
        "max_degree": np.max(degree),
        "clustering_coefficient": clustering_coeff,
        "n_components": len(components),
        "component_sizes": [len(comp) for comp in components],
        "largest_component_size": max(len(comp) for comp in components) if components else 0,
    }


def calculate_contact_persistence(contact_maps: List[np.ndarray]) -> np.ndarray:
    """Calculate contact persistence across multiple structures.

    Args:
        contact_maps: List of contact maps from different structures

    Returns:
        Persistence matrix showing how often each contact appears

    Example:
        >>> maps = [np.random.randint(0, 2, (5, 5)) for _ in range(3)]
        >>> persistence = calculate_contact_persistence(maps)
        >>> persistence.shape
        (5, 5)
    """
    if not contact_maps:
        return np.array([])

    n_residues = contact_maps[0].shape[0]
    persistence = np.zeros((n_residues, n_residues))

    for contact_map in contact_maps:
        persistence += contact_map

    # Normalize by number of structures
    persistence /= len(contact_maps)

    return persistence


def identify_disulfide_bonds(
    atoms: List[Dict[str, Any]], coords: np.ndarray, distance_threshold: float = 2.5
) -> List[Dict[str, Any]]:
    """Identify disulfide bonds between cysteine residues.

    Args:
        atoms: List of atom dictionaries
        coords: Atomic coordinates
        distance_threshold: Maximum S-S distance for disulfide bond (Å)

    Returns:
        List of disulfide bond dictionaries

    Example:
        >>> # Assuming atom data exists
        >>> # disulfides = identify_disulfide_bonds(atoms, coords)
        >>> # isinstance(disulfides, list)
        >>> # True
    """
    disulfides = []

    # Find cysteine SG atoms
    cysteines = []

    for i, atom in enumerate(atoms):
        if atom["res_name"] == "CYS" and atom["name"] == "SG":
            cysteines.append({"atom_idx": i, "residue": (atom["chain_id"], atom["res_seq"]), "position": coords[i]})

    # Check all cysteine pairs
    for i, cys1 in enumerate(cysteines):
        for cys2 in cysteines[i + 1 :]:
            distance = np.linalg.norm(cys1["position"] - cys2["position"])

            if distance <= distance_threshold:
                disulfide = {
                    "residue1": cys1["residue"],
                    "residue2": cys2["residue"],
                    "distance": distance,
                    "atom1_idx": cys1["atom_idx"],
                    "atom2_idx": cys2["atom_idx"],
                }
                disulfides.append(disulfide)

    return disulfides


def classify_contact_types(atoms: List[Dict[str, Any]], coords: np.ndarray) -> Dict[str, List[Dict[str, Any]]]:
    """Classify all types of contacts in a protein structure.

    Args:
        atoms: List of atom dictionaries
        coords: Atomic coordinates

    Returns:
        Dictionary mapping contact types to contact lists

    Example:
        >>> # Assuming atom data exists
        >>> # contacts = classify_contact_types(atoms, coords)
        >>> # isinstance(contacts, dict)
        >>> # True
    """
    # Get residue ranges (simplified)
    residue_ranges = []
    current_residue = None
    start_idx = 0

    for i, atom in enumerate(atoms):
        res_key = (atom["chain_id"], atom["res_seq"])
        if current_residue != res_key:
            if current_residue is not None:
                residue_ranges.append((start_idx, i))
            current_residue = res_key
            start_idx = i

    if current_residue is not None:
        residue_ranges.append((start_idx, len(atoms)))

    # Calculate contact map
    contact_map = calculate_residue_contacts(coords, residue_ranges)

    # Identify specific contact types
    h_bonds = identify_hydrogen_bonds(atoms, coords)
    salt_bridges = identify_salt_bridges(atoms, coords)
    hydrophobic = identify_hydrophobic_contacts(atoms, coords)
    disulfides = identify_disulfide_bonds(atoms, coords)

    return {
        "residue_contacts": contact_map,
        "hydrogen_bonds": h_bonds,
        "salt_bridges": salt_bridges,
        "hydrophobic_contacts": hydrophobic,
        "disulfide_bonds": disulfides,
        "summary": {
            "total_h_bonds": len(h_bonds),
            "total_salt_bridges": len(salt_bridges),
            "total_hydrophobic": len(hydrophobic),
            "total_disulfides": len(disulfides),
            "total_residue_contacts": int(np.sum(contact_map) / 2),
        },
    }


def compute_ca_contact_pairs(
    coords: list,
    threshold: float = 8.0,
) -> list:
    """Compute C-alpha contact pairs from a list of coordinates.

    Given a list of (x, y, z) tuples representing C-alpha positions,
    returns all pairs (i, j) with i < j whose Euclidean distance
    is within the threshold.

    Uses vectorized numpy distance computation for O(n^2) memory but
    fast execution, with scipy.spatial.distance.pdist when available.

    Args:
        coords: List of (x, y, z) tuples or lists for C-alpha atoms
        threshold: Distance threshold in Angstroms

    Returns:
        List of (i, j) tuples representing contact pairs

    Example:
        >>> coords = [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (10.0, 0.0, 0.0)]
        >>> pairs = compute_ca_contact_pairs(coords, threshold=4.0)
        >>> (0, 1) in pairs
        True
        >>> (0, 2) in pairs
        False
    """
    n = len(coords)
    if n < 2:
        return []

    arr = np.asarray(coords, dtype=np.float64)

    try:
        from scipy.spatial.distance import pdist, squareform

        condensed = pdist(arr)
        # Convert condensed index back to (i, j) pairs
        pairs = []
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if condensed[k] <= threshold:
                    pairs.append((i, j))
                k += 1
        return pairs
    except ImportError:
        pass

    # Fallback: numpy broadcasting
    diff = arr[:, np.newaxis, :] - arr[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff ** 2, axis=2))

    # Extract upper triangle where distance <= threshold
    i_idx, j_idx = np.where(np.triu(dist_matrix <= threshold, k=1))
    return list(zip(i_idx.tolist(), j_idx.tolist()))
