"""Protein structure analysis utilities.

This module provides tools for analyzing protein 3D structures,
including domain identification, surface analysis, and structural motifs.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy import ndimage
from scipy.sparse import csr_matrix, diags, triu
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh
from scipy.spatial import cKDTree

from metainformant.core.utils import logging

from .general import (
    _CA_BOND_RANGE,
    DEFAULT_SASA_POINTS,
    VAN_DER_WAALS_RADII,
    WATER_PROBE_RADIUS,
    _as_points,
    _dihedral_degrees,
    calculate_atomic_sasa,
    calculate_solvent_accessible_surface_area,
    calculate_tm_score,
    compute_rmsd_kabsch,
    van_der_waals_radius,
)

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


# ---------------------------------------------------------------------------------------------------------------------
# Domain identification from the residue contact graph
# ---------------------------------------------------------------------------------------------------------------------

#: Default C-alpha distance in angstrom within which two residues are in contact (the default of
#: :func:`calculate_contact_map`)
DOMAIN_CONTACT_DISTANCE = 8.0

#: Default smallest number of residues of a domain
DOMAIN_MIN_SIZE = 40

#: Default largest normalised cut (see :func:`identify_domains`) at which a bisection is accepted as a domain boundary
DOMAIN_MAX_NORMALIZED_CUT = 0.10

# Graphs up to this many residues are decomposed with a dense eigensolver, larger ones with ARPACK
_DENSE_EIGEN_LIMIT = 1000


def _record_type(atom: Dict[str, Any]) -> str:
    """Record type (``ATOM`` or ``HETATM``) of an atom dict from either PDB loader of this package."""
    return str(atom.get("record_name") or atom.get("record_type") or "ATOM").strip().upper()


def _insertion_code(atom: Dict[str, Any]) -> str:
    return str(atom.get("icode") or atom.get("ins_code") or "").strip()


def _ca_residues(atoms: Sequence[Dict[str, Any]]) -> List[Tuple[Tuple[str, Any, str], np.ndarray]]:
    """The C-alpha atom of every protein residue, in file order, as ``((chain, number, insertion code), position)``.

    Only ATOM records count, so a calcium ion (atom name ``CA`` in a HETATM record) is not a C-alpha, and only the
    first alternate location of a residue is used.

    Raises:
        ValueError: If a C-alpha atom has no ``x``, ``y`` or ``z`` coordinate
    """
    seen = set()
    residues: List[Tuple[Tuple[str, Any, str], np.ndarray]] = []
    for atom in atoms:
        if _record_type(atom) != "ATOM" or str(atom.get("name", "")).strip().upper() != "CA":
            continue
        key = (str(atom.get("chain_id", "")), atom["res_seq"], _insertion_code(atom))
        if key in seen:
            continue
        seen.add(key)
        try:
            position = np.array([atom["x"], atom["y"], atom["z"]], dtype=float)
        except KeyError as error:
            raise ValueError(
                f"atom records need x, y and z coordinates; residue {key} has no {error.args[0]!r}"
            ) from None
        residues.append((key, position))
    return residues


def _residue_contact_graph(positions: np.ndarray, contact_distance: float) -> csr_matrix:
    """Symmetric 0/1 adjacency matrix of the residues whose C-alpha atoms are at most ``contact_distance`` apart."""
    n_residues = len(positions)
    pairs = cKDTree(positions).query_pairs(contact_distance, output_type="ndarray")
    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    columns = np.concatenate([pairs[:, 1], pairs[:, 0]])
    return csr_matrix((np.ones(len(rows)), (rows, columns)), shape=(n_residues, n_residues))


def _spectral_order(adjacency: csr_matrix) -> np.ndarray:
    """Residue order along the Fiedler vector of a connected contact graph (spectral ordering).

    The Fiedler vector is the eigenvector of the second smallest eigenvalue of the normalised Laplacian
    ``I - D^-1/2 A D^-1/2``, scaled by ``D^-1/2``; it varies slowly inside a densely connected group of residues and
    jumps across a narrow connection, so sorting residues by it puts the two sides of a bottleneck on opposite ends.
    The eigensolver starts from a fixed vector, so the result is deterministic.
    """
    n_residues = adjacency.shape[0]
    inverse_root_degree = 1.0 / np.sqrt(np.asarray(adjacency.sum(axis=1)).ravel())
    if n_residues <= _DENSE_EIGEN_LIMIT:
        scaled = adjacency.toarray() * inverse_root_degree[:, np.newaxis] * inverse_root_degree[np.newaxis, :]
        _, vectors = np.linalg.eigh(scaled)
        second = vectors[:, -2]
    else:
        scaling = diags(inverse_root_degree)
        start = np.random.default_rng(0).standard_normal(n_residues)
        _, vectors = eigsh(scaling @ adjacency @ scaling, k=2, which="LA", v0=start)
        second = vectors[:, 0]  # eigsh returns the eigenvalues in ascending order; the first is the second largest
    return np.asarray(np.lexsort((np.arange(n_residues), second * inverse_root_degree)))


def _best_normalized_cut(adjacency: csr_matrix, order: np.ndarray, min_size: int) -> Tuple[int, float]:
    """Size ``k`` of the first side of the best split of ``order`` into two prefixes, and its normalised cut.

    The normalised cut of a split into ``S`` and the rest is ``cut / vol(S) + cut / vol(rest)``, with ``cut`` the number
    of contacts between the two sides and ``vol`` the total number of contacts of the residues of a side (Shi & Malik
    2000). Both sides must hold at least ``min_size`` residues; ``(0, inf)`` is returned when no split qualifies.
    """
    n_residues = len(order)
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    rank = np.empty(n_residues, dtype=int)
    rank[order] = np.arange(n_residues)
    edges = triu(adjacency, k=1).tocoo()
    internal = np.cumsum(np.bincount(np.maximum(rank[edges.row], rank[edges.col]), minlength=n_residues))
    volume = np.cumsum(degree[order])
    cut = volume - 2.0 * internal
    sizes = np.arange(1, n_residues + 1)
    valid = (sizes >= min_size) & (n_residues - sizes >= min_size) & (volume > 0) & (volume < volume[-1])
    if not valid.any():
        return 0, float("inf")
    normalized = np.full(n_residues, np.inf)
    normalized[valid] = cut[valid] / volume[valid] + cut[valid] / (volume[-1] - volume[valid])
    best = int(np.argmin(normalized))
    return best + 1, float(normalized[best])


def _decompose_contact_graph(
    members: np.ndarray, adjacency: csr_matrix, positions: np.ndarray, min_size: int, max_cut: float
) -> List[np.ndarray]:
    """Split a set of residues (indices into ``adjacency``) into domains by recursive spectral bisection.

    Residues that are not connected through contacts form separate pieces. A piece of at least ``min_size`` residues is
    bisected further while the best bisection costs at most ``max_cut``; a smaller piece is no domain and joins the
    domain that holds the residue closest to it in space.
    """
    if len(members) < 2 * min_size:
        return [members]
    subgraph = adjacency[members][:, members]
    n_pieces, labels = connected_components(subgraph, directed=False)
    if n_pieces > 1:
        pieces = [members[labels == label] for label in range(n_pieces)]
        large = [piece for piece in pieces if len(piece) >= min_size]
        if not large:
            return [members]
        domains = [
            domain
            for piece in large
            for domain in _decompose_contact_graph(piece, adjacency, positions, min_size, max_cut)
        ]
        owner = np.concatenate([np.full(len(domain), label) for label, domain in enumerate(domains)])
        placed = np.concatenate(domains)
        nearest_tree = cKDTree(positions[placed])
        for piece in (piece for piece in pieces if len(piece) < min_size):
            distance, closest = nearest_tree.query(positions[piece])
            label = int(owner[closest[int(np.argmin(distance))]])
            domains[label] = np.sort(np.concatenate([domains[label], piece]))
        return domains

    order = _spectral_order(subgraph)
    size, normalized_cut = _best_normalized_cut(subgraph, order, min_size)
    if size == 0 or normalized_cut > max_cut:
        return [members]
    first, second = np.sort(members[order[:size]]), np.sort(members[order[size:]])
    return _decompose_contact_graph(first, adjacency, positions, min_size, max_cut) + _decompose_contact_graph(
        second, adjacency, positions, min_size, max_cut
    )


def identify_domains(
    structure: Dict[str, Any],
    *,
    contact_distance: float = DOMAIN_CONTACT_DISTANCE,
    min_domain_size: int = DOMAIN_MIN_SIZE,
    max_normalized_cut: float = DOMAIN_MAX_NORMALIZED_CUT,
) -> List[Dict[str, Any]]:
    """Identify structural domains from the residue contact graph by recursive spectral bisection.

    Each chain is treated on its own. Its residues are the nodes of a graph in which two residues are linked when
    their C-alpha atoms are at most ``contact_distance`` apart. A domain is a group of residues that is far more
    densely connected inside than to the rest, so the method looks for the cheapest way to cut the graph in two:

    1. Residues that are not connected through contacts at all are separate pieces.
    2. A connected set of residues is ordered by the Fiedler vector of its normalised graph Laplacian (spectral
       ordering). Of the cuts that split this order into a first part and the rest, each part with at least
       ``min_domain_size`` residues, the one with the smallest normalised cut ``cut / vol(A) + cut / vol(B)`` is
       taken (Shi & Malik 2000), where ``cut`` counts the contacts between the parts and ``vol`` the contacts of a
       part's residues.
    3. The cut is accepted when its normalised cut is at most ``max_normalized_cut``, and both parts are then
       treated in the same way; otherwise the set is one domain.
    4. A piece of fewer than ``min_domain_size`` residues that is disconnected from the rest is no domain; it joins the
       domain holding the residue nearest to it.

    Domains need not be continuous in sequence: a domain inserted into another splits it into two ``segments``.
    The result depends on the coordinates only; residue numbers are used for labelling.

    The parameters are not fitted to a domain database. With the defaults, planted compact globules of 90-250
    residues keep a best bisection with a normalised cut of 0.18-0.32 and stay whole, whereas two compact groups of
    80 residues joined by a linker of 6-14 residues cost 0.002-0.005 and are split. Elongated single-domain
    structures (long helices, coiled coils, solenoids) have a cheap bisection, and are split too.

    Args:
        structure: Structure dictionary from PDB loading. Only ATOM records with the atom name ``CA`` (the first
            alternate location of each residue) are used, so every atom needs ``chain_id``, ``res_seq``, ``name``,
            ``x``, ``y`` and ``z``; the keys of both PDB loaders of this package work.
        contact_distance: C-alpha distance in angstrom up to which two residues are in contact
        min_domain_size: Smallest number of residues of a domain (also the smallest part of an accepted cut)
        max_normalized_cut: Largest normalised cut, between 0 and 2, at which a cut is accepted; lower values
            split less readily

    Returns:
        Domains in order of their first residue, each a dictionary with ``domain_id`` (1, 2, ... over the structure),
        ``chain_id``, ``residues`` (residue numbers in chain order), ``residue_ids`` (the same with insertion codes,
        such as ``"52A"``), ``start`` and ``end`` (lowest and highest residue number), ``size``, ``segments``
        (``{"start", "end"}`` residue numbers of each stretch of consecutive residues of the chain that belongs to the
        domain) and ``n_contacts`` (contacts between residues of the domain). A structure without atoms, or without
        protein residues, has no domains.

    Raises:
        ValueError: If a parameter is out of range, atoms lack coordinates, or the structure has coordinates but no
            ``atoms`` records

    Reference:
        Shi J, Malik J. Normalized cuts and image segmentation. IEEE Trans Pattern Anal Mach Intell 2000;22:888-905.
        The use of the Fiedler vector for graph bisection follows Fiedler (1973) and Pothen, Simon & Liou (1990).

    Example:
        >>> identify_domains({"atoms": []})
        []
    """
    if contact_distance <= 0 or not np.isfinite(contact_distance):
        raise ValueError(f"contact_distance must be positive, got {contact_distance!r}")
    if int(min_domain_size) != min_domain_size or min_domain_size < 2:
        raise ValueError(f"min_domain_size must be an integer of at least 2, got {min_domain_size!r}")
    if not 0 <= max_normalized_cut <= 2:
        raise ValueError(f"max_normalized_cut must lie in [0, 2], got {max_normalized_cut!r}")

    atoms = structure.get("atoms") or []
    if not atoms:
        coordinates = structure.get("coordinates")
        if coordinates is not None and len(coordinates) > 0:
            raise ValueError(
                "structure has coordinates but no 'atoms' records: chain and residue labels are needed to "
                "define domains"
            )
        return []

    chains: Dict[str, List[Tuple[Tuple[str, Any, str], np.ndarray]]] = {}
    for residue in _ca_residues(atoms):
        chains.setdefault(residue[0][0], []).append(residue)

    domains: List[Dict[str, Any]] = []
    for chain_id, residues in chains.items():
        positions = np.array([position for _, position in residues])
        adjacency = _residue_contact_graph(positions, contact_distance)
        parts = _decompose_contact_graph(
            np.arange(len(residues)), adjacency, positions, int(min_domain_size), float(max_normalized_cut)
        )
        for members in sorted((np.sort(part) for part in parts), key=lambda part: int(part[0])):
            numbers = [residues[i][0][1] for i in members]
            labels = [f"{residues[i][0][1]}{residues[i][0][2]}" for i in members]
            breaks = np.flatnonzero(np.diff(members) != 1) + 1
            segments = [
                {"start": numbers[int(chunk[0])], "end": numbers[int(chunk[-1])]}
                for chunk in np.split(np.arange(len(members)), breaks)
            ]
            domains.append(
                {
                    "domain_id": len(domains) + 1,
                    "chain_id": chain_id,
                    "residues": numbers,
                    "residue_ids": labels,
                    "start": min(numbers),
                    "end": max(numbers),
                    "size": len(numbers),
                    "segments": segments,
                    "n_contacts": int(round(adjacency[members][:, members].sum() / 2)),
                }
            )
    return domains


# ---------------------------------------------------------------------------------------------------------------------
# Solvent accessible surface area
# ---------------------------------------------------------------------------------------------------------------------


def calculate_surface_area(
    coords: np.ndarray,
    probe_radius: float = WATER_PROBE_RADIUS,
    radii: Union[None, float, Sequence[float], np.ndarray] = None,
    n_points: int = DEFAULT_SASA_POINTS,
) -> float:
    """Calculate the solvent accessible surface area of a set of atoms (Shrake-Rupley).

    This is :func:`metainformant.protein.structure.general.calculate_solvent_accessible_surface_area`: every atom is a
    sphere of its van der Waals radius plus the probe radius, ``n_points`` test points per sphere are checked against
    all other spheres, and the area is the accessible fraction of the sphere areas. An isolated atom of radius ``r``
    has exactly ``4 pi (r + probe_radius)^2``.

    Args:
        coords: Atom centres (n_atoms, 3) in angstrom
        probe_radius: Solvent probe radius (Å), 1.4 for water
        radii: Van der Waals radii in angstrom, one per atom or one for all; None gives every atom the carbon radius
            (1.70 Å), since bare coordinates carry no elements. Use :func:`calculate_structure_sasa` for a structure.
        n_points: Test points per atomic sphere; more points approximate the exact area more closely

    Returns:
        Surface area in Å² (0.0 for no atoms)

    Example:
        >>> single_carbon = np.zeros((1, 3))
        >>> abs(calculate_surface_area(single_carbon) - 4 * math.pi * (1.70 + 1.4) ** 2) < 1e-9
        True
    """
    return calculate_solvent_accessible_surface_area(coords, probe_radius=probe_radius, radii=radii, n_points=n_points)


def calculate_structure_sasa(
    structure: Dict[str, Any],
    *,
    probe_radius: float = WATER_PROBE_RADIUS,
    n_points: int = DEFAULT_SASA_POINTS,
    include_hetero: bool = False,
) -> Dict[str, Any]:
    """Solvent accessible surface area of a structure by the Shrake-Rupley algorithm, per atom and per residue.

    Radii are the van der Waals radii of the atoms' elements (``VAN_DER_WAALS_RADII`` of
    :mod:`metainformant.protein.structure.general`, Bondi 1964). The element is the ``element`` field of the atom;
    when it is empty the first letter of the atom name is used for ATOM records (``CA`` of a protein residue is
    carbon), while a HETATM record without an element is an error. Only the first alternate location of each atom is
    used. The area of an atom counts the other atoms present, so the result for ``include_hetero=False`` is the
    surface of the protein alone.

    Args:
        structure: Structure dictionary from PDB loading, with ``atoms`` (``name``, ``res_name``, ``chain_id``,
            ``res_seq``, ``x``, ``y``, ``z`` and normally ``element``)
        probe_radius: Probe radius in angstrom (1.4 for water)
        n_points: Test points per atomic sphere
        include_hetero: Whether HETATM records (waters, ions, ligands) take part; by default only ATOM records do

    Returns:
        Dictionary with ``total`` (Å²), ``per_atom`` (Å² of each used atom), ``atom_indices`` (indices of those atoms
        in ``structure["atoms"]``), ``per_residue`` (one ``{"chain_id", "res_seq", "res_name", "n_atoms", "sasa"}`` per
        residue in order of appearance), ``n_atoms``, ``probe_radius`` and ``n_points``

    Raises:
        ValueError: If a parameter is out of range, atoms lack coordinates, an element has no known radius, or the
            structure has coordinates but no ``atoms`` records

    Example:
        >>> atoms = [{"name": "C", "res_name": "ALA", "chain_id": "A", "res_seq": 1, "x": 0.0, "y": 0.0, "z": 0.0}]
        >>> result = calculate_structure_sasa({"atoms": atoms})
        >>> abs(result["total"] - 4 * math.pi * (1.70 + 1.4) ** 2) < 1e-9
        True
    """
    atoms = structure.get("atoms") or []
    if not atoms:
        coordinates = structure.get("coordinates")
        if coordinates is not None and len(coordinates) > 0:
            raise ValueError("structure has coordinates but no 'atoms' records: the elements are needed for the radii")
    seen = set()
    used: List[int] = []
    radii: List[float] = []
    for index, atom in enumerate(atoms):
        is_atom_record = _record_type(atom) == "ATOM"
        if not (is_atom_record or include_hetero):
            continue
        key = (
            str(atom.get("chain_id", "")),
            atom.get("res_seq"),
            _insertion_code(atom),
            str(atom.get("name", "")).strip(),
        )
        if key in seen:
            continue
        seen.add(key)
        element = str(atom.get("element") or "").strip()
        if not element:
            if not is_atom_record:
                raise ValueError(f"HETATM {key[3]!r} of residue {key[:3]} has no element, so its radius is unknown")
            element = key[3].lstrip("0123456789")[:1]
        used.append(index)
        radii.append(van_der_waals_radius(element))

    try:
        coords = np.array([[atoms[i]["x"], atoms[i]["y"], atoms[i]["z"]] for i in used], dtype=float).reshape(-1, 3)
    except KeyError as error:
        raise ValueError(f"atom records need x, y and z coordinates; an atom has no {error.args[0]!r}") from None
    areas = calculate_atomic_sasa(coords, radii=radii or None, probe_radius=probe_radius, n_points=n_points)

    per_residue: Dict[Tuple[str, Any, str], Dict[str, Any]] = {}
    for index, area in zip(used, areas):
        atom = atoms[index]
        residue = (str(atom.get("chain_id", "")), atom.get("res_seq"), _insertion_code(atom))
        entry = per_residue.setdefault(
            residue,
            {
                "chain_id": residue[0],
                "res_seq": residue[1],
                "res_name": str(atom.get("res_name", "")).strip(),
                "n_atoms": 0,
                "sasa": 0.0,
            },
        )
        entry["n_atoms"] += 1
        entry["sasa"] += float(area)
    return {
        "total": math.fsum(areas),
        "per_atom": [float(area) for area in areas],
        "atom_indices": used,
        "per_residue": list(per_residue.values()),
        "n_atoms": len(used),
        "probe_radius": float(probe_radius),
        "n_points": int(n_points),
    }


#: Motif types :func:`analyze_structural_motifs` searches for by default.
SUPPORTED_MOTIF_TYPES: Tuple[str, ...] = ("helix", "extended", "zinc_site")

# C-alpha trace geometry. An ideal right-handed alpha-helix (radius 2.3 A, rise 1.5 A, 100 degrees per residue) has
# CA(i)-CA(i+3) = 5.05 A, CA(i)-CA(i+4) = 6.20 A and a virtual torsion of +50 degrees over four consecutive C-alpha
# atoms; an extended strand has CA(i)-CA(i+2) of about 6.6 A and CA(i)-CA(i+3) of about 10 A. The windows below
# are those ideal values with a margin that admits real, distorted elements (distances in A, torsions in degrees).
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


def _in_window(values: np.ndarray, window: Tuple[float, float]) -> np.ndarray:
    return (values >= window[0]) & (values <= window[1])


def _ca_segments(atoms: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """C-alpha atoms of the protein residues, per chain and split wherever consecutive C-alpha atoms are not bonded.

    Only ATOM records count, so a calcium ion (atom name ``CA`` in a HETATM record) is not a C-alpha, and only
    the first alternate location of a residue is used.
    """
    chains: Dict[str, List[Tuple[int, np.ndarray]]] = {}
    for (chain_id, number, _), position in _ca_residues(atoms):
        chains.setdefault(chain_id, []).append((number, position))

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


_PAIRING_MODES = ("residue_id", "order")


def _ca_trace(structure: Dict[str, Any], name: str) -> Tuple[List[Tuple[str, Any, str]], np.ndarray]:
    """Residue labels and C-alpha positions of a structure: its C-alpha ATOM records, or else its coordinates.

    A structure without ``atoms`` records but with ``coordinates`` is read as a C-alpha trace: one point per residue,
    labelled by its position.

    Raises:
        ValueError: If the structure has no coordinates, or atom records but no C-alpha atoms
    """
    atoms = structure.get("atoms") or []
    if atoms:
        residues = _ca_residues(atoms)
        if not residues:
            raise ValueError(f"{name} has atom records but no C-alpha ATOM records to compare")
        return [key for key, _ in residues], np.array([position for _, position in residues])
    coordinates = structure.get("coordinates")
    if coordinates is None or len(coordinates) == 0:
        raise ValueError(f"{name} has no coordinates to compare")
    points = _as_points(coordinates, f"{name} coordinates")
    return [("", index, "") for index in range(len(points))], points


def match_ca_atoms(
    structure1: Dict[str, Any], structure2: Dict[str, Any], pair_by: str = "residue_id"
) -> Dict[str, Any]:
    """Pair the C-alpha atoms of two structures residue by residue.

    Only C-alpha atoms of ATOM records are used (a structure without ``atoms`` but with ``coordinates`` is read as a
    C-alpha trace of one point per residue).

    Args:
        structure1: Reference structure dictionary
        structure2: Structure dictionary to pair with it
        pair_by: ``"residue_id"`` pairs residues with the same chain, residue number and insertion code and leaves
            residues present in only one structure out; ``"order"`` pairs the i-th C-alpha atom of one with the i-th
            of the other, which needs equally many

    Returns:
        Dictionary with ``reference`` and ``mobile`` (the coordinates of the paired atoms, shape (n_pairs, 3), in
        the order of ``structure1``), ``labels`` (the (chain, number, insertion code) of each reference atom),
        ``n_residues_1`` and ``n_residues_2`` (C-alpha atoms in each structure)

    Raises:
        ValueError: If a structure has no coordinates or no C-alpha atoms, no residue is shared (``"residue_id"``),
            the lengths differ (``"order"``) or ``pair_by`` is unknown
    """
    if pair_by not in _PAIRING_MODES:
        raise ValueError(f"pair_by must be one of {', '.join(_PAIRING_MODES)}, got {pair_by!r}")
    labels1, trace1 = _ca_trace(structure1, "structure1")
    labels2, trace2 = _ca_trace(structure2, "structure2")

    if pair_by == "order":
        if len(trace1) != len(trace2):
            raise ValueError(
                f"pairing by order needs equally long C-alpha traces, got {len(trace1)} and {len(trace2)} residues"
            )
        first = np.arange(len(trace1))
        second = first
    else:
        position2 = {label: index for index, label in enumerate(labels2)}
        matched = [(index, position2[label]) for index, label in enumerate(labels1) if label in position2]
        if not matched:
            raise ValueError(
                "the structures have no residue in common (matched by chain, residue number and insertion code)"
            )
        first = np.array([pair[0] for pair in matched])
        second = np.array([pair[1] for pair in matched])

    return {
        "reference": trace1[first],
        "mobile": trace2[second],
        "labels": [labels1[index] for index in first],
        "n_residues_1": len(trace1),
        "n_residues_2": len(trace2),
    }


def _compare_ca_traces(structure1: Dict[str, Any], structure2: Dict[str, Any], pair_by: str) -> Dict[str, Any]:
    """Kabsch RMSD and TM-score of the paired C-alpha atoms of two structures (``calculate_structural_similarity``)."""
    pairs = match_ca_atoms(structure1, structure2, pair_by)
    reference, mobile = pairs["reference"], pairs["mobile"]
    scored = calculate_tm_score(reference, mobile, n_reference=pairs["n_residues_1"])
    return {
        "rmsd": compute_rmsd_kabsch(reference, mobile),
        "tm_score": scored["tm_score"],
        "d0": scored["d0"],
        "n_aligned": len(reference),
        "n_residues_1": pairs["n_residues_1"],
        "n_residues_2": pairs["n_residues_2"],
    }


def calculate_structural_similarity(
    structure1: Dict[str, Any], structure2: Dict[str, Any], *, pair_by: str = "residue_id"
) -> Dict[str, Any]:
    """Compare two structures by the Kabsch-superposed RMSD and the TM-score of their C-alpha atoms.

    The C-alpha atoms (ATOM records) of the two structures are paired; the pairs are superposed and scored:

    - ``rmsd``: the root-mean-square deviation of all pairs after the optimal rigid superposition (Kabsch algorithm
      with a reflection check, :func:`metainformant.protein.structure.general.compute_rmsd_kabsch`), in angstrom. It
      does not depend on how the structures are placed in space.
    - ``tm_score``: the TM-score of Zhang & Skolnick (2004), the maximum over superpositions of
      ``(1 / L) * sum 1 / (1 + (d_i / d0)^2)`` with ``d0 = 1.24 * (L - 15)^(1/3) - 1.8`` (0.5 for ``L <= 21``) and
      ``L`` the number of C-alpha atoms of ``structure1``. It lies in (0, 1] and is 1 for identical structures; its
      d0 scales with the length, so the score of unrelated proteins (about 0.17) does not depend on it. Residues of
      ``structure1`` without a partner score 0.

    No structural alignment is searched for: the pairing must be given by the data.

    Args:
        structure1: Reference structure dictionary from PDB loading (``atoms``; or only ``coordinates``, read as one
            point per residue)
        structure2: Structure dictionary to compare
        pair_by: ``"residue_id"`` pairs residues with the same chain, residue number and insertion code (residues
            present in only one structure are left out of the superposition); ``"order"`` pairs the i-th C-alpha atom
            of one with the i-th of the other, which needs equally many

    Returns:
        Dictionary with ``rmsd`` (angstrom), ``tm_score``, ``d0`` (angstrom), ``n_aligned`` (pairs), ``n_residues_1``
        and ``n_residues_2`` (C-alpha atoms of each structure)

    Raises:
        ValueError: If a structure has no coordinates or C-alpha atoms, fewer than 3 residues can be paired, no
            residue is shared (``pair_by="residue_id"``), the lengths differ (``pair_by="order"``) or ``pair_by`` is
            unknown

    Reference:
        Zhang Y, Skolnick J. Scoring function for automated assessment of protein structure template quality.
        Proteins 2004;57:702-710.

    Example:
        >>> def chain(shift):  # a straight chain of 30 C-alpha atoms, moved along x by ``shift``
        ...     return {"atoms": [{"name": "CA", "chain_id": "A", "res_seq": i, "x": 3.8 * i + shift, "y": 0, "z": 0}
        ...                       for i in range(1, 31)]}
        >>> result = calculate_structural_similarity(chain(0.0), chain(10.0))
        >>> round(result["tm_score"], 6), round(result["rmsd"], 6)
        (1.0, 0.0)
    """
    return _compare_ca_traces(structure1, structure2, pair_by)


# Van der Waals radii (A, Bondi) by element, shared with the SASA calculation; any other element gets the carbon radius.
_VDW_RADII = VAN_DER_WAALS_RADII
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


def calculate_structural_alignment_quality(
    structure1: Dict[str, Any], structure2: Dict[str, Any], *, pair_by: str = "residue_id"
) -> Dict[str, Any]:
    """Assess how well two structures with corresponding residues superpose, by the TM-score and the RMSD.

    This is :func:`calculate_structural_similarity` under the keys of a quality report. ``quality_score`` is the
    TM-score (0 to 1, higher is better; 1 for identical structures). ``rmsd`` is the Kabsch RMSD of all pairs.

    Args:
        structure1: Reference structure dictionary
        structure2: Structure dictionary to assess
        pair_by: ``"residue_id"`` or ``"order"``, as in :func:`calculate_structural_similarity`

    Returns:
        Dictionary with ``quality_score`` (= ``tm_score``), ``tm_score``, ``rmsd`` (angstrom), ``d0`` (angstrom),
        ``aligned_atoms`` (C-alpha atoms paired), ``structure1_atoms`` and ``structure2_atoms`` (C-alpha atoms of each)

    Raises:
        ValueError: As :func:`calculate_structural_similarity`

    Example:
        >>> def chain(shift):  # a straight chain of 30 C-alpha atoms, moved along x by ``shift``
        ...     return {"atoms": [{"name": "CA", "chain_id": "A", "res_seq": i, "x": 3.8 * i + shift, "y": 0, "z": 0}
        ...                       for i in range(1, 31)]}
        >>> quality = calculate_structural_alignment_quality(chain(0.0), chain(-4.0))
        >>> round(quality["quality_score"], 6), quality["aligned_atoms"]
        (1.0, 30)
    """
    comparison = _compare_ca_traces(structure1, structure2, pair_by)
    return {
        "quality_score": comparison["tm_score"],
        "tm_score": comparison["tm_score"],
        "rmsd": comparison["rmsd"],
        "d0": comparison["d0"],
        "aligned_atoms": comparison["n_aligned"],
        "structure1_atoms": comparison["n_residues_1"],
        "structure2_atoms": comparison["n_residues_2"],
    }
