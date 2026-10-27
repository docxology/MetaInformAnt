"""Regression tests: structural motif search and ligand-pocket detection work on the coordinates they are given.

Defect: ``analyze_structural_motifs`` carried the note "Placeholder for structural motif analysis" and, for any
structure with more than 20 coordinates, returned a ``potential_helix`` (confidence 0.5) for *every* window of ten
coordinates whatever their geometry. ``identify_ligand_binding_sites`` ("Simplified placeholder") returned, for any
structure with 20 or more coordinates, one "site" made of the atoms farther than the median distance from the
centroid, a volume of ten cubic Angstrom per such atom and a fixed confidence of 0.3. Both are now real geometric
searches (C-alpha trace geometry and zinc coordination; grid-based buriedness) and report nothing where nothing is
there.

Oracles: structures are built from analytic geometry -- an ideal right-handed alpha-helix (radius 2.3 A, rise 1.5 A,
100 degrees per residue), a pleated extended strand (3.3 A advance, +-0.95 A zigzag), a dense spherical shell of atoms
around an empty cavity whose accessible volume follows from the shell radius, the atomic radius and the probe
radius, and a zinc ion with ligands placed at chosen distances along tetrahedral directions. The expected results
follow from how the structures were built, and the cavity volume is also recomputed by brute force on an
independent, finer grid.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np
import pytest

from metainformant.protein.structure.analysis import analyze_structural_motifs, identify_ligand_binding_sites
from metainformant.protein.structure.io import parse_pdb_file
from metainformant.protein.structure.pdb import load_pdb_file

# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------


def helix(n: int, handed: int = +1) -> np.ndarray:
    """C-alpha trace of an ideal alpha-helix with ``n`` residues (right-handed for ``handed=+1``)."""
    k = np.arange(n)
    angle = np.radians(100.0) * k * handed
    return np.column_stack([2.3 * np.cos(angle), 2.3 * np.sin(angle), 1.5 * k])


def strand(n: int) -> np.ndarray:
    """C-alpha trace of a pleated extended strand along +x (3.8 A between neighbours)."""
    return np.array([[3.3 * i, 0.95 * (-1) ** i, 0.0] for i in range(n)])


def attach(previous_end: np.ndarray, block: np.ndarray, direction: Sequence[float]) -> np.ndarray:
    """Translate ``block`` so that its first point is 3.8 A from ``previous_end`` along ``direction``."""
    unit = np.asarray(direction, dtype=float)
    unit = unit / np.linalg.norm(unit)
    return block - block[0] + previous_end + 3.8 * unit


def ca_atoms(coords: np.ndarray, chain: str = "A", first_residue: int = 1, **extra: Any) -> List[Dict[str, Any]]:
    """Atom dictionaries (the keys of ``parse_pdb_file``) for a C-alpha trace."""
    return [
        {
            "record_name": "ATOM",
            "serial": i + 1,
            "name": "CA",
            "res_name": "ALA",
            "chain_id": chain,
            "res_seq": first_residue + i,
            "x": float(c[0]),
            "y": float(c[1]),
            "z": float(c[2]),
            "element": "C",
            **extra,
        }
        for i, c in enumerate(coords)
    ]


def structure_of(atoms: List[Dict[str, Any]]) -> Dict[str, Any]:
    return {"atoms": atoms, "coordinates": np.array([[a["x"], a["y"], a["z"]] for a in atoms])}


def strand_helix_strand() -> np.ndarray:
    """Residues 1-6 extended, 7-20 an alpha-helix, 21-26 extended; consecutive C-alpha atoms 3.8 A apart."""
    first = strand(6)
    middle = attach(first[-1], helix(14), [0.3, 0.8, 0.5])
    last = attach(middle[-1], strand(6)[:, [1, 0, 2]], [0.2, -0.7, 0.6])
    return np.vstack([first, middle, last])


EXPECTED_STRAND_HELIX_STRAND = [
    {"type": "extended", "chain_id": "A", "start_residue": 1, "end_residue": 6, "length": 6},
    {"type": "helix", "chain_id": "A", "start_residue": 7, "end_residue": 20, "length": 14},
    {"type": "extended", "chain_id": "A", "start_residue": 21, "end_residue": 26, "length": 6},
]


def fibonacci_sphere(n: int, radius: float, center: Sequence[float] = (0.0, 0.0, 0.0)) -> np.ndarray:
    """``n`` points spread evenly over a sphere."""
    k = np.arange(n) + 0.5
    polar = np.arccos(1 - 2 * k / n)
    azimuth = math.pi * (1 + 5**0.5) * k
    unit = np.column_stack([np.cos(azimuth) * np.sin(polar), np.sin(azimuth) * np.sin(polar), np.cos(polar)])
    return unit * radius + np.asarray(center, dtype=float)


def atoms_at(coords: np.ndarray, residue: int = 1, record: str = "ATOM", res_name: str = "ALA") -> List[Dict[str, Any]]:
    """Carbon atoms at ``coords``, all in one residue."""
    return [
        {
            "record_name": record,
            "serial": i + 1,
            "name": "C",
            "res_name": res_name,
            "chain_id": "A",
            "res_seq": residue,
            "x": float(c[0]),
            "y": float(c[1]),
            "z": float(c[2]),
            "element": "C",
        }
        for i, c in enumerate(coords)
    ]


def pdb_line(
    serial: int, name: str, residue: str, chain: str, number: int, xyz: Sequence[float], record: str = "ATOM  "
) -> str:
    element = name[0]
    return (
        f"{record}{serial:5d} {' ' + name:<4s} {residue:>3s} {chain}{number:4d}    "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}  "
    )


# ---------------------------------------------------------------------------
# Motifs: helices and extended segments
# ---------------------------------------------------------------------------


def test_the_built_chain_has_ideal_bonds_and_helix_geometry() -> None:
    """Guard the builders: 3.8 A bonds, and the textbook C-alpha distances of an alpha-helix."""
    coords = strand_helix_strand()
    bonds = np.linalg.norm(coords[1:] - coords[:-1], axis=1)
    ideal = helix(10)

    assert bonds.min() == pytest.approx(3.8, abs=0.05) and bonds.max() == pytest.approx(3.8, abs=0.05)
    assert np.linalg.norm(ideal[3] - ideal[0]) == pytest.approx(5.05, abs=0.01)  # CA(i)-CA(i+3)
    assert np.linalg.norm(ideal[4] - ideal[0]) == pytest.approx(6.20, abs=0.01)  # CA(i)-CA(i+4)


def test_a_helix_between_two_strands_is_located_exactly() -> None:
    motifs = analyze_structural_motifs(structure_of(ca_atoms(strand_helix_strand())))

    assert motifs == EXPECTED_STRAND_HELIX_STRAND


def test_a_random_coil_has_no_helix_and_none_of_the_old_fabricated_entries() -> None:
    """A 40-residue random walk used to yield a 'potential_helix' for each of its 30 windows of ten coordinates."""
    rng = np.random.default_rng(3)
    steps = rng.normal(size=(39, 3))
    steps = 3.8 * steps / np.linalg.norm(steps, axis=1, keepdims=True)
    coil = np.vstack([np.zeros(3), np.cumsum(steps, axis=0)])

    motifs = analyze_structural_motifs(structure_of(ca_atoms(coil)))

    assert [m for m in motifs if m["type"] == "helix"] == []
    assert all(m["type"] in ("helix", "extended", "zinc_site") for m in motifs)  # no 'potential_helix'
    assert all("confidence" not in m for m in motifs)  # the old fixed confidence of 0.5 is gone
    assert len(motifs) < 3  # at most the odd stretched-out stretch of a coil, not one entry per window


def test_a_cloud_of_points_that_is_not_a_chain_has_no_helix() -> None:
    cloud = np.random.default_rng(1).uniform(0, 30, size=(60, 3))

    assert analyze_structural_motifs(structure_of(ca_atoms(cloud))) == []


def test_a_left_handed_helix_is_not_reported_as_a_helix() -> None:
    mirror = np.vstack([strand(6), attach(strand(6)[-1], helix(14, handed=-1), [0.3, 0.8, 0.5])])

    motifs = analyze_structural_motifs(structure_of(ca_atoms(mirror)))

    assert [m["type"] for m in motifs] == ["extended"]
    assert (motifs[0]["start_residue"], motifs[0]["end_residue"]) == (1, 6)


def test_a_helix_is_found_wherever_it_sits_and_however_it_is_oriented() -> None:
    """The criteria are internal distances and torsions, so rigid motions leave the result unchanged."""
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]) @ np.array(
        [[1.0, 0.0, 0.0], [0.0, 0.6, -0.8], [0.0, 0.8, 0.6]]
    )
    moved = strand_helix_strand() @ rotation.T + np.array([120.0, -40.0, 7.5])

    assert analyze_structural_motifs(structure_of(ca_atoms(moved))) == EXPECTED_STRAND_HELIX_STRAND


def test_helices_shorter_than_six_residues_are_not_reported() -> None:
    five = strand_helix_strand()
    short = np.vstack([strand(6), attach(strand(6)[-1], helix(5), [0.3, 0.8, 0.5])])

    assert [m for m in analyze_structural_motifs(structure_of(ca_atoms(short))) if m["type"] == "helix"] == []
    assert [m for m in analyze_structural_motifs(structure_of(ca_atoms(five))) if m["type"] == "helix"]
    six = np.vstack([strand(6), attach(strand(6)[-1], helix(6), [0.3, 0.8, 0.5])])
    helices = [m for m in analyze_structural_motifs(structure_of(ca_atoms(six))) if m["type"] == "helix"]
    assert [(m["start_residue"], m["end_residue"]) for m in helices] == [(7, 12)]


def test_a_gap_in_the_chain_separates_two_helices_even_with_consecutive_numbering() -> None:
    first = helix(10)
    second = helix(10) + np.array([0.0, 0.0, first[-1][2] + 9.0])  # the next C-alpha is 9 A away: residues are missing
    coords = np.vstack([first, second])

    motifs = analyze_structural_motifs(structure_of(ca_atoms(coords)), motif_types=["helix"])

    assert [(m["start_residue"], m["end_residue"], m["length"]) for m in motifs] == [(1, 10, 10), (11, 20, 10)]


def test_each_chain_is_searched_on_its_own_and_reported_in_order() -> None:
    chain_a = ca_atoms(helix(12), chain="A")
    chain_b = ca_atoms(helix(9) + np.array([40.0, 0.0, 0.0]), chain="B")

    motifs = analyze_structural_motifs(structure_of(chain_a + chain_b))

    assert [(m["chain_id"], m["start_residue"], m["end_residue"]) for m in motifs] == [("A", 1, 12), ("B", 1, 9)]


def test_hetero_atoms_named_ca_are_not_part_of_the_trace() -> None:
    """A calcium ion has the atom name CA too; it must not enter the C-alpha trace."""
    atoms = ca_atoms(strand_helix_strand())
    calcium = {**atoms[0], "record_name": "HETATM", "res_name": "CA", "element": "CA", "serial": 900, "res_seq": 301}
    calcium.update({"x": 3.0, "y": 3.0, "z": 3.0})
    water = {**calcium, "name": "O", "res_name": "HOH", "element": "O", "serial": 901, "res_seq": 302}

    assert analyze_structural_motifs(structure_of(atoms + [calcium, water])) == EXPECTED_STRAND_HELIX_STRAND


def test_alternate_locations_do_not_count_a_residue_twice() -> None:
    atoms = ca_atoms(strand_helix_strand())
    duplicates = [{**a, "alt_loc": "B", "x": a["x"] + 0.4} for a in atoms[8:12]]

    assert analyze_structural_motifs(structure_of(atoms + duplicates)) == EXPECTED_STRAND_HELIX_STRAND


def test_residue_numbers_come_from_the_structure_not_from_positions() -> None:
    atoms = ca_atoms(strand_helix_strand(), first_residue=101)

    motifs = analyze_structural_motifs(structure_of(atoms))

    assert [(m["start_residue"], m["end_residue"]) for m in motifs] == [(101, 106), (107, 120), (121, 126)]


# ---------------------------------------------------------------------------
# Motifs: zinc sites
# ---------------------------------------------------------------------------

TETRAHEDRON = [np.array(v) / math.sqrt(3.0) for v in ([1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1])]


def ligand(
    name: str, residue: str, number: int, direction: np.ndarray, distance: float, serial: int, **extra: Any
) -> Dict[str, Any]:
    position = distance * direction
    return {
        "record_name": "ATOM",
        "serial": serial,
        "name": name,
        "res_name": residue,
        "chain_id": "A",
        "res_seq": number,
        "x": float(position[0]),
        "y": float(position[1]),
        "z": float(position[2]),
        "element": name[0],
        **extra,
    }


def zinc_ion(position: Sequence[float] = (0.0, 0.0, 0.0), number: int = 200, serial: int = 99) -> Dict[str, Any]:
    return {
        "record_name": "HETATM",
        "serial": serial,
        "name": "ZN",
        "res_name": "ZN",
        "chain_id": "A",
        "res_seq": number,
        "x": float(position[0]),
        "y": float(position[1]),
        "z": float(position[2]),
        "element": "ZN",
    }


def test_a_cys2_his2_site_is_reported_with_its_ligands_and_distances() -> None:
    atoms = [
        zinc_ion(),
        ligand("SG", "CYS", 5, TETRAHEDRON[0], 2.30, 1),
        ligand("SG", "CYS", 8, TETRAHEDRON[1], 2.30, 2),
        ligand("NE2", "HIS", 21, TETRAHEDRON[2], 2.05, 3),
        ligand("ND1", "HIS", 25, TETRAHEDRON[3], 2.05, 4),
        ligand("SG", "CYS", 40, np.array([1.0, 0.0, 0.0]), 3.5, 5),  # too far to coordinate
        ligand("NZ", "LYS", 41, np.array([0.0, 1.0, 0.0]), 2.5, 6),  # close, but not a zinc ligand type
        ligand("SG", "CYS", 42, np.array([0.0, 0.0, 1.0]), 2.2, 7, record_name="HETATM"),  # hetero atom: not protein
    ]

    sites = analyze_structural_motifs(structure_of(atoms), motif_types=["zinc_site"])

    assert len(sites) == 1
    site = sites[0]
    assert site["type"] == "zinc_site"
    assert site["coordination_pattern"] == "C2H2"
    assert site["n_ligand_residues"] == 4
    assert site["zinc_finger_like"] is True
    assert site["position"] == [0.0, 0.0, 0.0]
    assert [(lig["res_name"], lig["res_seq"], lig["atom"]) for lig in site["ligands"]] == [
        ("HIS", 21, "NE2"),
        ("HIS", 25, "ND1"),
        ("CYS", 5, "SG"),
        ("CYS", 8, "SG"),
    ]  # nearest first
    assert [lig["distance"] for lig in site["ligands"]] == pytest.approx([2.05, 2.05, 2.30, 2.30])


@pytest.mark.parametrize(
    ("residues", "pattern", "finger_like"),
    [
        (["CYS", "CYS", "CYS", "CYS"], "C4", True),
        (["CYS", "CYS", "CYS", "HIS"], "C3H", True),
        (["CYS", "HIS", "HIS", "HIS"], "CH3", True),
        (["CYS", "CYS", "HIS"], "C2H", False),  # three ligands only
        (["CYS", "CYS", "HIS", "ASP"], "C2HD", False),  # a carboxylate ligand: a catalytic-type site
        (["HIS", "HIS", "GLU", "ASP"], "H2DE", False),
    ],
)
def test_coordination_patterns_and_the_zinc_finger_signature(
    residues: List[str], pattern: str, finger_like: bool
) -> None:
    names = {"CYS": "SG", "HIS": "NE2", "ASP": "OD1", "GLU": "OE1"}
    atoms = [zinc_ion()] + [
        ligand(names[residue], residue, 10 + i, TETRAHEDRON[i], 2.1, i + 1) for i, residue in enumerate(residues)
    ]

    (site,) = analyze_structural_motifs(structure_of(atoms), motif_types=["zinc_site"])

    assert site["coordination_pattern"] == pattern
    assert site["zinc_finger_like"] is finger_like


def test_a_bidentate_carboxylate_counts_as_one_ligand_residue() -> None:
    atoms = [
        zinc_ion(),
        ligand("SG", "CYS", 5, TETRAHEDRON[0], 2.3, 1),
        ligand("NE2", "HIS", 9, TETRAHEDRON[1], 2.05, 2),
        ligand("OD1", "ASP", 12, TETRAHEDRON[2], 2.1, 3),
        ligand("OD2", "ASP", 12, TETRAHEDRON[3], 2.4, 4),  # the same residue binds with both oxygens
    ]

    (site,) = analyze_structural_motifs(structure_of(atoms), motif_types=["zinc_site"])

    assert len(site["ligands"]) == 4  # four ligand atoms ...
    assert site["n_ligand_residues"] == 3  # ... of three residues
    assert site["coordination_pattern"] == "CHD"
    assert site["zinc_finger_like"] is False


def test_each_zinc_ion_gets_its_own_site_and_unbound_ions_are_left_out() -> None:
    far = np.array([50.0, 0.0, 0.0])
    atoms = [zinc_ion(number=200, serial=98), zinc_ion(far, number=201, serial=99), zinc_ion((0.0, 80.0, 0.0), 202, 97)]
    atoms += [ligand("SG", "CYS", 5 + i, TETRAHEDRON[i], 2.3, i + 1) for i in range(4)]
    atoms += [
        {**ligand("NE2", "HIS", 30 + i, TETRAHEDRON[i], 2.05, 10 + i), "x": 50.0 + 2.05 * TETRAHEDRON[i][0]}
        | {"y": 2.05 * TETRAHEDRON[i][1], "z": 2.05 * TETRAHEDRON[i][2]}
        for i in range(4)
    ]

    sites = analyze_structural_motifs(structure_of(atoms), motif_types=["zinc_site"])

    assert [(s["ion_res_seq"], s["coordination_pattern"]) for s in sites] == [(200, "C4"), (201, "H4")]


def test_zinc_sites_follow_the_secondary_structure_elements() -> None:
    atoms = ca_atoms(strand_helix_strand()) + [zinc_ion((100.0, 100.0, 100.0))]
    atoms += [
        {**ligand("SG", "CYS", 700 + i, TETRAHEDRON[i], 2.3, 500 + i), "x": 100.0 + 2.3 * TETRAHEDRON[i][0]}
        | {"y": 100.0 + 2.3 * TETRAHEDRON[i][1], "z": 100.0 + 2.3 * TETRAHEDRON[i][2]}
        for i in range(4)
    ]

    motifs = analyze_structural_motifs(structure_of(atoms))

    assert [m["type"] for m in motifs] == ["extended", "helix", "extended", "zinc_site"]


# ---------------------------------------------------------------------------
# Motifs: arguments, empty input, both PDB loaders
# ---------------------------------------------------------------------------


def test_motif_types_select_what_is_searched_for() -> None:
    structure = structure_of(ca_atoms(strand_helix_strand()))

    helices = analyze_structural_motifs(structure, motif_types=["helix"])
    extended = analyze_structural_motifs(structure, motif_types=("extended",))
    as_text = analyze_structural_motifs(structure, motif_types="helix")

    assert helices == [EXPECTED_STRAND_HELIX_STRAND[1]] == as_text
    assert extended == [EXPECTED_STRAND_HELIX_STRAND[0], EXPECTED_STRAND_HELIX_STRAND[2]]


def test_helix_turn_helix_is_not_implemented_instead_of_guessed() -> None:
    structure = structure_of(ca_atoms(strand_helix_strand()))

    with pytest.raises(NotImplementedError, match="helix_turn_helix is not implemented"):
        analyze_structural_motifs(structure, motif_types=["helix", "helix_turn_helix"])


def test_an_unknown_motif_type_is_a_value_error() -> None:
    with pytest.raises(ValueError, match="Unknown motif type 'beta_barrel'"):
        analyze_structural_motifs(structure_of(ca_atoms(strand_helix_strand())), motif_types=["beta_barrel"])


def test_empty_structures_have_no_motifs_and_bare_coordinates_are_refused() -> None:
    assert analyze_structural_motifs({"atoms": []}) == []
    assert analyze_structural_motifs({}) == []
    assert analyze_structural_motifs({"atoms": [], "coordinates": np.empty((0, 3))}) == []
    with pytest.raises(ValueError, match="no 'atoms' records"):
        analyze_structural_motifs({"coordinates": strand_helix_strand()})


def test_both_pdb_loaders_give_the_same_motifs(tmp_path: Path) -> None:
    coords = strand_helix_strand()
    lines = [pdb_line(i + 1, "CA", "ALA", "A", i + 1, c) for i, c in enumerate(coords)]
    lines.append(pdb_line(900, "ZN", "ZN", "A", 300, (200.0, 0.0, 0.0), record="HETATM"))
    path = tmp_path / "chain.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")

    from_io = analyze_structural_motifs(parse_pdb_file(path))
    from_pdb = analyze_structural_motifs(load_pdb_file(path))

    assert from_io == from_pdb == EXPECTED_STRAND_HELIX_STRAND  # the lone ion has no ligands, so no site


def test_the_module_no_longer_carries_the_placeholder_markers() -> None:
    import metainformant.protein.structure.analysis as analysis

    source = Path(analysis.__file__).read_text().lower()

    assert "placeholder" not in source
    assert "potential_helix" not in source


# ---------------------------------------------------------------------------
# Pockets
# ---------------------------------------------------------------------------

SHELL_RADIUS = 7.0
CARBON_RADIUS = 1.70
WATER_PROBE = 1.40


def shell_structure(
    radius: float = SHELL_RADIUS, n_atoms: int = 600, center: Sequence[float] = (0, 0, 0), residue: int = 1
):
    return structure_of(atoms_at(fibonacci_sphere(n_atoms, radius, center), residue=residue))


def accessible_volume_by_brute_force(
    coords: np.ndarray, center: np.ndarray, reach: float, minimum_distance: float
) -> float:
    """Volume of the points inside a ball around ``center`` that are at least ``minimum_distance`` from every atom.

    Computed on its own grid of 0.25 A, independently of the code under test.
    """
    step = 0.25
    axis = np.arange(-reach, reach + step / 2, step)
    grid = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    grid = grid[np.linalg.norm(grid, axis=1) <= reach] + center
    nearest = np.empty(len(grid))
    for start in range(0, len(grid), 20000):
        chunk = grid[start : start + 20000]
        nearest[start : start + 20000] = np.min(np.linalg.norm(chunk[:, None, :] - coords[None, :, :], axis=2), axis=1)
    return float(np.count_nonzero(nearest >= minimum_distance)) * step**3


def test_a_closed_cavity_is_found_at_its_centre_with_the_expected_volume() -> None:
    structure = shell_structure()

    pockets = identify_ligand_binding_sites(structure, ray_length=16.0)

    assert len(pockets) == 1
    pocket = pockets[0]
    assert pocket["rank"] == 1
    assert pocket["center"] == pytest.approx([0.0, 0.0, 0.0], abs=0.3)
    # The accessible ball reaches to the radius at which a probe still clears the shell atoms.
    analytic = 4 / 3 * math.pi * (SHELL_RADIUS - CARBON_RADIUS - WATER_PROBE) ** 3  # 248.5 A^3
    assert pocket["volume"] == pytest.approx(analytic, rel=0.08)
    brute_force = accessible_volume_by_brute_force(
        structure["coordinates"], np.zeros(3), reach=6.5, minimum_distance=CARBON_RADIUS + WATER_PROBE
    )
    assert pocket["volume"] == pytest.approx(brute_force, rel=0.06)
    assert pocket["n_grid_points"] == round(pocket["volume"])  # 1 A cells
    assert pocket["mean_buriedness"] == pytest.approx(7.0)  # every direction is blocked on both sides


def test_the_pocket_is_lined_by_the_atoms_around_it() -> None:
    structure = shell_structure()

    (pocket,) = identify_ligand_binding_sites(structure, ray_length=16.0)

    assert len(pocket["lining_atoms"]) >= 0.95 * len(structure["atoms"])
    assert pocket["lining_atoms"] == sorted(pocket["lining_atoms"])
    assert pocket["lining_residues"] == [{"chain_id": "A", "res_seq": 1, "res_name": "ALA"}]


def test_a_solid_ball_and_a_flat_slab_have_no_pocket() -> None:
    """The old code returned a 'site' for any structure of 20 or more coordinates."""
    grid = np.arange(-8.0, 8.1, 1.5)
    ball = np.array([(x, y, z) for x in grid for y in grid for z in grid if x * x + y * y + z * z <= 64.0])
    slab = np.array([(x, y, 0.0) for x in np.arange(-10, 10.1, 1.5) for y in np.arange(-10, 10.1, 1.5)])

    assert len(ball) > 20 and len(slab) > 20
    assert identify_ligand_binding_sites(structure_of(atoms_at(ball))) == []
    assert identify_ligand_binding_sites(structure_of(atoms_at(slab))) == []


def test_two_cavities_are_ranked_by_volume_and_keep_their_own_lining() -> None:
    small = atoms_at(fibonacci_sphere(600, 7.0, (0, 0, 0)), residue=1)
    large = atoms_at(fibonacci_sphere(900, 9.0, (40, 0, 0)), residue=2)

    pockets = identify_ligand_binding_sites(structure_of(small + large), ray_length=20.0)

    assert [p["rank"] for p in pockets] == [1, 2]
    assert pockets[0]["volume"] > pockets[1]["volume"]
    assert pockets[0]["center"] == pytest.approx([40.0, 0.0, 0.0], abs=0.3)  # the larger cavity is ranked first
    assert pockets[1]["center"] == pytest.approx([0.0, 0.0, 0.0], abs=0.3)
    assert [r["res_seq"] for r in pockets[0]["lining_residues"]] == [2]
    assert [r["res_seq"] for r in pockets[1]["lining_residues"]] == [1]
    assert pockets[0]["volume"] == pytest.approx(4 / 3 * math.pi * (9.0 - 3.1) ** 3, rel=0.08)
    assert pockets[1]["volume"] == pytest.approx(4 / 3 * math.pi * (7.0 - 3.1) ** 3, rel=0.08)


def test_the_probe_radius_sets_how_much_of_the_cavity_is_accessible() -> None:
    structure = shell_structure()
    for probe in (0.0, 1.0, 1.4, 2.0):
        (pocket,) = identify_ligand_binding_sites(structure, probe_radius=probe, ray_length=16.0)

        expected = 4 / 3 * math.pi * (SHELL_RADIUS - CARBON_RADIUS - probe) ** 3
        assert pocket["volume"] == pytest.approx(expected, rel=0.12), probe


def test_opening_the_wall_lowers_the_buriedness_and_a_wide_opening_removes_the_pocket() -> None:
    """Rays that leave through the opening are not buried: buriedness falls as the cap removed from the shell grows.

    Removing wall atoms also adds accessible space, so the volume is not monotonic for a small hole; what must
    hold is that buriedness falls, the pocket retreats towards the closed side and a hemisphere-plus opening has none.
    """
    coords = fibonacci_sphere(900, SHELL_RADIUS)
    polar = np.degrees(np.arccos(coords[:, 2] / SHELL_RADIUS))
    pocket_by_opening = {}
    for opening in (0, 25, 60, 80, 100):  # degrees of polar angle around +z that are left without wall
        pockets = identify_ligand_binding_sites(structure_of(atoms_at(coords[polar >= opening])), ray_length=16.0)
        pocket_by_opening[opening] = pockets[0] if pockets else None

    assert pocket_by_opening[100] is None  # more than a hemisphere open: nothing is buried enough
    buriedness = [pocket_by_opening[o]["mean_buriedness"] for o in (0, 25, 60, 80)]
    assert buriedness[0] == pytest.approx(7.0)
    assert all(earlier > later for earlier, later in zip(buriedness, buriedness[1:]))
    closed_volume = pocket_by_opening[0]["volume"]
    assert pocket_by_opening[60]["volume"] < 0.8 * closed_volume
    assert pocket_by_opening[80]["volume"] < 0.2 * closed_volume
    assert pocket_by_opening[60]["center"][2] < -0.5  # the buried part is on the side away from the opening (+z)


def test_waters_ions_and_ligands_do_not_fill_the_pocket() -> None:
    structure = shell_structure()
    bound = [
        *atoms_at(np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]]), residue=500, record="HETATM", res_name="LIG"),
        *atoms_at(np.array([[0.0, 2.0, 0.0]]), residue=501, record="HETATM", res_name="HOH"),
    ]
    with_ligand = structure_of(structure["atoms"] + bound)

    assert identify_ligand_binding_sites(with_ligand, ray_length=16.0) == identify_ligand_binding_sites(
        structure, ray_length=16.0
    )


def test_atoms_that_are_all_hetero_atoms_leave_nothing_to_analyse() -> None:
    assert identify_ligand_binding_sites(structure_of(atoms_at(fibonacci_sphere(600, 7.0), record="HETATM"))) == []


def test_bare_coordinates_are_analysed_with_carbon_radii() -> None:
    coordinates = fibonacci_sphere(600, SHELL_RADIUS)

    (pocket,) = identify_ligand_binding_sites({"coordinates": coordinates}, ray_length=16.0)

    assert pocket["center"] == pytest.approx([0.0, 0.0, 0.0], abs=0.3)
    assert pocket["lining_residues"] == []
    assert max(pocket["lining_atoms"]) < len(coordinates)  # indices into the coordinates
    (with_atoms,) = identify_ligand_binding_sites(shell_structure(), ray_length=16.0)
    assert pocket["volume"] == with_atoms["volume"]


def test_the_volume_threshold_and_buriedness_threshold_filter_pockets() -> None:
    structure = shell_structure()

    assert identify_ligand_binding_sites(structure, ray_length=16.0, min_volume=300.0) == []
    assert len(identify_ligand_binding_sites(structure, ray_length=16.0, min_volume=240.0)) == 1
    # Rays of 8 A cannot cross a 14 A cavity: far from the centre fewer than seven directions are buried.
    short = identify_ligand_binding_sites(structure, ray_length=8.0, min_buriedness=7, min_volume=0.0)
    full = identify_ligand_binding_sites(structure, ray_length=8.0, min_buriedness=5, min_volume=0.0)
    assert (
        short[0]["volume"] < full[0]["volume"] < identify_ligand_binding_sites(structure, ray_length=16.0)[0]["volume"]
    )


def test_the_result_is_deterministic_and_plain_data() -> None:
    structure = shell_structure()

    first = identify_ligand_binding_sites(structure, ray_length=16.0)
    second = identify_ligand_binding_sites(structure, ray_length=16.0)

    assert first == second
    pocket = first[0]
    assert set(pocket) == {
        "rank",
        "center",
        "volume",
        "n_grid_points",
        "mean_buriedness",
        "lining_atoms",
        "lining_residues",
    }
    assert all(isinstance(value, float) for value in pocket["center"])
    assert all(isinstance(index, int) for index in pocket["lining_atoms"])


def test_empty_input_and_invalid_parameters() -> None:
    assert identify_ligand_binding_sites({"atoms": []}) == []
    assert identify_ligand_binding_sites({}) == []
    structure = shell_structure()
    for bad in ({"grid_spacing": 0.0}, {"ray_length": -1.0}, {"probe_radius": -0.1}, {"min_volume": -1.0}):
        with pytest.raises(ValueError, match="must"):
            identify_ligand_binding_sites(structure, **bad)
    for bad_buriedness in (0, 8):
        with pytest.raises(ValueError, match="min_buriedness"):
            identify_ligand_binding_sites(structure, min_buriedness=bad_buriedness)


def test_a_grid_that_would_not_fit_in_memory_is_refused_before_it_is_built() -> None:
    with pytest.raises(ValueError, match="increase grid_spacing"):
        identify_ligand_binding_sites(shell_structure(), grid_spacing=0.01)
