"""Tests: solvent accessible surface areas are computed by the Shrake-Rupley algorithm.

Defect: ``calculate_solvent_accessible_surface_area`` and ``calculate_surface_area`` returned the area of the convex
hull of the atom centres (``calculate_surface_area`` scaled by ``1 + probe_radius / 10``, and for fewer than four
points the area of a sphere through the farthest point). The hull knows neither the atomic radii nor the probe nor any
concavity, so a single carbon atom gave 24.6 A^2 where a water probe rolling on it traces 120.8 A^2. Both now
implement Shrake & Rupley (1973) with Bondi van der Waals radii and a 1.4 A probe.

Oracles are analytic: an isolated atom of van der Waals radius r exposes ``4 pi (r + probe)^2``; for two overlapping
spheres of expanded radii R1 and R2 at distance d, sphere 1 loses a spherical cap, so it exposes
``2 pi R1 (R1 + x1)`` with ``x1 = (d^2 + R1^2 - R2^2) / (2 d)`` (and likewise sphere 2). A buried atom exposes nothing,
and the result is cross-checked against Biopython's independent implementation when Biopython is installed.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List

import numpy as np
import pytest

from metainformant.protein.structure import analysis, general
from metainformant.protein.structure.analysis import calculate_surface_area
from metainformant.protein.structure.general import calculate_solvent_accessible_surface_area

PROBE = 1.4


def sphere_area(radius: float, probe: float = PROBE) -> float:
    return 4.0 * math.pi * (radius + probe) ** 2


def exposed_cap_area(big_r: float, other_r: float, distance: float) -> float:
    """Area of the sphere of radius ``big_r`` left outside a sphere of radius ``other_r`` ``distance`` away."""
    x = (distance**2 + big_r**2 - other_r**2) / (2.0 * distance)
    return 2.0 * math.pi * big_r * (big_r + x)


def icosahedron(radius: float) -> np.ndarray:
    phi = (1.0 + 5.0**0.5) / 2.0
    vertices = []
    for a in (-1.0, 1.0):
        for b in (-phi, phi):
            vertices += [(0.0, a, b), (a, b, 0.0), (b, 0.0, a)]
    unit = np.array(vertices)
    return radius * unit / np.linalg.norm(unit, axis=1, keepdims=True)


# ---------------------------------------------------------------------------
# Analytic oracles
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("radius", [1.2, 1.55, 1.7, 2.75])
@pytest.mark.parametrize("probe", [0.0, 1.4, 2.0])
@pytest.mark.parametrize("n_points", [1, 50, 960])
def test_an_isolated_sphere_exposes_exactly_the_area_of_its_probe_sphere(radius: float, probe: float, n_points: int):
    area = calculate_solvent_accessible_surface_area(
        np.zeros((1, 3)), probe_radius=probe, radii=radius, n_points=n_points
    )

    assert area == pytest.approx(sphere_area(radius, probe), rel=1e-12)


def test_the_example_of_the_old_hull_code_a_single_carbon_atom() -> None:
    area = calculate_surface_area(np.zeros((1, 3)))

    assert area == pytest.approx(4 * math.pi * 3.1**2)  # 120.76 A^2; the convex-hull code gave 24.6
    assert area > 100.0


def test_atoms_that_do_not_touch_add_up() -> None:
    centres = np.array([[0.0, 0.0, 0.0], [30.0, 0.0, 0.0], [0.0, 40.0, 0.0]])
    radii = [1.7, 1.55, 1.52]

    areas = general.calculate_atomic_sasa(centres, radii=radii)

    assert areas == pytest.approx([sphere_area(r) for r in radii], rel=1e-12)
    assert calculate_solvent_accessible_surface_area(centres, radii=radii) == pytest.approx(
        sum(sphere_area(r) for r in radii), rel=1e-12
    )


@pytest.mark.parametrize("distance", [1.0, 2.5, 3.8, 5.5])
def test_two_overlapping_spheres_match_the_spherical_cap_formula_along_the_pole(distance: float) -> None:
    """The axis through the centres is the pole of the point set, where the golden spiral is most regular."""
    r1, r2 = 1.7, 1.55
    big1, big2 = r1 + PROBE, r2 + PROBE
    centres = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, distance]])

    areas = general.calculate_atomic_sasa(centres, radii=[r1, r2], n_points=20000)

    assert areas[0] == pytest.approx(exposed_cap_area(big1, big2, distance), rel=2e-3)
    assert areas[1] == pytest.approx(exposed_cap_area(big2, big1, distance), rel=2e-3)
    assert areas.sum() == pytest.approx(
        exposed_cap_area(big1, big2, distance) + exposed_cap_area(big2, big1, distance), rel=2e-3
    )


def test_two_overlapping_spheres_match_the_cap_formula_along_an_oblique_axis() -> None:
    r1, r2, distance = 1.7, 1.52, 2.9
    big1, big2 = r1 + PROBE, r2 + PROBE
    direction = np.array([1.0, 2.0, -0.5]) / np.linalg.norm([1.0, 2.0, -0.5])
    centres = np.array([[3.0, -2.0, 7.0], np.array([3.0, -2.0, 7.0]) + distance * direction])

    areas = general.calculate_atomic_sasa(centres, radii=[r1, r2], n_points=20000)

    assert areas[0] == pytest.approx(exposed_cap_area(big1, big2, distance), rel=5e-3)
    assert areas[1] == pytest.approx(exposed_cap_area(big2, big1, distance), rel=5e-3)


@pytest.mark.parametrize("n_points, tolerance", [(100, 4e-2), (960, 8e-3), (5000, 3e-3)])
def test_the_error_of_the_estimate_shrinks_with_the_number_of_points(n_points: int, tolerance: float) -> None:
    r1, r2, distance = 1.7, 1.55, 2.5
    exact = exposed_cap_area(r1 + PROBE, r2 + PROBE, distance) + exposed_cap_area(r2 + PROBE, r1 + PROBE, distance)
    direction = np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)
    centres = np.array([np.zeros(3), distance * direction])

    area = calculate_solvent_accessible_surface_area(centres, radii=[r1, r2], n_points=n_points)

    assert abs(area - exact) / exact < tolerance


def test_spheres_just_out_of_reach_are_untouched_and_just_in_reach_lose_a_sliver() -> None:
    r, big = 1.7, 1.7 + PROBE
    apart = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2 * big + 0.01]])
    touching = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2 * big - 0.2]])

    assert general.calculate_atomic_sasa(apart, radii=r, n_points=4000) == pytest.approx(
        [sphere_area(r)] * 2, rel=1e-12
    )
    areas = general.calculate_atomic_sasa(touching, radii=r, n_points=4000)
    assert (areas < sphere_area(r)).all()
    assert areas[0] == pytest.approx(exposed_cap_area(big, big, 2 * big - 0.2), rel=5e-3)


def test_a_sphere_inside_another_is_buried_and_leaves_the_outer_one_whole() -> None:
    centres = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    areas = general.calculate_atomic_sasa(centres, radii=[3.6, 1.6], probe_radius=PROBE)  # expanded radii 5.0 and 3.0

    assert areas[1] == 0.0  # 0.5 + 3.0 < 5.0: the small sphere lies inside the large one
    assert areas[0] == pytest.approx(sphere_area(3.6), rel=1e-12)


def test_an_atom_in_a_shell_of_neighbours_is_buried() -> None:
    centres = np.vstack([np.zeros((1, 3)), icosahedron(1.5)])  # twelve neighbours 1.5 A from the middle one

    areas = general.calculate_atomic_sasa(centres, radii=1.7)

    assert areas[0] == 0.0
    assert (areas[1:] > 0.0).all()  # the neighbours still show their outer faces


def test_the_area_does_not_depend_on_position_or_on_the_order_of_the_atoms() -> None:
    rng = np.random.default_rng(8)
    centres = rng.uniform(0, 10, size=(30, 3))
    radii = rng.choice([1.55, 1.7, 1.8], size=30)
    reference = general.calculate_atomic_sasa(centres, radii=radii)

    order = rng.permutation(30)
    shuffled = general.calculate_atomic_sasa(centres[order], radii=radii[order])
    shifted = general.calculate_atomic_sasa(centres + np.array([100.0, -50.0, 25.0]), radii=radii)

    assert shuffled == pytest.approx(reference[order], abs=1e-9)
    assert shifted == pytest.approx(reference, abs=0.5)  # at most a few of the 960 points per atom sit on a boundary


def test_the_default_radius_is_that_of_carbon_and_one_value_can_stand_for_all() -> None:
    centres = np.random.default_rng(1).uniform(0, 9, size=(25, 3))

    default = calculate_solvent_accessible_surface_area(centres)
    explicit = calculate_solvent_accessible_surface_area(centres, radii=general.VAN_DER_WAALS_RADII["C"])
    per_atom = calculate_solvent_accessible_surface_area(centres, radii=[1.7] * 25)

    assert default == explicit == per_atom


def test_calculate_surface_area_is_the_same_calculation() -> None:
    centres = np.random.default_rng(3).uniform(0, 12, size=(35, 3))

    assert calculate_surface_area(centres, probe_radius=1.0, n_points=200) == calculate_solvent_accessible_surface_area(
        centres, probe_radius=1.0, n_points=200
    )


# ---------------------------------------------------------------------------
# Cross-check against an independent implementation
# ---------------------------------------------------------------------------


def test_the_areas_agree_with_biopythons_shrake_rupley() -> None:
    pytest.importorskip("Bio")
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.Residue import Residue
    from Bio.PDB.SASA import ShrakeRupley
    from Bio.PDB.Structure import Structure

    rng = np.random.default_rng(21)
    n_atoms = 120
    centres = rng.uniform(0, 14, size=(n_atoms, 3))
    elements = rng.choice(["C", "N", "O", "S"], size=n_atoms, p=[0.55, 0.15, 0.25, 0.05])

    structure = Structure("probe")
    model = Model(0)
    chain = Chain("A")
    for serial, (xyz, element) in enumerate(zip(centres, elements), start=1):
        residue = Residue((" ", serial, " "), "ALA", " ")
        residue.add(Atom(f"{element}1", xyz.astype(float), 20.0, 1.0, " ", f" {element}1 ", serial, element=element))
        chain.add(residue)
    model.add(chain)
    structure.add(model)
    ShrakeRupley(probe_radius=1.4, n_points=960).compute(structure, level="A")
    reference = np.array([atom.sasa for atom in structure.get_atoms()])

    radii = [general.van_der_waals_radius(str(element)) for element in elements]
    mine = general.calculate_atomic_sasa(centres, radii=radii, probe_radius=1.4, n_points=960)

    assert reference.sum() > 0
    assert mine == pytest.approx(reference, abs=0.5)  # a few of 960 points may sit exactly on a boundary
    assert mine.sum() == pytest.approx(reference.sum(), rel=3e-3)


# ---------------------------------------------------------------------------
# Points, radii and inputs
# ---------------------------------------------------------------------------


def test_the_golden_spiral_points_are_unit_vectors_spread_over_the_whole_sphere() -> None:
    points = general.golden_spiral_points(960)

    assert points.shape == (960, 3)
    assert np.linalg.norm(points, axis=1) == pytest.approx(1.0, abs=1e-12)
    assert np.linalg.norm(points.mean(axis=0)) < 2e-3  # centred on the origin
    assert np.all(np.diff(points[:, 2]) < 0)  # equal-area bands in z, from the north pole down
    assert points[0, 2] == pytest.approx(1 - 1 / 960)
    assert general.golden_spiral_points(1)[0] == pytest.approx([1.0, 0.0, 0.0])  # a single point sits on the equator
    assert np.array_equal(general.golden_spiral_points(960), points)  # deterministic


def test_the_radii_are_the_bondi_values_by_element() -> None:
    assert [general.van_der_waals_radius(e) for e in ("H", "C", "N", "O", "S", "P")] == [
        1.20,
        1.70,
        1.55,
        1.52,
        1.80,
        1.80,
    ]
    assert general.van_der_waals_radius("zn") == general.van_der_waals_radius(" ZN ") == 1.39
    with pytest.raises(ValueError, match="no van der Waals radius for element 'Xx'"):
        general.van_der_waals_radius("Xx")
    with pytest.raises(ValueError, match="element ''"):
        general.van_der_waals_radius("")


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"probe_radius": -0.1}, "probe_radius"),
        ({"probe_radius": float("nan")}, "probe_radius"),
        ({"n_points": 0}, "n_points"),
        ({"n_points": 12.5}, "n_points"),
        ({"radii": [1.7, 1.7]}, "one value per atom"),
        ({"radii": -1.0}, "positive"),
        ({"radii": 0.0}, "positive"),
    ],
)
def test_invalid_parameters_are_rejected(kwargs: Dict[str, Any], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        general.calculate_atomic_sasa(np.zeros((3, 3)) + np.arange(3)[:, None] * 5.0, **kwargs)


def test_invalid_coordinates_are_rejected_and_no_atoms_have_no_area() -> None:
    with pytest.raises(ValueError, match="shape"):
        general.calculate_atomic_sasa(np.zeros((4, 2)))
    with pytest.raises(ValueError, match="finite"):
        general.calculate_atomic_sasa(np.array([[0.0, 0.0, np.nan]]))
    assert calculate_solvent_accessible_surface_area(np.zeros((0, 3))) == 0.0
    assert general.calculate_atomic_sasa(np.zeros((0, 3))).shape == (0,)


# ---------------------------------------------------------------------------
# Structure dictionaries
# ---------------------------------------------------------------------------


def atom(name: str, xyz, element: str = "", res_name: str = "ALA", res_seq: int = 1, record: str = "ATOM", **extra):
    return {
        "record_name": record,
        "name": name,
        "res_name": res_name,
        "chain_id": "A",
        "res_seq": res_seq,
        "x": float(xyz[0]),
        "y": float(xyz[1]),
        "z": float(xyz[2]),
        "element": element,
        **extra,
    }


def test_structure_sasa_uses_the_radius_of_each_elements_and_sums_per_residue() -> None:
    atoms: List[Dict[str, Any]] = [
        atom("N", (0.0, 0.0, 0.0), "N", res_seq=1),
        atom("CA", (40.0, 0.0, 0.0), "C", res_seq=1),
        atom("O", (0.0, 40.0, 0.0), "O", res_seq=2),
        atom("SG", (0.0, 0.0, 40.0), "S", res_name="CYS", res_seq=3),
    ]

    result = analysis.calculate_structure_sasa({"atoms": atoms})

    expected = [sphere_area(1.55), sphere_area(1.70), sphere_area(1.52), sphere_area(1.80)]
    assert result["per_atom"] == pytest.approx(expected, rel=1e-12)
    assert result["total"] == pytest.approx(sum(expected), rel=1e-12)
    assert result["atom_indices"] == [0, 1, 2, 3]
    assert [r["res_seq"] for r in result["per_residue"]] == [1, 2, 3]
    assert result["per_residue"][0]["sasa"] == pytest.approx(expected[0] + expected[1], rel=1e-12)
    assert result["per_residue"][0]["n_atoms"] == 2
    assert result["per_residue"][2]["res_name"] == "CYS"
    assert (result["n_atoms"], result["probe_radius"], result["n_points"]) == (4, 1.4, 960)


def test_structure_sasa_leaves_out_hetero_atoms_unless_asked_and_counts_their_shadow() -> None:
    protein = [atom("C", (0.0, 0.0, 0.0), "C")]
    water = atom("O", (0.0, 0.0, 3.0), "O", res_name="HOH", res_seq=500, record="HETATM")

    alone = analysis.calculate_structure_sasa({"atoms": protein + [water]})
    together = analysis.calculate_structure_sasa({"atoms": protein + [water]}, include_hetero=True)

    assert alone["total"] == pytest.approx(sphere_area(1.70), rel=1e-12)  # the protein surface alone
    assert alone["n_atoms"] == 1
    assert together["n_atoms"] == 2
    assert together["per_atom"][0] < alone["per_atom"][0]  # the water covers part of the carbon


def test_structure_sasa_takes_the_first_alternate_location_of_each_atom() -> None:
    atoms = [
        atom("C", (0.0, 0.0, 0.0), "C", alt_loc="A"),
        atom("C", (0.2, 0.1, 0.0), "C", alt_loc="B"),  # the same atom in another position: not a second atom
    ]

    result = analysis.calculate_structure_sasa({"atoms": atoms})

    assert result["n_atoms"] == 1
    assert result["total"] == pytest.approx(sphere_area(1.70), rel=1e-12)


def test_structure_sasa_reads_the_element_from_the_atom_name_only_for_protein_atoms() -> None:
    protein_ca = atom("CA", (0.0, 0.0, 0.0), element="")  # a protein C-alpha: carbon
    assert analysis.calculate_structure_sasa({"atoms": [protein_ca]})["total"] == pytest.approx(
        sphere_area(1.70), rel=1e-12
    )

    ion_without_element = atom("ZN", (0.0, 0.0, 0.0), element="", res_name="ZN", record="HETATM")
    with pytest.raises(ValueError, match="no element"):
        analysis.calculate_structure_sasa({"atoms": [ion_without_element]}, include_hetero=True)


def test_structure_sasa_refuses_an_element_without_a_radius_instead_of_guessing() -> None:
    with pytest.raises(ValueError, match="no van der Waals radius for element 'U'"):
        analysis.calculate_structure_sasa({"atoms": [atom("U1", (0.0, 0.0, 0.0), "U")]})


def test_structure_sasa_of_nothing_and_of_bare_coordinates() -> None:
    empty = analysis.calculate_structure_sasa({"atoms": []})
    assert empty["total"] == 0.0 and empty["per_residue"] == [] and empty["n_atoms"] == 0

    with pytest.raises(ValueError, match="no 'atoms' records"):
        analysis.calculate_structure_sasa({"coordinates": np.zeros((3, 3))})
    with pytest.raises(ValueError, match="coordinates"):
        analysis.calculate_structure_sasa({"atoms": [{"name": "C", "element": "C", "res_seq": 1}]})
