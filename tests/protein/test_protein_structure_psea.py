"""Tests: secondary structure elements come from the P-SEA criteria on the C-alpha trace.

Defect: ``identify_secondary_structure_elements`` reported a "helix" with an invented confidence of 0.8 at every
position where CA(i)-CA(i+3) lay between 5 and 7 A, so an ideal helix of 12 residues gave eight overlapping elements
of five atoms each and no strand was ever found. It now implements P-SEA (Labesse et al. 1997): the distances
CA(i)-CA(i+2), CA(i)-CA(i+3) and CA(i)-CA(i+4) and the virtual bond angle and torsion are compared with the windows
of an alpha-helix and a beta-strand.

Oracles: the traces are built from standard geometry, independently of the module under test -- the ideal
alpha-helix (radius 2.3 A, rise 1.5 A, 100 degrees per residue), a chain grown from internal coordinates (3.8 A bonds,
virtual bond angle and torsion), and the pleated extended strand (3.3 A advance, +-0.95 A zigzag). The builders are
checked against distances, angles and torsions measured by separate code below, and the windows of P-SEA are checked
against the geometry of the ideal elements instead of being copied into the expectations.
"""

from __future__ import annotations

import math
from typing import List, Sequence

import numpy as np
import pytest

from metainformant.protein.structure import general
from metainformant.protein.structure.general import identify_secondary_structure_elements

# ---------------------------------------------------------------------------
# Geometry builders and independent measurements
# ---------------------------------------------------------------------------


def ideal_helix(n: int, handedness: int = +1) -> np.ndarray:
    """C-alpha trace of an ideal alpha-helix: radius 2.3 A, rise 1.5 A per residue, 100 degrees per residue."""
    k = np.arange(n)
    angle = np.radians(100.0) * k * handedness
    return np.column_stack([2.3 * np.cos(angle), 2.3 * np.sin(angle), 1.5 * k])


def pleated_strand(n: int) -> np.ndarray:
    """C-alpha trace of a pleated extended strand along x: 3.3 A advance and a +-0.95 A zigzag (3.8 A bonds)."""
    return np.array([[3.3 * i, 0.95 * (-1) ** i, 0.0] for i in range(n)])


def grow_chain(angles: Sequence[float], torsions: Sequence[float], bond: float = 3.8) -> np.ndarray:
    """Chain grown from internal coordinates.

    ``angles[i]`` (i >= 2) is the virtual bond angle at residue i-1 and ``torsions[i]`` (i >= 3) the virtual torsion
    of residues i-3..i, both in degrees. The first entries are ignored.
    """
    points = [np.zeros(3), np.array([bond, 0.0, 0.0])]
    theta = math.radians(angles[2])
    points.append(points[1] + bond * np.array([-math.cos(theta), math.sin(theta), 0.0]))
    for i in range(3, len(angles)):
        a, b, c = points[-3], points[-2], points[-1]
        theta, tau = math.radians(angles[i]), math.radians(torsions[i])
        along = (c - b) / np.linalg.norm(c - b)
        normal = np.cross(b - a, along)
        normal /= np.linalg.norm(normal)
        binormal = np.cross(normal, along)
        local = np.array(
            [-bond * math.cos(theta), bond * math.sin(theta) * math.cos(tau), bond * math.sin(theta) * math.sin(tau)]
        )
        points.append(c + local[0] * along + local[1] * binormal + local[2] * normal)
    return np.array(points)


def uniform_chain(n: int, angle: float, torsion: float) -> np.ndarray:
    return grow_chain([angle] * n, [torsion] * n)


def measured_angle(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> float:
    u, v = p0 - p1, p2 - p1
    return math.degrees(math.acos(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))))


def measured_torsion(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
    m1 = np.cross(n1, b1 / np.linalg.norm(b1))
    return -math.degrees(math.atan2(np.dot(m1, n2), np.dot(n1, n2)))  # IUPAC sign: right-handed twist is positive


def descriptors(ca: np.ndarray, i: int) -> dict:
    """The P-SEA descriptors of residue i, measured with the code above."""
    return {
        "d2": float(np.linalg.norm(ca[i] - ca[i + 2])),
        "d3": float(np.linalg.norm(ca[i] - ca[i + 3])),
        "d4": float(np.linalg.norm(ca[i] - ca[i + 4])),
        "angle": measured_angle(ca[i], ca[i + 1], ca[i + 2]),
        "torsion": measured_torsion(ca[i], ca[i + 1], ca[i + 2], ca[i + 3]),
    }


def in_window(value: float, window: Sequence[float], circular: bool = False) -> bool:
    centre, half_width = window
    delta = value - centre
    if circular:
        delta = (delta + 180.0) % 360.0 - 180.0
    return abs(delta) <= half_width


def rotation_from_seed(seed: int) -> np.ndarray:
    """A proper rotation matrix from a seeded QR decomposition."""
    q, r = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def states_of(ca: np.ndarray) -> str:
    return "".join(general.assign_secondary_structure_psea(ca))


# ---------------------------------------------------------------------------
# The builders and the windows
# ---------------------------------------------------------------------------


def test_the_builders_have_the_geometry_they_claim() -> None:
    helix = ideal_helix(12)
    assert np.allclose(np.linalg.norm(np.diff(helix, axis=0), axis=1), 3.8, atol=0.05)
    first = descriptors(helix, 0)
    assert first["d3"] == pytest.approx(5.05, abs=0.01)  # the textbook CA(i)-CA(i+3) of an alpha-helix
    assert first["d4"] == pytest.approx(6.20, abs=0.01)  # and CA(i)-CA(i+4)
    assert first["torsion"] == pytest.approx(50.0, abs=1.0)  # right-handed: positive virtual torsion
    assert descriptors(ideal_helix(12, handedness=-1), 0)["torsion"] == pytest.approx(-50.0, abs=1.0)

    strand = uniform_chain(12, angle=120.0, torsion=-170.0)
    assert np.allclose(np.linalg.norm(np.diff(strand, axis=0), axis=1), 3.8)
    assert measured_angle(strand[3], strand[4], strand[5]) == pytest.approx(120.0)
    assert measured_torsion(strand[2], strand[3], strand[4], strand[5]) == pytest.approx(-170.0)


@pytest.mark.parametrize(
    "ca, windows_name",
    [
        (ideal_helix(12), "PSEA_HELIX_WINDOWS"),
        (uniform_chain(12, angle=120.0, torsion=-170.0), "PSEA_STRAND_WINDOWS"),
        (uniform_chain(12, angle=124.0, torsion=-170.0), "PSEA_STRAND_WINDOWS"),
        (pleated_strand(12), "PSEA_STRAND_WINDOWS"),
    ],
    ids=["ideal-helix", "strand-120", "strand-124", "pleated-strand"],
)
def test_ideal_geometry_lies_inside_the_published_windows(ca: np.ndarray, windows_name: str) -> None:
    """The windows are tied to the physics: every descriptor of an ideal element falls inside its window."""
    windows = getattr(general, windows_name)
    measured = descriptors(ca, 3)
    for name in ("d2", "d3", "d4", "angle"):
        assert in_window(measured[name], windows[name]), (name, measured[name], windows[name])
    assert in_window(measured["torsion"], windows["torsion"], circular=True), measured["torsion"]


# ---------------------------------------------------------------------------
# Ideal elements are found whole
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [5, 6, 9, 12, 30])
def test_an_ideal_alpha_helix_is_one_helix_element_covering_every_residue(n: int) -> None:
    elements = identify_secondary_structure_elements(ideal_helix(n))

    assert elements == [{"type": "helix", "start_residue": 0, "end_residue": n - 1, "length": n}]
    assert states_of(ideal_helix(n)) == "H" * n


@pytest.mark.parametrize(
    "ca",
    [
        uniform_chain(10, angle=120.0, torsion=-170.0),
        uniform_chain(10, angle=124.0, torsion=-170.0),
        uniform_chain(10, angle=120.0, torsion=180.0),
        pleated_strand(10),
    ],
    ids=["strand-120", "strand-124", "planar-zigzag", "pleated-strand"],
)
def test_an_ideal_beta_strand_is_one_strand_element_covering_every_residue(ca: np.ndarray) -> None:
    elements = identify_secondary_structure_elements(ca)

    assert elements == [{"type": "strand", "start_residue": 0, "end_residue": 9, "length": 10}]


def test_elements_shorter_than_one_window_are_not_reported() -> None:
    assert identify_secondary_structure_elements(ideal_helix(4)) == []
    assert identify_secondary_structure_elements(pleated_strand(4)) == []
    assert states_of(ideal_helix(4)) == "CCCC"
    assert identify_secondary_structure_elements(np.zeros((0, 3))) == []
    assert identify_secondary_structure_elements(ideal_helix(1)) == []


def test_a_helix_loop_strand_chain_is_segmented_into_its_three_parts() -> None:
    """12 helical residues, a four-residue loop that fits neither conformation, 10 extended residues."""
    angles = [91.0] * 12 + [100.0, 105.0, 110.0, 100.0] + [120.0] * 10
    torsions = [50.0] * 12 + [-80.0, 70.0, -120.0, 100.0] + [-170.0] * 10
    ca = grow_chain(angles, torsions)
    assert np.allclose(np.linalg.norm(np.diff(ca, axis=0), axis=1), 3.8)

    states = states_of(ca)

    assert states[:12] == "H" * 12  # the helix is found whole
    assert states[12:14] == "CC"  # the loop is coil where neither conformation fits
    assert states[14:] == "E" * 12  # the strand block, with two loop residues in its first window
    assert identify_secondary_structure_elements(ca) == [
        {"type": "helix", "start_residue": 0, "end_residue": 11, "length": 12},
        {"type": "strand", "start_residue": 14, "end_residue": 25, "length": 12},
    ]


def test_a_chain_without_regular_geometry_has_no_elements() -> None:
    loop = grow_chain([105.0] * 20, [95.0, -60.0] * 10)
    assert states_of(loop) == "C" * 20
    assert identify_secondary_structure_elements(loop) == []


def test_the_old_confidence_and_overlapping_windows_are_gone() -> None:
    elements = identify_secondary_structure_elements(ideal_helix(12))

    assert all("confidence" not in element for element in elements)
    assert len(elements) == 1  # the old code returned one 5-atom element per window


# ---------------------------------------------------------------------------
# Chain breaks, invariance and inputs
# ---------------------------------------------------------------------------


def test_a_gap_in_the_chain_separates_elements() -> None:
    first = ideal_helix(10)
    second = ideal_helix(10) + np.array([0.0, 0.0, 30.0])  # 20 A away from the end of the first: not bonded
    trace = np.vstack([first, second])

    assert identify_secondary_structure_elements(trace) == [
        {"type": "helix", "start_residue": 0, "end_residue": 9, "length": 10},
        {"type": "helix", "start_residue": 10, "end_residue": 19, "length": 10},
    ]


def test_assignment_does_not_depend_on_position_or_orientation() -> None:
    angles = [91.0] * 12 + [100.0, 105.0, 110.0, 100.0] + [120.0] * 10
    torsions = [50.0] * 12 + [-80.0, 70.0, -120.0, 100.0] + [-170.0] * 10
    ca = grow_chain(angles, torsions)
    moved = ca @ rotation_from_seed(4).T + np.array([120.0, -45.0, 8.0])

    assert states_of(moved) == states_of(ca)
    assert identify_secondary_structure_elements(moved) == identify_secondary_structure_elements(ca)


def test_the_distance_criterion_does_not_tell_a_mirror_image_helix_apart() -> None:
    """Documented property: d(i,i+2), d(i,i+3) and d(i,i+4) are the same for a left-handed helix."""
    assert states_of(ideal_helix(12, handedness=-1)) == "H" * 12


def test_ca_indices_pick_the_ca_atoms_out_of_all_atoms() -> None:
    helix = ideal_helix(12)
    atoms = np.empty((4 * len(helix), 3))
    atoms[0::4] = helix  # the CA atoms at every fourth position
    rng = np.random.default_rng(2)
    for offset in (1, 2, 3):  # other backbone atoms: they are not part of the C-alpha trace
        atoms[offset::4] = helix + rng.normal(scale=0.7, size=helix.shape)

    elements = identify_secondary_structure_elements(atoms, list(range(0, len(atoms), 4)))

    assert elements == [{"type": "helix", "start_residue": 0, "end_residue": 11, "length": 12}]


def test_a_helix_read_from_its_other_end_is_the_same_helix() -> None:
    reversed_order = ideal_helix(12)[::-1]

    assert identify_secondary_structure_elements(reversed_order) == [
        {"type": "helix", "start_residue": 0, "end_residue": 11, "length": 12}
    ]


@pytest.mark.parametrize(
    "bad",
    [np.zeros((5, 2)), np.zeros(9), np.array([[0.0, 0.0, np.nan]] * 6), np.array([[0.0, np.inf, 0.0]] * 6)],
    ids=["two-columns", "flat", "nan", "inf"],
)
def test_invalid_coordinates_are_rejected(bad: np.ndarray) -> None:
    with pytest.raises(ValueError):
        general.assign_secondary_structure_psea(bad)
    with pytest.raises(ValueError):
        identify_secondary_structure_elements(bad)


def test_invalid_ca_indices_are_rejected() -> None:
    coords = ideal_helix(8)
    with pytest.raises(ValueError, match="ca_indices"):
        identify_secondary_structure_elements(coords, [0, 1, 99])
    with pytest.raises(ValueError, match="ca_indices"):
        identify_secondary_structure_elements(coords, [-1, 0, 1])
    with pytest.raises(ValueError, match="integers"):
        identify_secondary_structure_elements(coords, [0.5, 1.5])
    with pytest.raises(ValueError, match="integers"):
        identify_secondary_structure_elements(coords, [True, False, True])


def test_the_states_use_the_letters_of_the_secondary_structure_module() -> None:
    from metainformant.protein.structure.secondary import calculate_ss_composition, identify_ss_elements

    states: List[str] = general.assign_secondary_structure_psea(ideal_helix(12))

    assert calculate_ss_composition(states)["helix"] == 1.0
    assert identify_ss_elements(states, min_length=5) == [{"type": "H", "start": 0, "end": 11, "length": 12}]
