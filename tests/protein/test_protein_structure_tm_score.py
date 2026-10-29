"""Tests: structural similarity is the Kabsch-superposed RMSD and the TM-score, not an arbitrary rescaling.

Defect: ``calculate_structural_similarity`` took the RMSD of the coordinates *without superposing them* and turned it
into ``1 / (1 + rmsd / 10)``; ``calculate_structural_alignment_quality`` did the same with ``max(0, 1 - rmsd / 20)``.
Two copies of one structure in different poses were therefore "dissimilar", and the scales (10 and 20 A) had no
meaning. Both now report the RMSD after the optimal rigid superposition (Kabsch, SVD with a reflection check) and the
TM-score of Zhang & Skolnick (2004), with ``d0 = 1.24 (L - 15)^(1/3) - 1.8`` and the floor of 0.5 A for ``L <= 21``.

Oracles: the d0 values and score sums below are worked out by hand (the cube root of 8, 27, 64 and 125 is 2, 3, 4 and
5, so d0 is exact for L = 23, 42, 79 and 140); identical structures give RMSD 0 and TM-score 1 in any pose; a copy with
one residue moved 50 A away scores ``(29 + 1 / (1 + (50 / d0)^2)) / 30``, because the other 29 residues superpose
exactly; and the score is compared with a numerical optimiser and with an independent Kabsch implementation.
"""

from __future__ import annotations

from typing import Any, Dict

import numpy as np
import pytest

from metainformant.protein.structure import general
from metainformant.protein.structure.analysis import (
    calculate_structural_alignment_quality,
    calculate_structural_similarity,
)
from metainformant.protein.structure.general import compute_rmsd_kabsch


def chain_like(n: int, seed: int, step: float = 3.8) -> np.ndarray:
    """A smooth random walk with ``step`` between neighbours, like the C-alpha trace of a chain."""
    rng = np.random.default_rng(seed)
    vectors = rng.normal(size=(n, 3))
    vectors /= np.linalg.norm(vectors, axis=1, keepdims=True)
    directions = np.cumsum(vectors, axis=0) * 0.35 + vectors
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    return np.vstack([np.zeros(3), np.cumsum(directions[:-1] * step, axis=0)])


def rotation_from_seed(seed: int) -> np.ndarray:
    q, r = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def pose(coords: np.ndarray, seed: int) -> np.ndarray:
    """``coords`` rotated and translated by a seeded rigid motion."""
    shift = np.random.default_rng(seed + 100).uniform(-80, 80, size=3)
    return coords @ rotation_from_seed(seed).T + shift


def ca_structure(coords: np.ndarray, chain: str = "A", first: int = 1) -> Dict[str, Any]:
    return {
        "atoms": [
            {
                "record_name": "ATOM",
                "name": "CA",
                "res_name": "ALA",
                "chain_id": chain,
                "res_seq": first + i,
                "x": float(c[0]),
                "y": float(c[1]),
                "z": float(c[2]),
            }
            for i, c in enumerate(coords)
        ]
    }


# ---------------------------------------------------------------------------
# d0 and the score formula, by hand
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "length, d0",
    [
        (23, 1.24 * 2 - 1.8),  # (23 - 15) = 8 = 2^3
        (42, 1.24 * 3 - 1.8),  # 27 = 3^3
        (79, 1.24 * 4 - 1.8),  # 64 = 4^3
        (140, 1.24 * 5 - 1.8),  # 125 = 5^3
    ],
)
def test_d0_follows_the_published_formula_exactly_where_the_cube_root_is_whole(length: int, d0: float) -> None:
    assert general.tm_score_d0(length) == pytest.approx(d0, abs=1e-12)


def test_d0_is_floored_at_half_an_angstrom_for_up_to_21_residues() -> None:
    assert [general.tm_score_d0(n) for n in (1, 5, 15, 16, 20, 21)] == [
        0.5
    ] * 6  # the formula is negative or below 0.5 there
    assert general.tm_score_d0(22) == pytest.approx(1.24 * 7 ** (1 / 3) - 1.8)  # 0.572: the formula takes over
    assert general.tm_score_d0(22) > 0.5
    assert general.tm_score_d0(100) == pytest.approx(3.6521, abs=1e-4)
    with pytest.raises(ValueError):
        general.tm_score_d0(0)


def test_the_score_of_given_distances_is_the_hand_computed_sum() -> None:
    d0 = general.tm_score_d0(23)  # 0.68
    # terms 1/(1 + (d/d0)^2): 1, 1/2, 1/5
    assert general.tm_score_from_distances([0.0, d0, 2 * d0], 23) == pytest.approx((1 + 0.5 + 0.2) / 23, abs=1e-12)
    assert general.tm_score_from_distances([0.0] * 10, 40) == pytest.approx(0.25)  # unaligned residues count as zero
    assert general.tm_score_from_distances([], 40) == 0.0
    assert general.tm_score_from_distances([0.0] * 11 + [general.tm_score_d0(22)] * 11, 22) == pytest.approx(
        0.75, abs=1e-12
    )


def test_distances_that_cannot_be_scored_are_rejected() -> None:
    with pytest.raises(ValueError):
        general.tm_score_from_distances([1.0, -1.0], 10)
    with pytest.raises(ValueError):
        general.tm_score_from_distances([1.0, float("nan")], 10)
    with pytest.raises(ValueError, match="cannot exceed"):
        general.tm_score_from_distances([1.0] * 11, 10)


# ---------------------------------------------------------------------------
# TM-score of superposable structures
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [3, 10, 25, 80])
def test_identical_structures_have_tm_score_one_in_any_pose(n: int) -> None:
    ref = chain_like(n, seed=n)

    result = general.calculate_tm_score(ref, pose(ref, seed=n + 1))

    assert result["tm_score"] == pytest.approx(1.0, abs=1e-9)
    assert (result["n_aligned"], result["n_reference"]) == (n, n)
    assert result["d0"] == pytest.approx(general.tm_score_d0(n))


def test_the_superposition_is_returned_and_is_a_proper_rotation() -> None:
    ref = chain_like(40, seed=2)
    mobile = pose(ref + np.random.default_rng(3).normal(scale=0.4, size=ref.shape), seed=5)

    result = general.calculate_tm_score(ref, mobile)
    rotation, translation = result["rotation"], result["translation"]

    assert rotation @ rotation.T == pytest.approx(np.eye(3), abs=1e-10)
    assert np.linalg.det(rotation) == pytest.approx(1.0, abs=1e-10)
    superposed = mobile @ rotation + translation
    assert general.tm_score_from_distances(np.linalg.norm(superposed - ref, axis=1), 40) == pytest.approx(
        result["tm_score"], abs=1e-12
    )


def test_the_score_is_invariant_to_rotating_or_translating_either_structure() -> None:
    ref = chain_like(60, seed=11)
    mobile = ref + np.random.default_rng(12).normal(scale=1.2, size=ref.shape)
    base = general.calculate_tm_score(ref, mobile)
    base_rmsd = compute_rmsd_kabsch(ref, mobile)

    for seed in (1, 2, 3):
        moved = general.calculate_tm_score(pose(ref, seed), pose(mobile, seed + 10))
        assert moved["tm_score"] == pytest.approx(base["tm_score"], abs=1e-8)
        assert compute_rmsd_kabsch(pose(ref, seed), pose(mobile, seed + 10)) == pytest.approx(base_rmsd, abs=1e-8)


def test_one_residue_moved_far_away_costs_exactly_its_own_term() -> None:
    """The 29 other residues superpose exactly, so TM = (29 + 1 / (1 + (D / d0)^2)) / 30 by hand."""
    ref = chain_like(30, seed=4)
    mobile = ref.copy()
    mobile[7] += np.array([50.0, 0.0, 0.0])
    d0 = 1.24 * 15 ** (1 / 3) - 1.8  # L = 30

    result = general.calculate_tm_score(ref, mobile)

    assert result["tm_score"] == pytest.approx((29 + 1 / (1 + (50.0 / d0) ** 2)) / 30, abs=1e-7)
    assert result["rotation"] == pytest.approx(np.eye(3), abs=1e-6)
    assert result["translation"] == pytest.approx(np.zeros(3), abs=1e-6)


def test_residues_without_a_partner_lower_the_score_by_their_share() -> None:
    ref = chain_like(30, seed=6)

    full = general.calculate_tm_score(ref, pose(ref, 3))
    half = general.calculate_tm_score(ref, pose(ref, 3), n_reference=60)

    assert full["tm_score"] == pytest.approx(1.0, abs=1e-9)
    assert half["tm_score"] == pytest.approx(0.5, abs=1e-9)  # 30 perfect pairs of 60 reference residues
    assert half["d0"] == pytest.approx(general.tm_score_d0(60))
    assert (half["n_aligned"], half["n_reference"]) == (30, 60)


def test_the_score_never_exceeds_the_aligned_fraction() -> None:
    ref = chain_like(50, seed=8)
    mobile = ref + np.random.default_rng(9).normal(scale=0.5, size=ref.shape)

    assert 0.0 < general.calculate_tm_score(ref, mobile, n_reference=80)["tm_score"] <= 50 / 80


# ---------------------------------------------------------------------------
# The search, noise, partial similarity, mirror images
# ---------------------------------------------------------------------------


def test_more_noise_means_a_lower_score_and_a_larger_rmsd() -> None:
    ref = chain_like(70, seed=14)
    scores, rmsds = [], []
    for scale in (0.3, 1.0, 2.5, 6.0):
        mobile = ref + np.random.default_rng(15).normal(size=ref.shape) * scale
        scores.append(general.calculate_tm_score(ref, mobile)["tm_score"])
        rmsds.append(compute_rmsd_kabsch(ref, mobile))

    assert scores == sorted(scores, reverse=True) and len(set(scores)) == 4
    assert rmsds == sorted(rmsds)
    assert scores[0] > 0.9 and scores[-1] < 0.35


def test_the_search_finds_the_conserved_core_that_the_rmsd_superposition_misses() -> None:
    """Half of the residues are identical, half are scattered: the best superposition fits the identical half."""
    rng = np.random.default_rng(10)
    ref = rng.normal(size=(100, 3)) * 10
    mobile = ref.copy()
    mobile[50:] += rng.normal(size=(50, 3)) * 30

    result = general.calculate_tm_score(ref, mobile)
    # the score of the superposition that minimises the RMSD over all 100 residues, for comparison
    from metainformant.protein.structure.general import align_structures_kabsch

    aligned, _, _ = align_structures_kabsch(ref, mobile)
    rmsd_superposition_score = general.tm_score_from_distances(np.linalg.norm(aligned - ref, axis=1), 100)

    assert 0.5 <= result["tm_score"] < 0.6  # the 50 identical residues give 0.5, the rest only a little
    assert result["tm_score"] > rmsd_superposition_score + 0.1


def test_the_score_is_close_to_what_a_numerical_optimiser_reaches() -> None:
    from scipy.optimize import minimize
    from scipy.spatial.transform import Rotation

    ref = chain_like(60, seed=42)
    mobile = pose(ref + np.random.default_rng(43).normal(scale=1.0, size=ref.shape), seed=7)
    result = general.calculate_tm_score(ref, mobile)
    d0 = result["d0"]

    def negative_score(x: np.ndarray) -> float:
        moved = mobile @ Rotation.from_rotvec(x[:3]).as_matrix().T + x[3:]
        return float(-np.sum(1.0 / (1.0 + (np.linalg.norm(moved - ref, axis=1) / d0) ** 2)) / 60)

    start = np.concatenate([Rotation.from_matrix(result["rotation"].T).as_rotvec(), result["translation"]])
    refined = minimize(
        negative_score, start, method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-12, "maxiter": 5000}
    )

    assert -refined.fun >= result["tm_score"] - 1e-9  # the TM-score of an actual superposition is a lower bound ...
    assert -refined.fun - result["tm_score"] < 3e-3  # ... and the search lands within 0.3% of the local optimum


def test_a_mirror_image_does_not_superpose_because_reflections_are_not_allowed() -> None:
    k = np.arange(20)
    helix = np.column_stack([2.3 * np.cos(np.radians(100) * k), 2.3 * np.sin(np.radians(100) * k), 1.5 * k])
    mirror = helix * np.array([-1.0, 1.0, 1.0])

    result = general.calculate_tm_score(helix, mirror)

    assert np.linalg.det(result["rotation"]) == pytest.approx(1.0, abs=1e-10)
    assert result["tm_score"] < 0.5
    assert compute_rmsd_kabsch(helix, mirror) > 1.0  # the best proper rotation leaves a large residual
    turned = helix @ rotation_from_seed(3).T
    assert compute_rmsd_kabsch(helix, turned) < 1e-8  # whereas a rotated copy superposes exactly


def test_unrelated_structures_score_low() -> None:
    rng = np.random.default_rng(31)
    first = rng.normal(size=(100, 3)) * 10
    second = rng.normal(size=(100, 3)) * 10

    assert general.calculate_tm_score(first, second)["tm_score"] < 0.3


def test_the_kabsch_rmsd_agrees_with_an_independent_implementation() -> None:
    pytest.importorskip("Bio")
    from Bio.SVDSuperimposer import SVDSuperimposer

    ref = chain_like(80, seed=50)
    mobile = pose(ref + np.random.default_rng(51).normal(scale=1.5, size=ref.shape), seed=3)
    superimposer = SVDSuperimposer()
    superimposer.set(ref, mobile)
    superimposer.run()

    assert compute_rmsd_kabsch(ref, mobile) == pytest.approx(superimposer.get_rms(), abs=1e-8)


def test_invalid_arrays_are_rejected() -> None:
    ref = chain_like(10, seed=1)
    with pytest.raises(ValueError, match="same shape"):
        general.calculate_tm_score(ref, ref[:9])
    with pytest.raises(ValueError, match="at least 3"):
        general.calculate_tm_score(ref[:2], ref[:2])
    with pytest.raises(ValueError, match="finite"):
        general.calculate_tm_score(ref, np.where(np.arange(10)[:, None] == 3, np.nan, ref))
    with pytest.raises(ValueError, match="shape"):
        general.calculate_tm_score(ref[:, :2], ref[:, :2])
    with pytest.raises(ValueError, match="n_reference"):
        general.calculate_tm_score(ref, ref, n_reference=9)


# ---------------------------------------------------------------------------
# Structure dictionaries
# ---------------------------------------------------------------------------


def test_the_same_structure_in_another_pose_has_rmsd_zero_and_tm_one() -> None:
    """The former code returned the RMSD of the unsuperposed coordinates, tens of angstrom here."""
    coords = chain_like(50, seed=21)
    first, second = ca_structure(coords), ca_structure(pose(coords, 9))

    result = calculate_structural_similarity(first, second)

    assert result["rmsd"] == pytest.approx(0.0, abs=1e-8)
    assert result["tm_score"] == pytest.approx(1.0, abs=1e-9)
    assert (result["n_aligned"], result["n_residues_1"], result["n_residues_2"]) == (50, 50, 50)
    assert result["d0"] == pytest.approx(general.tm_score_d0(50))
    assert "similarity" not in result  # the arbitrary 1 / (1 + rmsd / 10) is gone


def test_noise_lowers_the_tm_score_and_raises_the_rmsd_of_structures() -> None:
    coords = chain_like(70, seed=22)
    noisy_coords = coords + np.random.default_rng(23).normal(size=coords.shape) * 1.5

    similarity = calculate_structural_similarity(ca_structure(coords), ca_structure(noisy_coords))

    assert similarity["rmsd"] == pytest.approx(compute_rmsd_kabsch(coords, noisy_coords))
    assert similarity["tm_score"] == pytest.approx(general.calculate_tm_score(coords, noisy_coords)["tm_score"])
    assert 0.5 < similarity["tm_score"] < 1.0 and similarity["rmsd"] > 1.0


def test_residues_missing_from_the_second_structure_count_as_unaligned() -> None:
    coords = chain_like(30, seed=24)
    second = ca_structure(coords)
    second["atoms"] = [atom for i, atom in enumerate(second["atoms"]) if not 10 <= i < 20]  # residues 11-20 are missing

    result = calculate_structural_similarity(ca_structure(coords), second)

    assert (result["n_aligned"], result["n_residues_1"], result["n_residues_2"]) == (20, 30, 20)
    assert result["rmsd"] == pytest.approx(0.0, abs=1e-8)
    assert result["tm_score"] == pytest.approx(
        20 / 30, abs=1e-9
    )  # 20 perfect pairs, normalised by the 30 residues of structure1


def test_pairing_by_order_ignores_the_labels_and_by_residue_id_needs_them() -> None:
    coords = chain_like(40, seed=25)
    first = ca_structure(coords, chain="A", first=1)
    renumbered = ca_structure(pose(coords, 4), chain="B", first=101)

    with pytest.raises(ValueError, match="no residue in common"):
        calculate_structural_similarity(first, renumbered)
    by_order = calculate_structural_similarity(first, renumbered, pair_by="order")

    assert by_order["tm_score"] == pytest.approx(1.0, abs=1e-9)
    with pytest.raises(ValueError, match="equally long"):
        calculate_structural_similarity(first, ca_structure(coords[:30]), pair_by="order")
    with pytest.raises(ValueError, match="pair_by"):
        calculate_structural_similarity(first, first, pair_by="alignment")


def test_chain_and_insertion_code_are_part_of_the_residue_label() -> None:
    coords = chain_like(30, seed=26)
    first = ca_structure(coords)
    second = ca_structure(coords, chain="B")  # the same numbers in another chain: nothing is shared
    with pytest.raises(ValueError, match="no residue in common"):
        calculate_structural_similarity(first, second)

    with_icode = ca_structure(coords)
    with_icode["atoms"][5]["icode"] = "A"  # residue 6 becomes 6A: it no longer pairs with residue 6
    assert calculate_structural_similarity(first, with_icode)["n_aligned"] == 29


def test_only_c_alpha_atoms_of_protein_residues_are_compared() -> None:
    coords = chain_like(30, seed=27)
    first = ca_structure(coords)
    padded = ca_structure(pose(coords, 6))
    for atom in list(padded["atoms"]):
        padded["atoms"].append({**atom, "name": "CB", "x": atom["x"] + 1.0})  # side chain atoms
    padded["atoms"].append(
        {
            "record_name": "HETATM",
            "name": "CA",
            "res_name": "CA",
            "chain_id": "A",
            "res_seq": 1,
            "x": 9.0,
            "y": 9.0,
            "z": 9.0,
        }
    )

    result = calculate_structural_similarity(first, padded)

    assert result["n_residues_2"] == 30
    assert result["tm_score"] == pytest.approx(1.0, abs=1e-9)


def test_bare_coordinates_are_read_as_a_c_alpha_trace() -> None:
    coords = chain_like(30, seed=28)

    result = calculate_structural_similarity({"coordinates": coords}, {"coordinates": pose(coords, 2)})

    assert result["tm_score"] == pytest.approx(1.0, abs=1e-9)
    assert result["n_aligned"] == 30


def test_structures_that_cannot_be_compared_are_an_error_not_a_made_up_score() -> None:
    """The former code returned ``{"rmsd": inf, "similarity": 0.0}`` for a missing structure."""
    good = ca_structure(chain_like(30, seed=29))
    with pytest.raises(ValueError, match="structure2 has no coordinates"):
        calculate_structural_similarity(good, {})
    with pytest.raises(ValueError, match="structure1 has no coordinates"):
        calculate_structural_similarity({"atoms": []}, good)
    with pytest.raises(ValueError, match="no C-alpha"):
        calculate_structural_similarity(
            good,
            {
                "atoms": [
                    {"record_name": "ATOM", "name": "N", "chain_id": "A", "res_seq": 1, "x": 0.0, "y": 0.0, "z": 0.0}
                ]
            },
        )
    with pytest.raises(ValueError, match="at least 3"):
        calculate_structural_similarity(ca_structure(chain_like(2, seed=1)), ca_structure(chain_like(2, seed=1)))


def test_the_alignment_quality_is_the_tm_score_and_no_longer_a_rescaled_rmsd() -> None:
    coords = chain_like(60, seed=30)
    first = ca_structure(coords)
    moved = ca_structure(pose(coords, 8))
    noisy = ca_structure(coords + np.random.default_rng(32).normal(size=coords.shape) * 2.0)

    perfect = calculate_structural_alignment_quality(first, moved)
    imperfect = calculate_structural_alignment_quality(first, noisy)
    similarity = calculate_structural_similarity(first, noisy)

    assert perfect["quality_score"] == pytest.approx(
        1.0, abs=1e-9
    )  # max(0, 1 - rmsd / 20) of the unsuperposed pose was lower
    assert imperfect["quality_score"] == imperfect["tm_score"] == pytest.approx(similarity["tm_score"])
    assert imperfect["rmsd"] == pytest.approx(similarity["rmsd"])
    assert imperfect["quality_score"] != pytest.approx(max(0.0, 1.0 - imperfect["rmsd"] / 20.0), abs=1e-3)
    assert (imperfect["aligned_atoms"], imperfect["structure1_atoms"], imperfect["structure2_atoms"]) == (60, 60, 60)
    assert imperfect["d0"] == pytest.approx(general.tm_score_d0(60))


def test_the_old_arbitrary_scalings_are_not_what_is_reported() -> None:
    coords = chain_like(40, seed=33)
    noisy_coords = coords + np.random.default_rng(34).normal(size=coords.shape) * 1.0
    result = calculate_structural_similarity(ca_structure(coords), ca_structure(noisy_coords))
    unsuperposed = float(np.sqrt(np.mean(np.sum((coords - noisy_coords) ** 2, axis=1))))

    assert result["rmsd"] < unsuperposed + 1e-12  # superposition can only lower the RMSD
    assert result["tm_score"] != pytest.approx(1.0 / (1.0 + result["rmsd"] / 10.0), abs=1e-3)
    assert set(result) == {"rmsd", "tm_score", "d0", "n_aligned", "n_residues_1", "n_residues_2"}


def test_structures_of_hundreds_of_residues_are_scored() -> None:
    coords = chain_like(300, seed=35)
    noisy = coords + np.random.default_rng(36).normal(scale=1.0, size=coords.shape)

    result = general.calculate_tm_score(coords, noisy)

    assert 0.5 < result["tm_score"] < 1.0
    assert result["d0"] == pytest.approx(general.tm_score_d0(300))
