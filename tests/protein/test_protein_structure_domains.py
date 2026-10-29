"""Tests: domains are found from the residue contact graph, not from gaps in the residue numbers.

Defect: ``identify_domains`` started a new domain wherever the chain identifier changed or the residue number jumped
by more than 10. It never looked at a coordinate, so two compact folds joined by a linker were one "domain", while a
single compact fold with a deletion in its numbering was two. It now cuts the contact graph of the C-alpha atoms by
recursive spectral bisection (Fiedler vector of the normalised Laplacian, best normalised cut).

Oracles: planted structures whose domains are known by construction -- compact groups of points at protein-like
density (minimum separation 3.8 A, drawn with a seeded generator) joined by a linker laid out in a line, a single
compact group, an inserted domain whose residue numbers interrupt the host domain, and lattice blocks big enough to
use the sparse eigensolver. The contact counts are recomputed by brute force.
"""

from __future__ import annotations

from typing import Any, Dict, List, Sequence

import numpy as np
import pytest

from metainformant.protein.structure.analysis import identify_domains

SURFACE_PADDING = 1.0


def blob(n: int, radius: float, rng: np.random.Generator, center: Sequence[float] = (0.0, 0.0, 0.0)) -> np.ndarray:
    """``n`` points in a ball, at least 3.8 A apart (a compact fold at protein-like density)."""
    points: List[np.ndarray] = []
    attempts = 0
    while len(points) < n:
        attempts += 1
        assert attempts < 200_000, "could not place the points"
        direction = rng.normal(size=3)
        candidate = direction / np.linalg.norm(direction) * radius * rng.random() ** (1.0 / 3.0)
        if all(np.linalg.norm(candidate - other) >= 3.8 for other in points):
            points.append(candidate)
    return np.array(points) + np.asarray(center, dtype=float)


def linker(start_x: float, end_x: float, n: int) -> np.ndarray:
    """``n`` points evenly spaced on the x axis strictly between ``start_x`` and ``end_x``."""
    x = start_x + (end_x - start_x) * (np.arange(1, n + 1) / (n + 1))
    return np.column_stack([x, np.zeros(n), np.zeros(n)])


def residues(coords: np.ndarray, chain: str = "A", numbers: Sequence[int] | None = None, **extra: Any):
    """C-alpha atom dictionaries (the keys of ``parse_pdb_file``)."""
    numbers = list(range(1, len(coords) + 1)) if numbers is None else list(numbers)
    return [
        {
            "record_name": "ATOM",
            "serial": i + 1,
            "name": "CA",
            "res_name": "ALA",
            "chain_id": chain,
            "res_seq": int(number),
            "x": float(c[0]),
            "y": float(c[1]),
            "z": float(c[2]),
            "element": "C",
            **extra,
        }
        for i, (number, c) in enumerate(zip(numbers, coords))
    ]


def two_groups(seed: int, linker_length: int = 10, n_first: int = 80, n_second: int = 80) -> np.ndarray:
    """Two compact groups of 13.5 A radius whose surfaces are bridged by ``linker_length`` residues."""
    rng = np.random.default_rng(seed)
    radius = 13.5
    gap = 3.8 * linker_length
    first = blob(n_first, radius, rng)
    second = blob(n_second, radius, rng, center=(2 * radius + gap, 0.0, 0.0))
    return np.vstack([first, linker(radius, radius + gap, linker_length), second])


def index_sets(domains: List[Dict[str, Any]]) -> List[List[int]]:
    return [d["residues"] for d in domains]


def brute_force_contacts(coords: np.ndarray, members: Sequence[int], cutoff: float = 8.0) -> int:
    chosen = coords[[m - 1 for m in members]]
    distances = np.linalg.norm(chosen[:, None, :] - chosen[None, :, :], axis=2)
    return int((np.triu(distances <= cutoff, k=1)).sum())


def rotation_from_seed(seed: int) -> np.ndarray:
    q, r = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


# ---------------------------------------------------------------------------
# Planted domains are found
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("seed, linker_length", [(1, 6), (2, 10), (3, 14), (4, 10)])
def test_two_compact_groups_joined_by_a_linker_are_two_domains(seed: int, linker_length: int) -> None:
    coords = two_groups(seed, linker_length)
    n_total = len(coords)

    domains = identify_domains({"atoms": residues(coords)})

    assert len(domains) == 2
    first, second = domains
    assert (first["domain_id"], second["domain_id"]) == (1, 2)
    assert first["chain_id"] == second["chain_id"] == "A"
    # every residue is in exactly one domain, and the cut lies in the linker
    assert sorted(first["residues"] + second["residues"]) == list(range(1, n_total + 1))
    assert set(range(1, 81)) <= set(first["residues"])
    assert set(range(80 + linker_length + 1, n_total + 1)) <= set(second["residues"])
    boundary = first["end"]
    assert 80 <= boundary <= 80 + linker_length
    assert first["residues"] == list(range(1, boundary + 1)) and second["residues"] == list(
        range(boundary + 1, n_total + 1)
    )
    assert first["segments"] == [{"start": 1, "end": boundary}]
    assert second["segments"] == [{"start": boundary + 1, "end": n_total}]
    assert (first["size"], second["size"]) == (len(first["residues"]), len(second["residues"]))
    assert (first["start"], second["end"]) == (1, n_total)


def test_the_contacts_of_a_domain_are_counted_inside_it() -> None:
    coords = two_groups(2, 10)

    first, second = identify_domains({"atoms": residues(coords)})

    assert first["n_contacts"] == brute_force_contacts(coords, first["residues"])
    assert second["n_contacts"] == brute_force_contacts(coords, second["residues"])
    assert first["n_contacts"] > 100 and second["n_contacts"] > 100


@pytest.mark.parametrize("n, radius, seed", [(90, 14.0, 5), (120, 15.5, 6)])
def test_a_single_compact_group_is_one_domain(n: int, radius: float, seed: int) -> None:
    coords = blob(n, radius, np.random.default_rng(seed))

    domains = identify_domains({"atoms": residues(coords)})

    assert len(domains) == 1
    assert domains[0]["residues"] == list(range(1, n + 1))
    assert domains[0]["segments"] == [{"start": 1, "end": n}]


def test_three_groups_in_a_row_are_three_domains() -> None:
    rng = np.random.default_rng(9)
    radius, gap = 13.5, 3.8 * 8
    pieces = []
    for index in range(3):
        center_x = index * (2 * radius + gap)
        pieces.append(blob(70, radius, rng, center=(center_x, 0.0, 0.0)))
        if index < 2:
            pieces.append(linker(center_x + radius, center_x + radius + gap, 8))
    coords = np.vstack(pieces)

    domains = identify_domains({"atoms": residues(coords)})

    assert [d["domain_id"] for d in domains] == [1, 2, 3]
    assert [len(d["segments"]) for d in domains] == [1, 1, 1]
    assert set(range(1, 71)) <= set(domains[0]["residues"])
    assert set(range(79, 149)) <= set(domains[1]["residues"])
    assert set(range(157, 227)) <= set(domains[2]["residues"])
    assert sum(d["size"] for d in domains) == len(coords)


def test_a_domain_inserted_into_another_makes_the_host_discontinuous() -> None:
    """Host domain: residues 1-30 and 71-100; inserted domain: residues 31-70, packed against it with a gap."""
    rng = np.random.default_rng(3)
    host = blob(60, 12.5, rng)
    inserted = blob(40, 10.5, rng, center=(12.5 + 10.5 + 6.0, 0.0, 0.0))
    coords = np.vstack([host[:30], inserted, host[30:]])

    domains = identify_domains({"atoms": residues(coords)})

    assert len(domains) == 2
    host_domain = next(d for d in domains if d["size"] == 60)
    inserted_domain = next(d for d in domains if d["size"] == 40)
    assert host_domain["segments"] == [{"start": 1, "end": 30}, {"start": 71, "end": 100}]
    assert inserted_domain["segments"] == [{"start": 31, "end": 70}]
    assert (host_domain["start"], host_domain["end"]) == (1, 100)
    assert host_domain["residues"] == list(range(1, 31)) + list(range(71, 101))


def test_residue_numbers_play_no_part_in_the_partition() -> None:
    coords = two_groups(2, 10)
    plain = identify_domains({"atoms": residues(coords)})
    spread = identify_domains({"atoms": residues(coords, numbers=[1000 + 37 * i for i in range(len(coords))])})

    assert [d["size"] for d in spread] == [d["size"] for d in plain]
    for renumbered, original in zip(spread, plain):
        assert renumbered["residues"] == [1000 + 37 * (r - 1) for r in original["residues"]]


def test_a_gap_in_the_numbering_of_a_compact_fold_does_not_make_two_domains() -> None:
    """The former rule (a jump of more than 10) split this single fold in two."""
    coords = blob(90, 14.0, np.random.default_rng(5))
    numbers = list(range(1, 46)) + list(range(200, 245))  # a 154-residue gap in the numbering

    domains = identify_domains({"atoms": residues(coords, numbers=numbers)})

    assert len(domains) == 1
    assert domains[0]["size"] == 90 and (domains[0]["start"], domains[0]["end"]) == (1, 244)


def test_the_result_does_not_depend_on_position_or_orientation() -> None:
    coords = two_groups(3, 14)
    moved = coords @ rotation_from_seed(7).T + np.array([55.0, -20.0, 310.0])

    assert index_sets(identify_domains({"atoms": residues(moved)})) == index_sets(
        identify_domains({"atoms": residues(coords)})
    )


def test_each_chain_is_decomposed_on_its_own_even_where_chains_touch() -> None:
    coords = two_groups(2, 10)
    atoms = residues(coords, chain="A") + residues(coords + np.array([0.0, 4.0, 0.0]), chain="B")  # B lies on A

    domains = identify_domains({"atoms": atoms})

    assert [d["chain_id"] for d in domains] == ["A", "A", "B", "B"]
    assert [d["domain_id"] for d in domains] == [1, 2, 3, 4]
    assert [d["size"] for d in domains[:2]] == [d["size"] for d in domains[2:]]


def test_a_small_disconnected_piece_joins_the_domain_nearest_to_it() -> None:
    rng = np.random.default_rng(11)
    first = blob(80, 13.5, rng)
    second = blob(80, 13.5, rng, center=(70.0, 0.0, 0.0))
    fragment = np.array([[70.0 + 13.5 + 25.0 + 3.8 * i, 0.0, 0.0] for i in range(10)])  # 25 A beyond the second group

    domains = identify_domains({"atoms": residues(np.vstack([first, second, fragment]))})

    assert [d["size"] for d in domains] == [80, 90]
    assert domains[1]["residues"] == list(range(81, 171))


def test_only_c_alpha_atoms_of_protein_residues_define_the_graph() -> None:
    coords = two_groups(2, 10)
    atoms = residues(coords)
    extras = []
    for atom in atoms[:60]:
        extras.append({**atom, "name": "CB", "x": atom["x"] + 1.5})  # side-chain atoms: not on the C-alpha trace
        extras.append({**atom, "alt_loc": "B", "y": atom["y"] + 0.3})  # a second alternate location of the same residue
    extras.append({**atoms[0], "record_name": "HETATM", "res_name": "CA", "res_seq": 900, "x": 500.0})  # calcium ion
    structure = {"atoms": atoms + extras}

    assert index_sets(identify_domains(structure)) == index_sets(identify_domains({"atoms": atoms}))


def test_insertion_codes_keep_residues_with_the_same_number_apart() -> None:
    coords = blob(90, 14.0, np.random.default_rng(5))
    atoms = residues(coords)
    atoms[45]["res_seq"] = atoms[44]["res_seq"]
    atoms[45]["icode"] = "A"

    domains = identify_domains({"atoms": atoms})

    assert len(domains) == 1 and domains[0]["size"] == 90
    assert domains[0]["residue_ids"][44:46] == ["45", "45A"]
    assert domains[0]["residues"][44:46] == [45, 45]


def test_both_pdb_loaders_give_the_same_domains(tmp_path) -> None:
    from metainformant.protein.structure.io import parse_pdb_file
    from metainformant.protein.structure.pdb import load_pdb_file

    coords = two_groups(2, 10)
    lines = [
        f"ATOM  {i + 1:5d}  CA  ALA A{i + 1:4d}    {c[0]:8.3f}{c[1]:8.3f}{c[2]:8.3f}{1.0:6.2f}{20.0:6.2f}           C  "
        for i, c in enumerate(coords)
    ]
    path = tmp_path / "two_groups.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")

    from_io = identify_domains(parse_pdb_file(path))
    from_pdb = identify_domains(load_pdb_file(path))

    assert [d["residues"] for d in from_io] == [d["residues"] for d in from_pdb]
    assert len(from_io) == 2


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------


def test_the_parameters_control_how_readily_domains_are_split() -> None:
    coords = two_groups(2, 10)
    atoms = residues(coords)

    assert len(identify_domains({"atoms": atoms}, max_normalized_cut=0.0)) == 1  # the linker is not free
    assert len(identify_domains({"atoms": atoms}, min_domain_size=100)) == 1  # no part would hold 100 residues
    assert len(identify_domains({"atoms": atoms}, contact_distance=12.0)) == 2  # the groups are far from touching
    single = residues(blob(120, 15.5, np.random.default_rng(6)))
    assert len(identify_domains({"atoms": single})) == 1
    assert len(identify_domains({"atoms": single}, max_normalized_cut=2.0)) >= 2  # a cut is accepted whatever it costs


def test_a_structure_too_small_to_hold_two_domains_is_one() -> None:
    for n in (1, 2, 30, 79):
        coords = np.column_stack([3.8 * np.arange(n), np.zeros(n), np.zeros(n)])
        domains = identify_domains({"atoms": residues(coords)})
        assert len(domains) == 1 and domains[0]["size"] == n


@pytest.mark.parametrize(
    "kwargs",
    [
        {"contact_distance": 0.0},
        {"contact_distance": float("nan")},
        {"min_domain_size": 1},
        {"min_domain_size": 10.5},
        {"max_normalized_cut": -0.1},
        {"max_normalized_cut": 2.5},
    ],
)
def test_invalid_parameters_are_rejected(kwargs: Dict[str, Any]) -> None:
    with pytest.raises(ValueError):
        identify_domains({"atoms": residues(two_groups(2, 10))}, **kwargs)


def test_structures_without_usable_input() -> None:
    assert identify_domains({"atoms": []}) == []
    assert identify_domains({}) == []
    ligand_only = [
        {
            "record_name": "HETATM",
            "name": "C1",
            "res_name": "LIG",
            "chain_id": "A",
            "res_seq": 1,
            "x": 0.0,
            "y": 0.0,
            "z": 0.0,
        }
    ]
    assert identify_domains({"atoms": ligand_only}) == []  # no protein residues, so no domains

    with pytest.raises(ValueError, match="no 'atoms' records"):
        identify_domains({"coordinates": np.zeros((5, 3))})
    with pytest.raises(ValueError, match="x, y and z"):
        identify_domains({"atoms": [{"name": "CA", "chain_id": "A", "res_seq": 1}]})


def test_residue_numbers_alone_define_nothing_any_more() -> None:
    """The structure of the former test: five residues with a jump of more than 10 and no coordinates."""
    atoms = [{"chain_id": "A", "res_seq": number} for number in (1, 2, 3, 50, 51)]

    assert identify_domains({"atoms": atoms}) == []  # no C-alpha coordinates, so nothing to partition


# ---------------------------------------------------------------------------
# Large graphs use the sparse eigensolver
# ---------------------------------------------------------------------------


def lattice_block(nx: int, ny: int, nz: int, z0: float) -> np.ndarray:
    return np.array([[3.8 * i, 3.8 * j, z0 + 3.8 * k] for i in range(nx) for j in range(ny) for k in range(nz)])


def test_a_structure_of_nearly_two_thousand_residues_is_split_at_its_neck() -> None:
    first = lattice_block(10, 10, 9, 0.0)
    second = lattice_block(10, 10, 9, 3.8 * 17)
    neck = np.array([[0.0, 0.0, 3.8 * (9 + i)] for i in range(8)])
    coords = np.vstack([first, neck, second])
    coords = coords + np.random.default_rng(5).normal(scale=0.3, size=coords.shape)  # break the lattice symmetry
    assert len(coords) > 1000

    domains = identify_domains({"atoms": residues(coords)})
    again = identify_domains({"atoms": residues(coords)})

    assert len(domains) == 2
    assert domains[0]["residues"] == list(range(1, domains[0]["size"] + 1))
    assert 900 <= domains[0]["size"] <= 908 and domains[0]["size"] + domains[1]["size"] == len(coords)
    assert index_sets(again) == index_sets(domains)  # the sparse solver starts from a fixed vector: reproducible
