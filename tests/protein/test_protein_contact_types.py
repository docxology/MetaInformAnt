"""Tests for contact classification on hand-built structures in real PDB format.

Atoms are written as fixed-column PDB ``ATOM`` records and parsed with the
repository's own ``parse_pdb_atoms``, so the contact functions run on the same
atom dictionaries and coordinate arrays they receive in production.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.protein.structure.contacts import (
    classify_contact_types,
    identify_disulfide_bonds,
    identify_hydrogen_bonds,
    identify_hydrophobic_contacts,
    identify_salt_bridges,
)
from metainformant.protein.structure.pdb import parse_pdb_atoms

# (atom name, residue name, residue number, (x, y, z) in angstrom, element)
Record = tuple[str, str, int, tuple[float, float, float], str]

# Geometry of CONTACT_RICH, centre-to-centre distances in angstrom:
#   LYS 1 NZ (0, 0, 0)    ASP 1 carboxylate centre (3.25, 0.5, 0)  -> salt bridge, hypot(3.25, 0.5)
#   LEU 3 sidechain centre (20.25, 0.5, 0)  VAL 4 sidechain centre (23.25, 0.5, 0)  -> 3.0
#   CYS 5 SG (40, 0, 0)   CYS 6 SG (42.05, 0, 0)                      -> disulfide, 2.05
#   NZ (0, 0, 0)  OD1 (3, 0, 0)                                       -> hydrogen bond, 3.0
CONTACT_RICH: list[Record] = [
    ("N", "LYS", 1, (-3.0, 0.0, 0.0), "N"),
    ("CA", "LYS", 1, (-2.0, 0.0, 0.0), "C"),
    ("NZ", "LYS", 1, (0.0, 0.0, 0.0), "N"),
    ("CA", "ASP", 2, (8.0, 0.0, 0.0), "C"),
    ("OD1", "ASP", 2, (3.0, 0.0, 0.0), "O"),
    ("OD2", "ASP", 2, (3.5, 1.0, 0.0), "O"),
    ("CB", "LEU", 3, (20.0, 0.0, 0.0), "C"),
    ("CG", "LEU", 3, (20.5, 1.0, 0.0), "C"),
    ("CB", "VAL", 4, (23.0, 0.0, 0.0), "C"),
    ("CG1", "VAL", 4, (23.5, 1.0, 0.0), "C"),
    ("SG", "CYS", 5, (40.0, 0.0, 0.0), "S"),
    ("SG", "CYS", 6, (42.05, 0.0, 0.0), "S"),
]


def _structure(records: list[Record]) -> tuple[list[dict], np.ndarray]:
    """Parse fixed-column PDB records into the atom dictionaries and coordinate array."""
    lines = [
        (
            f"ATOM  {serial:5d} {name:<4s} {res_name:>3s} A{res_seq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2s}"
        )
        for serial, (name, res_name, res_seq, (x, y, z), element) in enumerate(records, start=1)
    ]
    atoms = parse_pdb_atoms("\n".join(lines))
    coords = np.array([[atom["x"], atom["y"], atom["z"]] for atom in atoms])
    return atoms, coords


def test_salt_bridge_between_lysine_and_aspartate_uses_group_centres():
    atoms, coords = _structure(CONTACT_RICH)

    bridges = identify_salt_bridges(atoms, coords)

    assert len(bridges) == 1
    bridge = bridges[0]
    assert bridge["positive_residue"] == ("A", 1, "LYS")
    assert bridge["negative_residue"] == ("A", 2, "ASP")
    assert bridge["distance"] == pytest.approx(math.hypot(3.25, 0.5))
    assert bridge["positive_atoms"] == [2]
    assert bridge["negative_atoms"] == [4, 5]


def test_salt_bridge_beyond_threshold_is_not_reported():
    far_apart = [
        ("NZ", "LYS", 1, (0.0, 0.0, 0.0), "N"),
        ("OD1", "ASP", 2, (6.0, 0.0, 0.0), "O"),
        ("OD2", "ASP", 2, (6.5, 1.0, 0.0), "O"),
    ]
    atoms, coords = _structure(far_apart)

    assert identify_salt_bridges(atoms, coords) == []
    assert len(identify_salt_bridges(atoms, coords, distance_threshold=7.0)) == 1


def test_hydrophobic_contact_between_leucine_and_valine_sidechains():
    atoms, coords = _structure(CONTACT_RICH)

    contacts = identify_hydrophobic_contacts(atoms, coords)

    assert len(contacts) == 1
    assert contacts[0]["residue1"] == ("A", 3, "LEU")
    assert contacts[0]["residue2"] == ("A", 4, "VAL")
    assert contacts[0]["distance"] == pytest.approx(3.0)


def test_disulfide_bond_between_cysteine_sulfurs():
    atoms, coords = _structure(CONTACT_RICH)

    bonds = identify_disulfide_bonds(atoms, coords)

    assert len(bonds) == 1
    assert bonds[0]["residue1"] == ("A", 5)
    assert bonds[0]["residue2"] == ("A", 6)
    assert bonds[0]["distance"] == pytest.approx(2.05)


def test_hydrogen_bond_between_nitrogen_donor_and_oxygen_acceptor():
    atoms, coords = _structure(CONTACT_RICH)

    bonds = identify_hydrogen_bonds(atoms, coords)

    assert [(b["donor_atom"]["name"], b["acceptor_atom"]["name"]) for b in bonds] == [("NZ", "OD1")]
    assert bonds[0]["distance"] == pytest.approx(3.0)


def test_classify_contact_types_returns_lists_map_and_summary():
    atoms, coords = _structure(CONTACT_RICH)

    result = classify_contact_types(atoms, coords)

    assert set(result) == {
        "residue_contacts",
        "hydrogen_bonds",
        "salt_bridges",
        "hydrophobic_contacts",
        "disulfide_bonds",
        "summary",
    }
    # Only the contact lists are lists: the residue map is an array and the summary a dict of counts.
    for key in ("hydrogen_bonds", "salt_bridges", "hydrophobic_contacts", "disulfide_bonds"):
        assert isinstance(result[key], list) and len(result[key]) == 1
    contact_map = result["residue_contacts"]
    assert isinstance(contact_map, np.ndarray) and contact_map.shape == (6, 6)
    assert np.array_equal(contact_map, contact_map.T)
    assert result["summary"] == {
        "total_h_bonds": 1,
        "total_salt_bridges": 1,
        "total_hydrophobic": 1,
        "total_disulfides": 1,
        "total_residue_contacts": 3,
    }
