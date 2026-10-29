"""Tests: structure files, statistics and contacts report what the data hold and invent nothing.

Defects:

- ``calculate_pdb_statistics`` reported ``residue_types`` as the number of atoms of a residue name divided by 5
  ("Approximate"), which is wrong for every residue that does not have five atoms.
- ``pdb.validate_pdb_file`` read the first 10 kB and reported the atoms in that part as ``n_atoms`` of the file.
- ``get_pdb_sequence`` gave an "X" to every water and ligand, collapsed residues 52 and 52A and sorted by number
  instead of keeping file order.
- ``io._parse_atom_line`` turned a blank chain into "A" and blank coordinates or residue numbers into 0, so a damaged
  line became an atom at the origin; ``_format_atom_line`` wrote an alanine C-alpha at the origin for any field
  that the atom dictionary lacked.
- ``parse_mmcif_file`` read fixed column positions whatever the header said and dropped rows silently; its partner
  ``convert_pdb_to_mmcif`` wrote one value per line, which no CIF reader (including that one) can read back, and
  called every record ``ATOM``.
- ``identify_hydrogen_bonds`` documented an ``angle_threshold`` that it never used, returned an invented ``strength``
  and ``angle: None`` for each bond, and reported the amide nitrogen of every residue as hydrogen bonded to the carbonyl
  oxygen of the residue before it (2.3 A apart across the peptide bond), and nitrogen-oxygen pairs inside one residue.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pytest

from metainformant.protein.structure.contacts import classify_contact_types, identify_hydrogen_bonds
from metainformant.protein.structure.io import (
    convert_pdb_to_mmcif,
    parse_mmcif_file,
    parse_pdb_file,
    write_pdb_file,
)
from metainformant.protein.structure.pdb import (
    calculate_pdb_statistics,
    get_pdb_sequence,
    load_pdb_file,
    parse_pdb_atoms,
    validate_pdb_file,
)


def pdb_line(
    serial: int,
    name: str,
    res_name: str,
    number: int,
    xyz=(0.0, 0.0, 0.0),
    chain: str = "A",
    record: str = "ATOM",
    icode: str = " ",
    occupancy: float = 1.0,
    b_factor: float = 20.0,
    element: str = "C",
) -> str:
    """One fixed-column PDB coordinate record."""
    return (
        f"{record:<6}{serial:5d} {name:^4s} {res_name:>3s} {chain}{number:4d}{icode}   "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{occupancy:6.2f}{b_factor:6.2f}          {element:>2s}"
    )


# ---------------------------------------------------------------------------
# pdb.py
# ---------------------------------------------------------------------------


def residues_of(spec: Dict[str, int]) -> List[Dict[str, Any]]:
    """Atom dictionaries for residues ``{name: atoms per residue}``; each name gets three residues."""
    atoms: List[Dict[str, Any]] = []
    number = 0
    for res_name, atoms_per_residue in spec.items():
        for _ in range(3):
            number += 1
            for index in range(atoms_per_residue):
                atoms.append({"name": f"X{index}", "res_name": res_name, "chain_id": "A", "res_seq": number})
    return atoms


def test_residue_types_count_residues_not_atoms_over_five() -> None:
    atoms = residues_of({"GLY": 4, "TRP": 14})  # 12 atoms of glycine in 3 residues, 42 of tryptophan in 3

    stats = calculate_pdb_statistics({"atoms": atoms})

    assert stats["residue_types"] == {"GLY": 3, "TRP": 3}  # atoms // 5 gave 2 and 8
    assert stats["n_residues"] == 6 and stats["n_atoms"] == 54


def test_residues_with_an_insertion_code_are_distinct_residues(tmp_path: Path) -> None:
    path = tmp_path / "icode.pdb"
    path.write_text(
        "\n".join(
            [
                pdb_line(1, "CA", "GLY", 52, (0.0, 0.0, 0.0)),
                pdb_line(2, "CA", "ALA", 52, (3.8, 0.0, 0.0), icode="A"),
                pdb_line(3, "CA", "GLY", 53, (7.6, 0.0, 0.0)),
            ]
        )
    )

    data = load_pdb_file(path)

    assert data["n_residues"] == 3
    assert calculate_pdb_statistics(data)["residue_types"] == {"GLY": 2, "ALA": 1}
    assert parse_pdb_file(path)["n_residues"] == 3


def test_validation_counts_every_atom_of_a_large_file(tmp_path: Path) -> None:
    lines = ["HEADER    TEST".ljust(80)] + [pdb_line(i + 1, "CA", "ALA", i + 1, (i, 0.0, 0.0)) for i in range(600)]
    path = tmp_path / "large.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")
    assert path.stat().st_size > 40_000  # far more than the 10 kB that were read before

    result = validate_pdb_file(path)

    assert result["n_atoms"] == 600
    assert result["is_valid"] is True and result["issues"] == []


def test_validation_of_files_without_header_or_atoms(tmp_path: Path) -> None:
    no_header = tmp_path / "no_header.pdb"
    no_header.write_text(pdb_line(1, "CA", "ALA", 1) + "\n")
    no_atoms = tmp_path / "no_atoms.pdb"
    no_atoms.write_text("HEADER    NOTHING\nEND\n")

    assert validate_pdb_file(no_header)["issues"] == ["Missing HEADER record"]
    assert validate_pdb_file(no_atoms)["issues"] == ["No ATOM or HETATM records found"]
    assert validate_pdb_file(tmp_path / "missing.pdb")["is_valid"] is False


def test_the_sequence_is_that_of_the_polymer_atoms_in_file_order() -> None:
    atoms = parse_pdb_atoms(
        "\n".join(
            [
                pdb_line(1, "CA", "GLY", 1),
                pdb_line(2, "CA", "ALA", 52),
                pdb_line(3, "CA", "SER", 52, icode="A"),
                pdb_line(4, "CA", "VAL", 7),  # out of numeric order: the file order is the chain order
                pdb_line(5, "O", "HOH", 101, record="HETATM", element="O"),
                pdb_line(6, "FE", "HEM", 102, record="HETATM", element="FE"),
                pdb_line(7, "CA", "TRP", 1, chain="B"),
                pdb_line(8, "CA", "MSE", 2, chain="B", record="HETATM"),
            ]
        )
    )

    assert get_pdb_sequence({"atoms": atoms}) == "GASVW"  # no X for the water, the heme or the HETATM residue
    assert get_pdb_sequence({"atoms": atoms}, "A") == "GASV"
    assert get_pdb_sequence({"atoms": atoms}, "B") == "W"
    assert get_pdb_sequence({"atoms": []}) == ""


def test_the_sequence_reads_the_atoms_of_either_loader(tmp_path: Path) -> None:
    path = tmp_path / "two.pdb"
    path.write_text(
        "\n".join(
            [
                pdb_line(1, "CA", "GLY", 1),
                pdb_line(2, "CA", "ALA", 2, (3.8, 0.0, 0.0)),
                pdb_line(3, "O", "HOH", 9, (30.0, 0.0, 0.0), record="HETATM", element="O"),
            ]
        )
    )

    assert get_pdb_sequence(parse_pdb_file(path)) == get_pdb_sequence(load_pdb_file(path)) == "GA"


# ---------------------------------------------------------------------------
# io.py: PDB
# ---------------------------------------------------------------------------


def test_a_blank_chain_stays_blank_in_both_loaders(tmp_path: Path) -> None:
    path = tmp_path / "blank_chain.pdb"
    path.write_text(pdb_line(1, "CA", "ALA", 1, chain=" ") + "\n")

    assert parse_pdb_file(path)["atoms"][0]["chain_id"] == ""  # was "A", unlike pdb.load_pdb_file
    assert load_pdb_file(path)["atoms"][0]["chain_id"] == ""


def test_a_line_with_blank_coordinates_or_residue_number_is_skipped_not_moved_to_the_origin(tmp_path: Path) -> None:
    good = pdb_line(1, "CA", "ALA", 1, (11.0, 12.0, 13.0))
    blank_x = (
        pdb_line(2, "CA", "ALA", 2, (1.0, 2.0, 3.0))[:30] + " " * 8 + pdb_line(2, "CA", "ALA", 2, (1.0, 2.0, 3.0))[38:]
    )
    blank_number = (
        pdb_line(3, "CA", "ALA", 3, (4.0, 5.0, 6.0))[:22] + "    " + pdb_line(3, "CA", "ALA", 3, (4.0, 5.0, 6.0))[26:]
    )
    path = tmp_path / "damaged.pdb"
    path.write_text("\n".join([good, blank_x, blank_number]) + "\n")

    result = parse_pdb_file(path)

    assert result["n_atoms"] == 1
    assert result["atoms"][0]["x"] == 11.0
    assert np.count_nonzero(np.all(result["coordinates"] == 0.0, axis=1)) == 0  # nothing at (0, 0, 0)


def test_atoms_missing_required_fields_are_not_written_as_alanine_at_the_origin(tmp_path: Path) -> None:
    complete = {"name": "CA", "res_name": "GLY", "res_seq": 7, "x": 1.5, "y": -2.5, "z": 3.25, "chain_id": "B"}

    for missing in ("name", "res_name", "res_seq", "x", "y", "z"):
        incomplete = {key: value for key, value in complete.items() if key != missing}
        with pytest.raises(ValueError, match=missing):
            write_pdb_file({"atoms": [incomplete]}, tmp_path / "bad.pdb")


def test_written_atoms_read_back_unchanged_and_optional_fields_stay_blank(tmp_path: Path) -> None:
    atoms = [
        {
            "name": "CA",
            "res_name": "GLY",
            "res_seq": 7,
            "x": 1.5,
            "y": -2.5,
            "z": 3.25,
            "chain_id": "B",
            "element": "C",
        },
        {"name": "O", "res_name": "HOH", "res_seq": 9, "x": 10.0, "y": 0.5, "z": -1.0, "record_name": "HETATM"},
    ]
    path = tmp_path / "written.pdb"

    write_pdb_file({"atoms": atoms}, path)
    back = parse_pdb_file(path)["atoms"]

    assert [(a["record_name"], a["serial"], a["name"], a["res_name"], a["chain_id"], a["res_seq"]) for a in back] == [
        ("ATOM", 1, "CA", "GLY", "B", 7),
        ("HETATM", 2, "O", "HOH", "", 9),  # no chain, no element, and the serial is the position
    ]
    assert [(a["x"], a["y"], a["z"]) for a in back] == [(1.5, -2.5, 3.25), (10.0, 0.5, -1.0)]
    assert back[1]["element"] == ""


# ---------------------------------------------------------------------------
# io.py: mmCIF
# ---------------------------------------------------------------------------

CIF_SHUFFLED = """data_test
#
loop_
_atom_site.Cartn_z
_atom_site.Cartn_y
_atom_site.Cartn_x
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.pdbx_PDB_model_num
_atom_site.B_iso_or_equiv
3.0 2.0 1.0 CA ALA A 10 1 ? ATOM 1 C 1 12.5
6.0 5.0 4.0 "O5'" HOH B 20 . ? HETATM 2 O 1 30.0
9.0 8.0 7.0
  N GLY A 11 2 A ATOM 3
  N 1 8.25
30.0 20.0 10.0 CA ALA A 10 1 ? ATOM 4 C 2 99.0
#
"""


def test_mmcif_columns_are_read_by_name_whatever_their_order_quoting_or_line_breaks(tmp_path: Path) -> None:
    path = tmp_path / "shuffled.cif"
    path.write_text(CIF_SHUFFLED)

    result = parse_mmcif_file(path)
    atoms = result["atoms"]

    assert result["n_atoms"] == 3  # the second model is left out
    assert result["coordinates"].tolist() == [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
    assert [
        (a["record_name"], a["serial"], a["name"], a["res_name"], a["chain_id"], a["res_seq"], a["icode"])
        for a in atoms
    ] == [
        ("ATOM", 1, "CA", "ALA", "A", 10, ""),
        ("HETATM", 2, "O5'", "HOH", "B", 20, ""),  # the quoted name keeps its prime
        ("ATOM", 3, "N", "GLY", "A", 11, "A"),  # a row that spans two lines
    ]
    assert [a["temp_factor"] for a in atoms] == [12.5, 30.0, 8.25]
    assert math.isnan(atoms[0]["occupancy"])  # no occupancy column: unknown, not 1.0
    assert atoms[1]["label_seq_id"] is None and atoms[0]["label_seq_id"] == 1 and atoms[0]["Cartn_x"] == 1.0


def test_an_mmcif_that_cannot_be_read_is_an_error_not_an_empty_structure(tmp_path: Path) -> None:
    def write(text: str) -> Path:
        path = tmp_path / f"case{len(list(tmp_path.iterdir()))}.cif"
        path.write_text(text)
        return path

    header = "data_x\nloop_\n_atom_site.label_atom_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
    with pytest.raises(FileNotFoundError):
        parse_mmcif_file(tmp_path / "absent.cif")
    with pytest.raises(ValueError, match="No _atom_site loop"):
        parse_mmcif_file(write("data_x\nloop_\n_cell.length_a\n10.0\n"))
    with pytest.raises(ValueError, match="Cartn_x/y/z"):
        parse_mmcif_file(write("data_x\nloop_\n_atom_site.id\n_atom_site.label_atom_id\n1 CA\n"))
    with pytest.raises(ValueError, match="do not fill rows"):
        parse_mmcif_file(write(header + "CA 1.0 2.0 3.0\nN 1.0 2.0\n"))
    with pytest.raises(ValueError):  # a coordinate that is not a number
        parse_mmcif_file(write(header + "CA one 2.0 3.0\n"))
    with pytest.raises(ValueError, match="No atoms found"):
        parse_mmcif_file(write(header + "#\n"))


def test_a_pdb_file_converted_to_mmcif_reads_back_as_the_same_atoms(tmp_path: Path) -> None:
    pdb = tmp_path / "source.pdb"
    pdb.write_text(
        "\n".join(
            [
                pdb_line(1, "N", "ALA", 1, (11.104, 13.207, 9.479), element="N"),
                pdb_line(2, "CA", "ALA", 1, (12.56, 13.414, 9.567)),
                pdb_line(3, "CA", "GLY", 52, (15.56, 14.414, 11.567), icode="A", occupancy=0.5, b_factor=21.5),
                pdb_line(4, "O5'", "HOH", 7, (1.0, 2.0, 3.0), chain=" ", record="HETATM", element="O"),
            ]
        )
        + "\nEND\n"
    )
    cif = tmp_path / "converted.cif"

    convert_pdb_to_mmcif(pdb, cif)
    original = parse_pdb_file(pdb)["atoms"]
    back = parse_mmcif_file(cif)["atoms"]

    keys = (
        "record_name",
        "serial",
        "name",
        "res_name",
        "chain_id",
        "res_seq",
        "icode",
        "x",
        "y",
        "z",
        "occupancy",
        "temp_factor",
        "element",
    )
    assert [{key: atom[key] for key in keys} for atom in back] == [
        {key: atom[key] for key in keys} for atom in original
    ]
    assert [atom["record_name"] for atom in back] == ["ATOM", "ATOM", "ATOM", "HETATM"]  # was "ATOM" throughout
    text = cif.read_text()
    assert text.splitlines()[0] == "data_source"
    assert sum(1 for line in text.splitlines() if line.startswith(("ATOM", "HETATM"))) == 4  # one atom per line


def test_the_atoms_of_an_mmcif_work_with_the_analysis_functions(tmp_path: Path) -> None:
    from metainformant.protein.structure.analysis import analyze_structural_motifs

    k = np.arange(12)
    helix = np.column_stack([2.3 * np.cos(np.radians(100) * k), 2.3 * np.sin(np.radians(100) * k), 1.5 * k])
    pdb = tmp_path / "helix.pdb"
    pdb.write_text("\n".join(pdb_line(i + 1, "CA", "ALA", i + 1, tuple(c)) for i, c in enumerate(helix)) + "\n")
    cif = tmp_path / "helix.cif"
    convert_pdb_to_mmcif(pdb, cif)

    motifs = analyze_structural_motifs(parse_mmcif_file(cif))

    assert [(m["type"], m["start_residue"], m["end_residue"]) for m in motifs] == [("helix", 1, 12)]


# ---------------------------------------------------------------------------
# contacts.py: hydrogen bonds
# ---------------------------------------------------------------------------


def backbone(res_seq: int, origin, res_name: str = "ALA") -> List[Dict[str, Any]]:
    """N, CA, C and O of one residue at fixed offsets from ``origin`` (distances of an extended backbone)."""
    offsets = {"N": (0.0, 0.0, 0.0), "CA": (1.46, 0.0, 0.0), "C": (2.0, 1.4, 0.0), "O": (1.2, 2.4, 0.0)}
    return [
        {
            "name": name,
            "element": name[0],
            "res_name": res_name,
            "chain_id": "A",
            "res_seq": res_seq,
            "x": origin[0] + dx,
            "y": origin[1] + dy,
            "z": origin[2] + dz,
        }
        for name, (dx, dy, dz) in offsets.items()
    ]


def coordinates_of(atoms: List[Dict[str, Any]]) -> np.ndarray:
    return np.array([[a["x"], a["y"], a["z"]] for a in atoms])


def test_the_peptide_bond_and_a_single_residue_are_not_hydrogen_bonds() -> None:
    """N(i+1) lies 2.2 A from O(i) across the peptide bond; N(i) and O(i) of one residue are 2.8 A apart."""
    atoms = backbone(1, (0.0, 0.0, 0.0)) + backbone(2, (2.2, 2.4, 0.0))
    coords = coordinates_of(atoms)
    n2, o1 = atoms[4], atoms[3]
    assert np.linalg.norm(coords[4] - coords[3]) < 2.5  # the geometry of the trap: N(2) next to O(1)
    assert o1["name"] == "O" and n2["name"] == "N"

    assert identify_hydrogen_bonds(atoms, coords) == []  # the former code reported N2-O1 and N1-O1


def test_a_backbone_hydrogen_bond_between_distant_residues_is_found() -> None:
    donor = backbone(5, (0.0, 0.0, 0.0))
    acceptor = backbone(1, (-1.2, -2.4, 2.9))  # places O(1) at (0, 0, 2.9), 2.9 A from N(5)
    atoms = donor + acceptor
    coords = coordinates_of(atoms)

    bonds = identify_hydrogen_bonds(atoms, coords, distance_threshold=3.6)

    pairs = {
        (b["donor_atom"]["res_seq"], b["acceptor_atom"]["res_seq"], b["donor_atom"]["name"], b["acceptor_atom"]["name"])
        for b in bonds
    }
    assert (5, 1, "N", "O") in pairs
    assert all(b["distance"] <= 3.6 for b in bonds)


def test_a_side_chain_next_to_the_following_backbone_can_hydrogen_bond() -> None:
    """Only the N(i+1)-O(i) pair of the peptide bond is excluded: Asn OD1 of residue i may bind N of residue i+1."""
    atoms = [
        {"name": "OD1", "element": "O", "res_name": "ASN", "chain_id": "A", "res_seq": 4, "x": 0.0, "y": 0.0, "z": 0.0},
        {"name": "N", "element": "N", "res_name": "ALA", "chain_id": "A", "res_seq": 5, "x": 2.9, "y": 0.0, "z": 0.0},
    ]

    bonds = identify_hydrogen_bonds(atoms, coordinates_of(atoms))

    assert [(b["donor_atom"]["res_seq"], b["acceptor_atom"]["name"]) for b in bonds] == [(5, "OD1")]


def test_the_proline_nitrogen_has_no_hydrogen_to_donate() -> None:
    atoms = [
        {"name": "N", "element": "N", "res_name": "PRO", "chain_id": "A", "res_seq": 9, "x": 0.0, "y": 0.0, "z": 0.0},
        {"name": "OG", "element": "O", "res_name": "SER", "chain_id": "A", "res_seq": 3, "x": 2.9, "y": 0.0, "z": 0.0},
        {
            "name": "NE2",
            "element": "N",
            "res_name": "GLN",
            "chain_id": "A",
            "res_seq": 12,
            "x": 2.9,
            "y": 2.9,
            "z": 0.0,
        },
    ]

    bonds = identify_hydrogen_bonds(atoms, coordinates_of(atoms))

    assert [(b["donor_atom"]["name"], b["acceptor_atom"]["name"]) for b in bonds] == [("NE2", "OG")]


def test_the_bonds_carry_what_was_measured_and_nothing_invented() -> None:
    atoms = [
        {"name": "NZ", "element": "N", "res_name": "LYS", "chain_id": "A", "res_seq": 1, "x": 0.0, "y": 0.0, "z": 0.0},
        {"name": "OD1", "element": "O", "res_name": "ASP", "chain_id": "A", "res_seq": 9, "x": 0.0, "y": 0.0, "z": 3.0},
    ]

    bonds = identify_hydrogen_bonds(atoms, coordinates_of(atoms))

    assert len(bonds) == 1 and set(bonds[0]) == {"donor_atom", "acceptor_atom", "distance"}
    assert bonds[0]["distance"] == 3.0 and isinstance(bonds[0]["distance"], float)
    with pytest.raises(TypeError):  # the angle that was documented but never tested is no parameter any more
        identify_hydrogen_bonds(atoms, coordinates_of(atoms), 3.5, 120.0)  # type: ignore[call-arg]
    with pytest.raises(TypeError):
        identify_hydrogen_bonds(atoms, coordinates_of(atoms), angle_threshold=120.0)  # type: ignore[call-arg]


def test_atoms_without_an_element_column_are_read_by_name() -> None:
    """A PDB file may leave out columns 77-78; the protein atom names then give the element."""
    parsed = parse_pdb_atoms(
        "\n".join(
            [
                pdb_line(1, "NZ", "LYS", 1, (0.0, 0.0, 0.0), element="")[:76],
                pdb_line(2, "OD1", "ASP", 9, (0.0, 0.0, 3.0), element="")[:76],
            ]
        )
    )
    assert [atom["element"] for atom in parsed] == ["", ""]

    bonds = identify_hydrogen_bonds(parsed, coordinates_of(parsed))

    assert len(bonds) == 1


def test_classification_of_a_real_backbone_has_no_peptide_bond_artifacts() -> None:
    atoms = backbone(1, (0.0, 0.0, 0.0)) + backbone(2, (2.2, 2.4, 0.0)) + backbone(3, (4.4, 4.8, 0.0))

    result = classify_contact_types(atoms, coordinates_of(atoms))

    assert result["hydrogen_bonds"] == [] and result["summary"]["total_h_bonds"] == 0
