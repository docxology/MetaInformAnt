"""Regression tests: AlphaFold confidence figures are per residue, as documented, not per atom.

Defect: ``parse_alphafold_confidence`` documented "confidence scores per residue" but returned one B-factor per ATOM
record. AlphaFold writes a residue's pLDDT into every one of its atoms, so a tryptophan (14 heavy atoms) counted
fourteen times and a glycine (4) four times. ``get_alphafold_structure_quality`` then averaged over atoms
(weighting residues by size) and reported atom counts as ``high/medium/low_confidence_residues`` and in its
``confidence_distribution``; ``validate_alphafold_structure`` did the same for ``avg_confidence`` and never filled
``n_residues`` (always 0). A residue is now identified by chain, residue number and insertion code and counted once.

Oracles: a small model is built by hand with residues of different sizes, and every expected value is plain
arithmetic on the per-residue scores written into it.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import pytest

from metainformant.protein.structure.alphafold import (
    get_alphafold_structure_quality,
    parse_alphafold_confidence,
    validate_alphafold_structure,
)
from metainformant.protein.workflow.orchestration import assess_alphafold_quality

GLY = ["N", "CA", "C", "O"]  # 4 atoms
TRP = ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]  # 14 atoms
ALA = ["N", "CA", "C", "O", "CB"]  # 5 atoms
SER = ["N", "CA", "C", "O", "CB", "OG"]  # 6 atoms

Residue = Tuple[str, Sequence[str], float]  # residue name, atom names, B-factor of every atom


def atom_line(
    serial: int,
    name: str,
    residue: str,
    chain: str,
    number: int,
    b_factor: Optional[float],
    insertion_code: str = " ",
    record: str = "ATOM  ",
) -> str:
    """One PDB coordinate line with the fixed columns of the format (B-factor in columns 61-66)."""
    b_field = f"{b_factor:6.2f}" if b_factor is not None else "      "
    return (
        f"{record}{serial:5d} {' ' + name:<4s} {residue:>3s} {chain}{number:4d}{insertion_code}   "
        f"{1.0 * serial:8.3f}{0.5 * serial:8.3f}{0.25 * serial:8.3f}{1.0:6.2f}{b_field}          {name[0]:>2s}  "
    )


def write_model(path: Path, residues: Sequence[Residue], chain: str = "A", first_number: int = 1) -> Path:
    lines = []
    serial = 0
    for number, (residue, atoms, b_factor) in enumerate(residues, start=first_number):
        for name in atoms:
            serial += 1
            lines.append(atom_line(serial, name, residue, chain, number, b_factor))
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


# Four residues of very different size. Per-residue mean pLDDT (95 + 30 + 75.5 + 55) / 4 = 63.875; the mean over the
# 29 atoms would be (4*95 + 14*30 + 5*75.5 + 6*55) / 29 = 51.98.
RESIDUES: List[Residue] = [("GLY", GLY, 95.0), ("TRP", TRP, 30.0), ("ALA", ALA, 75.5), ("SER", SER, 55.0)]
ATOM_WEIGHTED_MEAN = (4 * 95.0 + 14 * 30.0 + 5 * 75.5 + 6 * 55.0) / 29


# ---------------------------------------------------------------------------
# parse_alphafold_confidence
# ---------------------------------------------------------------------------


def test_one_confidence_score_per_residue_in_file_order(tmp_path: Path) -> None:
    model = write_model(tmp_path / "model.pdb", RESIDUES)

    scores = parse_alphafold_confidence(model)

    assert scores == [95.0, 30.0, 75.5, 55.0]  # four residues, not 29 atoms
    assert len(scores) != 29


def test_a_residue_is_not_counted_once_per_atom(tmp_path: Path) -> None:
    few = write_model(tmp_path / "few.pdb", [("GLY", ["CA"], 80.0), ("ALA", ["CA"], 40.0)])
    many = write_model(tmp_path / "many.pdb", [("GLY", GLY * 3, 80.0), ("ALA", ALA * 5, 40.0)])

    assert parse_alphafold_confidence(few) == parse_alphafold_confidence(many) == [80.0, 40.0]


def test_residues_with_the_same_number_in_different_chains_stay_apart(tmp_path: Path) -> None:
    lines = []
    serial = 0
    for chain, score in (("A", 90.0), ("B", 60.0)):
        for number in (1, 2):
            for name in GLY:
                serial += 1
                lines.append(atom_line(serial, name, "GLY", chain, number, score - number))
    path = tmp_path / "chains.pdb"
    path.write_text("\n".join(lines) + "\n")

    assert parse_alphafold_confidence(path) == [89.0, 88.0, 59.0, 58.0]  # chain A residues 1, 2 then chain B


def test_insertion_codes_make_residues_distinct(tmp_path: Path) -> None:
    lines = []
    serial = 0
    for number, code, score in ((52, " ", 70.0), (52, "A", 71.0), (52, "B", 72.0), (53, " ", 73.0)):
        for name in ALA:
            serial += 1
            lines.append(atom_line(serial, name, "ALA", "A", number, score, insertion_code=code))
    path = tmp_path / "insertions.pdb"
    path.write_text("\n".join(lines) + "\n")

    assert parse_alphafold_confidence(path) == [70.0, 71.0, 72.0, 73.0]


def test_only_atom_records_of_the_first_model_count(tmp_path: Path) -> None:
    model = write_model(tmp_path / "model.pdb", RESIDUES)
    extra = [
        atom_line(900, "O", "HOH", "A", 301, 5.0, record="HETATM"),  # water: not a residue of the model
        "MODEL        2",
        atom_line(901, "CA", "ALA", "A", 1, 10.0),  # a second model repeating residue 1 with another value
        "ENDMDL",
    ]
    lines = model.read_text().splitlines()
    text = "\n".join(lines[:-1] + ["ENDMDL"] + extra + ["END"]) + "\n"
    model.write_text(text)

    assert parse_alphafold_confidence(model) == [95.0, 30.0, 75.5, 55.0]


def test_atoms_with_different_b_factors_are_averaged_per_residue(tmp_path: Path) -> None:
    """Not an AlphaFold model (experimental B-factors vary by atom): the residue value is the atom mean."""
    lines = [atom_line(i + 1, name, "ALA", "A", 1, b) for i, (name, b) in enumerate(zip(ALA, (10, 20, 30, 40, 50)))]
    lines += [atom_line(10 + i, name, "GLY", "A", 2, b) for i, (name, b) in enumerate(zip(GLY, (150, 250, 250, 150)))]
    path = tmp_path / "varying.pdb"
    path.write_text("\n".join(lines) + "\n")

    assert parse_alphafold_confidence(path) == [30.0, 200.0]  # no 0-100 restriction in this function


def test_a_missing_file_and_a_malformed_record_give_an_empty_list_and_a_log_entry(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    broken = tmp_path / "broken.pdb"
    broken.write_text(atom_line(1, "CA", "ALA", "A", 1, 50.0) + "\n" + atom_line(2, "CA", "GLY", "A", 2, None) + "\n")

    with caplog.at_level(logging.WARNING):
        assert parse_alphafold_confidence(tmp_path / "absent.pdb") == []
        assert parse_alphafold_confidence(broken) == []  # not a silently truncated list of the good residue

    assert "not found" in caplog.text and "no numeric B-factor" in caplog.text


# ---------------------------------------------------------------------------
# get_alphafold_structure_quality
# ---------------------------------------------------------------------------


def test_quality_is_computed_per_residue(tmp_path: Path) -> None:
    model = write_model(tmp_path / "model.pdb", RESIDUES)

    quality = get_alphafold_structure_quality(model)

    assert quality["plddt_score"] == pytest.approx(63.875)
    assert quality["plddt_score"] != pytest.approx(ATOM_WEIGHTED_MEAN, abs=1.0)
    # Residues >= 90: GLY; >= 70: ALA; below 70: TRP and SER. (The atom counts would have been 4, 5 and 20.)
    assert quality["high_confidence_residues"] == 1
    assert quality["medium_confidence_residues"] == 1
    assert quality["low_confidence_residues"] == 2
    # 0-50: TRP (30); 50-70: SER (55); 70-90: ALA (75.5); 90-100: GLY (95)
    assert quality["confidence_distribution"] == {"very_low": 1, "low": 1, "medium": 1, "high": 1}


def test_quality_does_not_change_when_a_residue_gains_atoms(tmp_path: Path) -> None:
    small = write_model(tmp_path / "small.pdb", RESIDUES)
    bigger = write_model(
        tmp_path / "bigger.pdb", [("GLY", GLY + ["X1", "X2", "X3", "X4", "X5", "X6"], 95.0)] + RESIDUES[1:]
    )

    assert get_alphafold_structure_quality(small) == get_alphafold_structure_quality(bigger)


def test_the_counts_add_up_to_the_number_of_residues(tmp_path: Path) -> None:
    model = write_model(tmp_path / "model.pdb", RESIDUES * 3)  # 12 residues, 87 atoms

    quality = get_alphafold_structure_quality(model)

    assert (
        quality["high_confidence_residues"] + quality["medium_confidence_residues"] + quality["low_confidence_residues"]
    ) == 12
    assert sum(quality["confidence_distribution"].values()) == 12


def test_quality_of_a_missing_file_is_still_the_empty_report(tmp_path: Path) -> None:
    quality = get_alphafold_structure_quality(tmp_path / "absent.pdb")

    assert quality["plddt_score"] == 0.0
    assert quality["confidence_distribution"] == {}
    assert quality["high_confidence_residues"] == quality["medium_confidence_residues"] == 0


def test_the_quality_rating_of_a_model_follows_the_residue_mean(tmp_path: Path) -> None:
    """``assess_alphafold_quality`` rates by ``plddt_score``: a small confident residue next to a large 60 one.

    Per residue the mean is (95 + 60) / 2 = 77.5, rated "confident" (70 or more); weighted by atoms it was
    (4*95 + 14*60) / 18 = 67.8, which rated the same model "low".
    """
    model = write_model(tmp_path / "model.pdb", [("GLY", GLY, 95.0), ("TRP", TRP, 60.0)])

    report = assess_alphafold_quality(model)

    assert report["confidence"]["plddt_score"] == pytest.approx(77.5)
    assert report["quality_rating"] == "confident"
    assert (4 * 95.0 + 14 * 60.0) / 18 < 70.0  # the atom-weighted figure would have fallen below the threshold


# ---------------------------------------------------------------------------
# validate_alphafold_structure
# ---------------------------------------------------------------------------


def test_validation_reports_residues_and_the_residue_mean_confidence(tmp_path: Path) -> None:
    model = write_model(tmp_path / "model.pdb", RESIDUES)

    validation = validate_alphafold_structure(model)

    assert validation["is_valid"] is True
    assert validation["n_atoms"] == 29
    assert validation["n_residues"] == 4  # was always 0
    assert validation["has_confidence_scores"] is True
    assert validation["avg_confidence"] == pytest.approx(63.875)
    assert validation["issues"] == []


def test_validation_counts_residues_across_chains_and_insertion_codes(tmp_path: Path) -> None:
    lines = [
        atom_line(1, "CA", "ALA", "A", 1, 80.0),
        atom_line(2, "CB", "ALA", "A", 1, 80.0),
        atom_line(3, "CA", "ALA", "B", 1, 60.0),
        atom_line(4, "CA", "ALA", "A", 1, 70.0, insertion_code="A"),
    ]
    path = tmp_path / "multi.pdb"
    path.write_text("\n".join(lines) + "\n")

    validation = validate_alphafold_structure(path)

    assert (validation["n_atoms"], validation["n_residues"]) == (4, 3)
    assert validation["avg_confidence"] == pytest.approx((80.0 + 60.0 + 70.0) / 3)


def test_validation_of_atoms_without_b_factors_still_counts_them(tmp_path: Path) -> None:
    path = tmp_path / "no_b.pdb"
    path.write_text(
        "\n".join([atom_line(1, "CA", "ALA", "A", 1, None), atom_line(2, "CA", "GLY", "A", 2, None)]) + "\n"
    )

    validation = validate_alphafold_structure(path)

    assert validation["is_valid"] is True
    assert (validation["n_atoms"], validation["n_residues"]) == (2, 2)
    assert validation["has_confidence_scores"] is False
    assert validation["avg_confidence"] == 0.0


def test_validation_problems_are_reported_as_issues(tmp_path: Path) -> None:
    empty = tmp_path / "empty.pdb"
    empty.write_text("HEADER    NOTHING\nEND\n")
    outside = write_model(tmp_path / "outside.pdb", [("ALA", ALA, 150.0), ("GLY", GLY, 180.0)])

    assert validate_alphafold_structure(empty)["issues"] == ["No ATOM records found"]
    assert validate_alphafold_structure(empty)["is_valid"] is False
    assert validate_alphafold_structure(tmp_path / "absent.pdb")["issues"] == ["PDB file not found"]
    assert validate_alphafold_structure(outside)["issues"] == ["Confidence scores outside expected range"]
