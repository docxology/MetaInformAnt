"""Regression tests: AlphaFold helpers read local model files and refuse to invent service data.

Defect: four functions of ``protein/structure/alphafold.py`` were placeholders:

* ``get_alphafold_metadata`` returned a fixed dictionary with ``confidence_score: None  # Would need to parse
  B-factor column``;
* ``find_alphafold_models_by_sequence`` and ``search_alphafold_by_keyword`` logged "not fully implemented" and
  returned ``[]``, which reads as "no model matches";
* ``get_alphafold_coverage`` returned numbers typed into the source ("Approximate as of 2024").

What can be known offline is now read from the local files (pLDDT from the B-factor column of a model, PAE from the
JSON file), and what needs the AlphaFold service raises ``NotImplementedError``. Nothing here touches the network.

Oracles: the model and PAE files are written by hand in the formats of the AlphaFold database and every expected
value is worked out below with plain arithmetic (including the residue-versus-atom weighting of the mean pLDDT).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pytest

from metainformant.protein.structure import alphafold
from metainformant.protein.structure.alphafold import (
    build_alphafold_url,
    find_alphafold_models_by_sequence,
    get_alphafold_coverage,
    get_alphafold_metadata,
    search_alphafold_by_keyword,
)

GLY_ATOMS = ["N", "CA", "C", "O"]
TRP_ATOMS = ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]
ALA_ATOMS = ["N", "CA", "C", "O", "CB"]
SER_ATOMS = ["N", "CA", "C", "O", "CB", "OG"]


def atom_line(
    serial: int, name: str, residue: str, chain: str, number: int, b_factor: float, record: str = "ATOM  "
) -> str:
    """One PDB coordinate line with the fixed columns of the format (B-factor in columns 61-66)."""
    element = name[0]
    return (
        f"{record}{serial:5d} {' ' + name:<4s} {residue:>3s} {chain}{number:4d}    "
        f"{1.0 * serial:8.3f}{0.5 * serial:8.3f}{0.25 * serial:8.3f}{1.0:6.2f}{b_factor:6.2f}          {element:>2s}  "
    )


def write_model(path: Path, residues: Sequence[Tuple[str, Sequence[str], float]], extra: Sequence[str] = ()) -> Path:
    """An AlphaFold-style model: every atom of a residue carries that residue's pLDDT as its B-factor."""
    lines = ["HEADER    ALPHAFOLD MONOMER PREDICTION", "MODEL        1"]
    serial = 0
    for number, (residue, atoms, plddt) in enumerate(residues, start=1):
        for name in atoms:
            serial += 1
            lines.append(atom_line(serial, name, residue, "A", number, plddt))
    lines.append("TER")
    lines.append("ENDMDL")
    lines.extend(extra)
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return path


# Four residues of very different size: per-residue mean (95 + 30 + 75.5 + 55) / 4 = 63.875, whereas the mean over
# all 29 atoms would be (4*95 + 14*30 + 5*75.5 + 6*55) / 29 = 51.98.
RESIDUES = [("GLY", GLY_ATOMS, 95.0), ("TRP", TRP_ATOMS, 30.0), ("ALA", ALA_ATOMS, 75.5), ("SER", SER_ATOMS, 55.0)]


# ---------------------------------------------------------------------------
# Service-dependent functions raise instead of returning made-up results
# ---------------------------------------------------------------------------


def test_the_sequence_search_is_not_implemented_instead_of_returning_no_models() -> None:
    with pytest.raises(NotImplementedError, match="find_alphafold_models_by_sequence is not implemented") as excinfo:
        find_alphafold_models_by_sequence("MKLVFEQHINDASSTWYCGP")

    assert "fetch_alphafold_model" in str(excinfo.value)  # what to do instead


def test_the_keyword_search_is_not_implemented_instead_of_returning_no_structures() -> None:
    with pytest.raises(NotImplementedError, match="search_alphafold_by_keyword is not implemented") as excinfo:
        search_alphafold_by_keyword("kinase", max_results=10)

    assert "fetch_alphafold_model" in str(excinfo.value)


def test_the_coverage_statistics_are_not_implemented_instead_of_hard_coded() -> None:
    with pytest.raises(NotImplementedError, match="get_alphafold_coverage is not implemented") as excinfo:
        get_alphafold_coverage()

    assert "not queried values" in str(excinfo.value)


def test_arguments_make_no_difference_to_the_not_implemented_functions() -> None:
    with pytest.raises(NotImplementedError):
        find_alphafold_models_by_sequence("", identity_threshold=0.0)
    with pytest.raises(NotImplementedError):
        search_alphafold_by_keyword("", max_results=0)


# ---------------------------------------------------------------------------
# get_alphafold_metadata
# ---------------------------------------------------------------------------


def test_without_files_only_what_follows_from_the_accession_is_returned() -> None:
    metadata = get_alphafold_metadata("P12345")

    assert metadata["uniprot_accession"] == "P12345"
    assert metadata["url"] == build_alphafold_url("P12345", version=4)
    assert (metadata["version"], metadata["source"]) == (4, "AlphaFold v4")
    for unavailable in ("confidence_score", "n_residues", "plddt_bands", "mean_pae", "max_pae"):
        assert metadata[unavailable] is None, unavailable


def test_the_version_decides_the_url_and_the_source_label() -> None:
    metadata = get_alphafold_metadata("Q9Y6N1", version=3)

    assert metadata["url"] == build_alphafold_url("Q9Y6N1", version=3)
    assert "model_v3" in metadata["url"]
    assert metadata["source"] == "AlphaFold v3"


def test_the_confidence_score_is_the_mean_plddt_over_residues_not_atoms(tmp_path: Path) -> None:
    model = write_model(tmp_path / "model.pdb", RESIDUES)

    metadata = get_alphafold_metadata("P12345", pdb_path=model)

    assert metadata["n_residues"] == 4
    assert metadata["confidence_score"] == pytest.approx(63.875)
    atom_weighted = (4 * 95.0 + 14 * 30.0 + 5 * 75.5 + 6 * 55.0) / 29
    assert metadata["confidence_score"] != pytest.approx(atom_weighted, abs=1.0)


def test_residues_are_counted_per_confidence_band(tmp_path: Path) -> None:
    model = write_model(tmp_path / "model.pdb", RESIDUES)

    metadata = get_alphafold_metadata("P12345", pdb_path=model)

    assert metadata["plddt_bands"] == {"very_high": 1, "confident": 1, "low": 1, "very_low": 1}


def test_band_edges_belong_to_the_higher_band(tmp_path: Path) -> None:
    residues = [("ALA", ALA_ATOMS, score) for score in (90.0, 89.99, 70.0, 69.99, 50.0, 49.99, 100.0, 0.0)]
    model = write_model(tmp_path / "edges.pdb", residues)

    bands = get_alphafold_metadata("P12345", pdb_path=model)["plddt_bands"]

    assert bands == {"very_high": 2, "confident": 2, "low": 2, "very_low": 2}
    assert sum(bands.values()) == 8


def test_only_atom_records_of_the_first_model_count(tmp_path: Path) -> None:
    water = atom_line(900, "O", "HOH", "A", 301, 5.0, record="HETATM")  # a low "B-factor" that is not pLDDT
    second_model = ["MODEL        2", atom_line(901, "CA", "ALA", "A", 1, 10.0), "ENDMDL"]
    model = write_model(tmp_path / "model.pdb", RESIDUES, extra=[water, *second_model])

    metadata = get_alphafold_metadata("P12345", pdb_path=model)

    assert metadata["n_residues"] == 4
    assert metadata["confidence_score"] == pytest.approx(63.875)


def test_a_file_that_is_not_a_model_is_rejected(tmp_path: Path) -> None:
    experimental = write_model(tmp_path / "experimental.pdb", [("ALA", ALA_ATOMS, 35.0), ("GLY", GLY_ATOMS, 180.0)])
    empty = tmp_path / "empty.pdb"
    empty.write_text("HEADER    NOTHING HERE\nEND\n")
    broken = tmp_path / "broken.pdb"
    broken.write_text(atom_line(1, "CA", "ALA", "A", 1, 50.0)[:60] + "  n/a \n")

    with pytest.raises(ValueError, match="outside 0-100"):
        get_alphafold_metadata("P12345", pdb_path=experimental)
    with pytest.raises(ValueError, match="no ATOM records"):
        get_alphafold_metadata("P12345", pdb_path=empty)
    with pytest.raises(ValueError, match="numeric B-factor"):
        get_alphafold_metadata("P12345", pdb_path=broken)


def test_a_missing_file_raises_instead_of_returning_empty_metadata(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        get_alphafold_metadata("P12345", pdb_path=tmp_path / "absent.pdb")
    with pytest.raises(FileNotFoundError):
        get_alphafold_metadata("P12345", pae_path=tmp_path / "absent.json")


# ---------------------------------------------------------------------------
# parse_alphafold_pae
# ---------------------------------------------------------------------------

# An asymmetric 3 x 3 error matrix: off-diagonal entries 1, 2, 3, 4, 5, 6 -> mean 21 / 6 = 3.5, largest value 6.
MATRIX = [[0.0, 1.0, 2.0], [3.0, 0.0, 4.0], [5.0, 6.0, 0.0]]


def write_json(path: Path, document: object) -> Path:
    path.write_text(json.dumps(document))
    return path


def sparse_document(matrix: Sequence[Sequence[float]]) -> Dict[str, List[float]]:
    """The older database layout: every (residue1, residue2, distance) triple, residues numbered from 1."""
    n = len(matrix)
    return {
        "residue1": [i + 1 for i in range(n) for _ in range(n)],
        "residue2": [j + 1 for _ in range(n) for j in range(n)],
        "distance": [matrix[i][j] for i in range(n) for j in range(n)],
        "max_predicted_aligned_error": 31.75,
    }


@pytest.mark.parametrize(
    "document",
    [
        [{"predicted_aligned_error": MATRIX, "max_predicted_aligned_error": 31.75}],  # current database layout
        {"predicted_aligned_error": MATRIX, "max_predicted_aligned_error": 31.75},  # a bare object (ColabFold)
    ],
)
def test_the_matrix_layout_is_parsed(tmp_path: Path, document: object) -> None:
    pae = alphafold.parse_alphafold_pae(write_json(tmp_path / "pae.json", document))

    np.testing.assert_array_equal(pae["pae_matrix"], np.array(MATRIX))
    assert pae["n_residues"] == 3
    assert pae["mean_pae"] == pytest.approx(3.5)
    assert pae["max_pae"] == 6.0
    assert pae["reported_max_pae"] == 31.75
    assert pae["format"] == "matrix"


def test_the_sparse_layout_gives_the_same_matrix(tmp_path: Path) -> None:
    pae = alphafold.parse_alphafold_pae(write_json(tmp_path / "pae.json", [sparse_document(MATRIX)]))

    np.testing.assert_array_equal(pae["pae_matrix"], np.array(MATRIX))
    assert (pae["mean_pae"], pae["max_pae"], pae["format"]) == (3.5, 6.0, "sparse")


def test_the_sparse_layout_does_not_depend_on_the_order_of_the_triples(tmp_path: Path) -> None:
    document = sparse_document(MATRIX)
    order = list(reversed(range(9)))
    shuffled = {key: [document[key][i] for i in order] for key in ("residue1", "residue2", "distance")}

    pae = alphafold.parse_alphafold_pae(write_json(tmp_path / "pae.json", shuffled))

    np.testing.assert_array_equal(pae["pae_matrix"], np.array(MATRIX))


def test_a_missing_reported_maximum_is_none_and_a_single_residue_has_no_mean(tmp_path: Path) -> None:
    pae = alphafold.parse_alphafold_pae(write_json(tmp_path / "pae.json", {"predicted_aligned_error": [[0.0]]}))

    assert pae["reported_max_pae"] is None
    assert pae["mean_pae"] is None
    assert (pae["n_residues"], pae["max_pae"]) == (1, 0.0)


@pytest.mark.parametrize(
    ("document", "message"),
    [
        ({"predicted_aligned_error": [[0.0, 1.0, 2.0], [3.0, 0.0, 4.0]]}, "square matrix"),
        ({"predicted_aligned_error": []}, "square matrix"),
        ({"predicted_aligned_error": [[0.0, -1.0], [1.0, 0.0]]}, "finite and non-negative"),
        ({"predicted_aligned_error": [[0.0, float("nan")], [1.0, 0.0]]}, "finite and non-negative"),
        ({"predicted_aligned_error": [[0.0, float("inf")], [1.0, 0.0]]}, "finite and non-negative"),
        ({"something_else": 1}, "neither predicted_aligned_error"),
        ([{"predicted_aligned_error": MATRIX}, {"predicted_aligned_error": MATRIX}], "expected a PAE object"),
        ("not an object", "expected a PAE object"),
        ({"residue1": [1, 1], "residue2": [1, 2], "distance": [0.0, 1.0]}, "every pair"),
        ({"residue1": [1, 1, 2, 2], "residue2": [1, 1, 2, 2], "distance": [0.0] * 4}, "more than once"),
        ({"residue1": [0, 1, 2, 2], "residue2": [1, 2, 1, 2], "distance": [0.0] * 4}, "every pair"),
        ({"residue1": [1, 2], "residue2": [1], "distance": [0.0, 1.0]}, "equal-length"),
    ],
)
def test_malformed_pae_files_are_rejected_with_a_reason(tmp_path: Path, document: object, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        alphafold.parse_alphafold_pae(write_json(tmp_path / "pae.json", document))


def test_invalid_json_and_missing_files_raise(tmp_path: Path) -> None:
    broken = tmp_path / "broken.json"
    broken.write_text("{not json")

    with pytest.raises(ValueError):
        alphafold.parse_alphafold_pae(broken)
    with pytest.raises(FileNotFoundError):
        alphafold.parse_alphafold_pae(tmp_path / "absent.json")


# ---------------------------------------------------------------------------
# Model and PAE together
# ---------------------------------------------------------------------------


def test_model_and_pae_together_fill_every_field(tmp_path: Path) -> None:
    model = write_model(tmp_path / "model.pdb", RESIDUES)
    four = [[0.0 if i == j else float(1 + i + j) for j in range(4)] for i in range(4)]
    pae = write_json(tmp_path / "pae.json", [{"predicted_aligned_error": four, "max_predicted_aligned_error": 31.75}])

    metadata = get_alphafold_metadata("P12345", pdb_path=model, pae_path=pae)

    # The 12 off-diagonal entries are 1 + i + j; the diagonal is zero. The largest off-diagonal entry is 1 + 3 + 2.
    off_diagonal = [1 + i + j for i in range(4) for j in range(4) if i != j]
    assert metadata["n_residues"] == 4
    assert metadata["confidence_score"] == pytest.approx(63.875)
    assert metadata["mean_pae"] == pytest.approx(sum(off_diagonal) / 12)
    assert metadata["max_pae"] == 6.0
    assert metadata["url"].endswith("AF-P12345-F1-model_v4.pdb")


def test_a_pae_matrix_of_another_size_than_the_model_is_rejected(tmp_path: Path) -> None:
    model = write_model(tmp_path / "model.pdb", RESIDUES)  # 4 residues
    pae = write_json(tmp_path / "pae.json", {"predicted_aligned_error": MATRIX})  # 3 residues

    with pytest.raises(ValueError, match="covers 3 residues but the model has 4"):
        get_alphafold_metadata("P12345", pdb_path=model, pae_path=pae)


def test_pae_alone_gives_the_residue_count_and_the_error_statistics(tmp_path: Path) -> None:
    pae = write_json(tmp_path / "pae.json", {"predicted_aligned_error": MATRIX})

    metadata = get_alphafold_metadata("P12345", pae_path=pae)

    assert metadata["n_residues"] == 3
    assert metadata["mean_pae"] == pytest.approx(3.5)
    assert metadata["confidence_score"] is None  # no model file: no pLDDT


def test_the_module_no_longer_carries_the_placeholder_markers() -> None:
    source = Path(alphafold.__file__).read_text().lower()

    assert "placeholder" not in source
    assert "not fully implemented" not in source
    assert "approximate as of" not in source


def test_the_optional_arguments_are_keyword_only() -> None:
    with pytest.raises(TypeError):
        get_alphafold_metadata("P12345", 4)  # type: ignore[misc]
