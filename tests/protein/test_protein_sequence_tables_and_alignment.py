"""Tests: sequence properties use complete published tables, and alignment results are what they say they are.

Defects:

- ``instability_index`` and ``function.prediction_core._compute_instability_index`` each carried an 80-entry "subset of
  the most impactful pairs" of the Guruprasad dipeptide table: 16 of the entries had wrong values and 118 dipeptides
  whose true weight is not 1.0 were silently weighted 1.0, so the index of a 45-residue test sequence was 41.27
  where the published table gives 52.56. The complete 400-entry table is used, from one place.
- ``molecular_weight`` summed residue masses and left out the water of the termini (a tripeptide of glycine came out
  at 171.16 Da instead of 189.17), while ``compute_physicochemical`` added it.
- The ER retrieval signal pattern ``KKXX$`` contained a literal ``X`` and could never match.
- ``read_fasta`` kept only the last of two records with the same ID.
- ``matrix_align`` reported ``"matrix": "BLOSUM62"`` whatever matrix it was given.
- ``multi_sequence_alignment`` claimed a guide tree, aligned each new sequence to the first row only, inserted gaps into
  the other rows at positions of the ungapped first sequence (so rows no longer spelled their sequences), and returned
  the rows in the order of insertion, not of the input.

Oracles: Biopython's ProtParam (instability index, molecular weight) and its copy of the Guruprasad table when it is
installed, hand-computed sums from the published weights, and invariants of a multiple alignment (each row without its
gaps is its input sequence, all rows have one length).
"""

from __future__ import annotations

from pathlib import Path
from typing import List

import numpy as np
import pytest

from metainformant.protein.function.prediction_core import (
    _compute_instability_index,
    compute_physicochemical,
    predict_localization,
)
from metainformant.protein.sequence import sequences as sequence_module
from metainformant.protein.sequence.alignment import BLOSUM62, matrix_align, multi_sequence_alignment
from metainformant.protein.sequence.sequences import instability_index, molecular_weight, read_fasta, write_fasta

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
SEQUENCES = [
    "MKTIIALSYIFCLVFADYKDDDDKAGQEELKRSTPPWHHNCGGMM",  # the sequence of the defect description
    "ACDEFGHIKLMNPQRSTVWY" * 3,
    "GGG",
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",  # ubiquitin
]


# ---------------------------------------------------------------------------
# Instability index
# ---------------------------------------------------------------------------


def test_the_dipeptide_table_is_complete() -> None:
    DIWV = sequence_module.DIWV
    assert sorted(DIWV) == sorted(a + b for a in AMINO_ACIDS for b in AMINO_ACIDS)  # all 400 ordered pairs
    values = list(DIWV.values())
    assert values.count(1.0) == 247
    assert sorted(set(values)) == [
        -15.91,
        -14.03,
        -9.37,
        -7.49,
        -6.54,
        -1.88,
        1.0,
        13.34,
        18.38,
        20.26,
        24.64,
        24.68,
        33.6,
        33.601,
        44.94,
        58.28,
    ]


def test_the_dipeptide_table_equals_the_one_of_biopython() -> None:
    pytest.importorskip("Bio")
    from Bio.SeqUtils.ProtParamData import DIWV as reference

    DIWV = sequence_module.DIWV
    assert all(DIWV[a + b] == reference[a][b] for a in AMINO_ACIDS for b in AMINO_ACIDS)


def test_the_index_is_ten_over_length_times_the_sum_of_the_pair_weights() -> None:
    assert instability_index("WM") == pytest.approx(10 / 2 * 24.68)  # one pair, WM = 24.68
    assert instability_index("AA") == pytest.approx(10 / 2 * 1.0)
    assert instability_index("ACD") == pytest.approx(10 / 3 * (44.94 + 20.26))  # AC = 44.94, CD = 20.26
    assert instability_index("W") == 0.0 and instability_index("") == 0.0


@pytest.mark.parametrize("sequence", SEQUENCES)
def test_the_index_equals_biopythons(sequence: str) -> None:
    pytest.importorskip("Bio")
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    expected = ProteinAnalysis(sequence).instability_index()

    assert instability_index(sequence) == pytest.approx(expected, abs=1e-9)
    assert _compute_instability_index(sequence) == pytest.approx(expected, abs=1e-9)


def test_random_sequences_agree_with_biopython_and_the_old_subset_does_not_matter_any_more() -> None:
    pytest.importorskip("Bio")
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    rng = np.random.default_rng(3)
    for _ in range(25):
        sequence = "".join(rng.choice(list(AMINO_ACIDS), size=int(rng.integers(5, 300))))
        assert instability_index(sequence) == pytest.approx(ProteinAnalysis(sequence).instability_index(), abs=1e-9)


def test_the_value_of_the_defect_report() -> None:
    sequence = SEQUENCES[0]

    assert instability_index(sequence) == pytest.approx(52.5578, abs=1e-3)  # the published table; the subset gave 41.27
    assert compute_physicochemical(sequence)["instability_index"] == pytest.approx(52.56, abs=0.006)
    assert (
        compute_physicochemical(sequence)["is_stable"] is False
    )  # > 40; the subset's 41.27 also said so, but not by this much


def test_pairs_with_a_non_standard_residue_count_as_neutral() -> None:
    assert instability_index("AXA") == pytest.approx(10 / 3 * 2.0)  # AX and XA have no weight: 1.0 each


# ---------------------------------------------------------------------------
# Molecular weight
# ---------------------------------------------------------------------------


def test_molecular_weight_includes_the_water_of_the_termini() -> None:
    assert molecular_weight("G") == pytest.approx(57.0519 + 18.01528)  # free glycine, 75.07 Da
    assert molecular_weight("GGG") == pytest.approx(3 * 57.0519 + 18.01528)
    assert molecular_weight("GGG") > 189.0


@pytest.mark.parametrize("sequence", SEQUENCES)
def test_molecular_weight_equals_biopythons(sequence: str) -> None:
    pytest.importorskip("Bio")
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    assert molecular_weight(sequence) == pytest.approx(ProteinAnalysis(sequence).molecular_weight(), rel=2e-5)


@pytest.mark.parametrize("sequence", SEQUENCES)
def test_the_two_molecular_weights_of_the_package_agree(sequence: str) -> None:
    assert compute_physicochemical(sequence)["molecular_weight"] == pytest.approx(molecular_weight(sequence), abs=0.02)


def test_nothing_is_added_to_a_sequence_without_residues() -> None:
    assert molecular_weight("") == 0.0
    assert molecular_weight("XXX") == 0.0  # no residue of known mass: not the weight of a water molecule


# ---------------------------------------------------------------------------
# ER retrieval signal
# ---------------------------------------------------------------------------


def test_a_c_terminal_dilysine_signal_is_found() -> None:
    sequence = "M" + "A" * 30 + "KKAA"

    signals = predict_localization(sequence)["signals_found"]

    assert [s["type"] for s in signals if s["type"] == "ER_retention"] == ["ER_retention"]
    er = next(s for s in signals if s["type"] == "ER_retention")
    assert (er["sequence"], er["position"]) == ("KKAA", len(sequence) - 4)
    assert predict_localization(sequence)["scores"]["endoplasmic_reticulum"] > 0


def test_a_sequence_without_the_signal_has_none() -> None:
    for sequence in ("M" + "A" * 30 + "KAKA", "M" + "A" * 30 + "AAKK", "M" + "A" * 34):
        assert all(s["type"] != "ER_retention" for s in predict_localization(sequence)["signals_found"]), sequence


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------


def test_two_records_with_one_id_are_an_error_not_a_lost_sequence(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.fasta"
    path.write_text(">p1 first\nACDE\n>p2\nGGGG\n>p1 again\nWWWW\n")

    with pytest.raises(ValueError, match="Duplicate sequence ID 'p1'"):
        read_fasta(path)


def test_distinct_ids_are_all_kept(tmp_path: Path) -> None:
    path = tmp_path / "ok.fasta"
    write_fasta({"a": "ACDE", "b": "GGGG", "c": "WWWW"}, path)

    assert read_fasta(path) == {"a": "ACDE", "b": "GGGG", "c": "WWWW"}


# ---------------------------------------------------------------------------
# matrix_align labels its matrix
# ---------------------------------------------------------------------------


def test_the_reported_matrix_is_the_one_used() -> None:
    custom = {a: {b: (3 if a == b else -3) for b in AMINO_ACIDS} for a in AMINO_ACIDS}

    for mode in ("global", "local"):
        assert matrix_align("ACDEF", "ACDEF", mode=mode)["matrix"] == "BLOSUM62"
        assert matrix_align("ACDEF", "ACDEF", matrix=BLOSUM62, mode=mode)["matrix"] == "BLOSUM62"
        result = matrix_align("ACDEF", "ACDEF", matrix=custom, mode=mode)
        assert result["matrix"] == "custom"
        assert result["score"] == 15  # five matches of 3 each: the custom matrix was really used
    assert matrix_align("", "ACD", matrix=custom)["matrix"] == "custom"  # also on the early return for an empty input


# ---------------------------------------------------------------------------
# Multiple sequence alignment
# ---------------------------------------------------------------------------


def sum_of_pairs(rows: List[str]) -> int:
    """Independent sum-of-pairs score: +1 match, -1 mismatch, -2 against a gap, 0 for two gaps."""
    total = 0
    for column in zip(*rows):
        for a in range(len(column)):
            for b in range(a + 1, len(column)):
                if column[a] == "-" and column[b] == "-":
                    continue
                if "-" in (column[a], column[b]):
                    total -= 2
                else:
                    total += 1 if column[a] == column[b] else -1
    return total


def assert_valid_alignment(sequences: List[str], result: dict) -> None:
    rows = result["aligned_sequences"]
    assert [row.replace("-", "") for row in rows] == sequences  # every row spells its input sequence, in input order
    assert len({len(row) for row in rows}) == 1 and result["alignment_length"] == len(rows[0])
    assert result["n_sequences"] == len(sequences)
    assert all(any(character != "-" for character in column) for column in zip(*rows))  # no column of only gaps
    assert result["score"] == sum_of_pairs(rows)


def test_a_missing_residue_becomes_a_gap_in_its_row_only() -> None:
    result = multi_sequence_alignment(["ACDEFG", "ACDFG", "ACDEFG"])

    assert result["aligned_sequences"] == ["ACDEFG", "ACD-FG", "ACDEFG"]
    assert result["score"] == 12  # five columns of three matches (3 each), and E/-/E: -2 + 1 - 2 = -3
    assert (result["center"], result["method"], result["alignment_length"]) == (0, "center_star", 6)


@pytest.mark.parametrize(
    "sequences",
    [
        ["ACGT", "AGT", "ACT", "ACCGT"],
        ["MKTAYIAKQR", "MKTAYIAKQRQISFVKSH", "AYIAKQRQISF", "MKTAKQR"],
        ["ACD", "ACD", "ACD"],
        ["AAAAAAAA", "CCCCCCCC", "AAAAAAAA"],
        ["A", "AC", "ACD", "ACDE"],
        ["WWWACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWYWWW", "ACDEFGHIKLMNPQRSTVWY"],
    ],
)
def test_every_row_spells_its_sequence_and_the_rows_are_in_input_order(sequences: List[str]) -> None:
    assert_valid_alignment(sequences, multi_sequence_alignment(sequences))


def test_insertions_in_different_sequences_are_all_kept() -> None:
    """Sequences 1 and 2 each carry an insertion the other lacks: the gaps of both stay in every row."""
    sequences = ["ACDEFGHIK", "ACDXXEFGHIK", "ACDEFGYYHIK"]

    result = multi_sequence_alignment(sequences)

    assert_valid_alignment(sequences, result)
    assert result["alignment_length"] == len("ACDXXEFGYYHIK")


def test_random_sequences_give_valid_alignments() -> None:
    rng = np.random.default_rng(11)
    base = "".join(rng.choice(list(AMINO_ACIDS), size=40))
    sequences = []
    for _ in range(6):
        residues = list(base)
        for _ in range(4):
            position = int(rng.integers(0, len(residues)))
            if rng.random() < 0.5:
                del residues[position]
            else:
                residues.insert(position, str(rng.choice(list(AMINO_ACIDS))))
        sequences.append("".join(residues))

    assert_valid_alignment(sequences, multi_sequence_alignment(sequences))


def test_trivial_inputs_and_unknown_methods() -> None:
    assert multi_sequence_alignment([]) == {
        "aligned_sequences": [],
        "score": 0,
        "n_sequences": 0,
        "alignment_length": 0,
        "method": "center_star",
        "center": None,
    }
    single = multi_sequence_alignment(["ACDE"])
    assert single["aligned_sequences"] == ["ACDE"] and single["n_sequences"] == 1 and single["alignment_length"] == 4
    with pytest.raises(ValueError, match="center_star"):
        multi_sequence_alignment(["ACD", "ACE"], method="progressive")
