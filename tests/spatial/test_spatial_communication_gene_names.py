"""Ligand-receptor scoring matches genes to expression columns by name, and by nothing else.

``compute_ligand_receptor_interactions`` used to take no gene names: a numeric ligand or receptor string
selected that column, and every other name was hashed (CRC32) onto one of the columns. With real gene symbols
the score of a pair therefore came from two arbitrary columns, whatever the data contained. The function now
needs the names of the expression columns (or a dataset that carries them), scores only the database pairs
whose two genes are among them, and raises when the names are missing or no pair is present. Real computation
only; the expected scores are recomputed here by hand from the named columns.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from metainformant.spatial.communication.cell_communication import (
    compute_ligand_receptor_interactions,
    default_lr_database,
)
from metainformant.spatial.io.visium import SpatialDataset

GENES = ["ACTB", "TGFB1", "GAPDH", "WNT5A", "TGFBR1", "FZD5", "B2M", "MALAT1"]
DATABASE = {
    "pairs": [
        {"ligand": "WNT5A", "receptor": "FZD5"},  # planted: sender cells -> receiver cells
        {"ligand": "TGFB1", "receptor": "TGFBR1"},  # decoy: expressed evenly in every cell type
    ]
}
PLANTED = ("WNT5A", "FZD5", "sender", "receiver")


def _planted_expression(seed: int = 0) -> tuple[np.ndarray, list[str]]:
    """90 cells of three types over ``GENES``; WNT5A is high only in senders and FZD5 only in receivers."""
    rng = np.random.default_rng(seed)
    labels = ["sender"] * 30 + ["receiver"] * 30 + ["bystander"] * 30
    expression = rng.gamma(2.0, 0.5, size=(90, len(GENES)))  # a low background in every gene
    expression[:30, GENES.index("WNT5A")] += 10.0
    expression[30:60, GENES.index("FZD5")] += 10.0
    return expression, labels


def _key(row: dict[str, object]) -> tuple[object, ...]:
    return row["ligand"], row["receptor"], row["source_type"], row["target_type"]


class TestPlantedPair:
    def test_the_planted_pair_is_the_strongest_and_significant(self) -> None:
        expression, labels = _planted_expression()

        result = compute_ligand_receptor_interactions(expression, labels, DATABASE, gene_names=GENES, seed=7)

        top = max(result["interactions"], key=lambda row: row["score"])
        assert _key(top) == PLANTED
        assert top["p_value"] < 0.05
        assert top["score"] > 5.0  # about (11 x 11) / (4.3 x 4.3) = 6.4 for the planted direction
        weaker = [row["score"] for row in result["interactions"] if _key(row) != PLANTED]
        assert max(weaker) < 2.0

    def test_the_decoy_pair_is_not_reported_as_a_strong_interaction(self) -> None:
        expression, labels = _planted_expression()

        result = compute_ligand_receptor_interactions(expression, labels, DATABASE, gene_names=GENES, seed=7)

        decoy = [row for row in result["interactions"] if row["ligand"] == "TGFB1"]
        assert len(decoy) == 9  # three source types x three target types
        assert max(row["score"] for row in decoy) < 2.0

    def test_scores_come_from_the_columns_the_names_point_to(self) -> None:
        expression, labels = _planted_expression()
        types = np.array(labels)

        result = compute_ligand_receptor_interactions(expression, labels, DATABASE, gene_names=GENES, seed=7)

        assert len(result["interactions"]) == 2 * 3 * 3
        for row in result["interactions"]:
            ligand = expression[:, GENES.index(row["ligand"])]
            receptor = expression[:, GENES.index(row["receptor"])]
            expected = (ligand[types == row["source_type"]].mean() * receptor[types == row["target_type"]].mean()) / (
                ligand.mean() * receptor.mean()
            )
            assert row["score"] == pytest.approx(expected)

    def test_the_column_order_does_not_matter(self) -> None:
        expression, labels = _planted_expression()
        order = np.random.default_rng(3).permutation(len(GENES))
        shuffled_names = [GENES[i] for i in order]

        original = compute_ligand_receptor_interactions(expression, labels, DATABASE, gene_names=GENES, seed=7)
        shuffled = compute_ligand_receptor_interactions(
            expression[:, order], labels, DATABASE, gene_names=shuffled_names, seed=7
        )

        assert [_key(row) for row in original["interactions"]] == [_key(row) for row in shuffled["interactions"]]
        np.testing.assert_allclose(
            [row["score"] for row in original["interactions"]], [row["score"] for row in shuffled["interactions"]]
        )
        assert [row["p_value"] for row in original["interactions"]] == [
            row["p_value"] for row in shuffled["interactions"]
        ]

    def test_a_numeric_symbol_is_a_name_and_not_a_column_index(self) -> None:
        # The column named "1" comes first and the column named "0" second. By name the ligand "0" is the
        # SECOND column; the old reading of a numeric string as an index would have taken the first.
        expression = np.array([[10.0, 1.0], [10.0, 1.0], [1.0, 10.0], [1.0, 10.0]])
        labels = ["A", "A", "B", "B"]
        database = {"pairs": [{"ligand": "0", "receptor": "1"}]}

        result = compute_ligand_receptor_interactions(expression, labels, database, gene_names=["1", "0"], seed=1)

        scores = {(row["source_type"], row["target_type"]): row["score"] for row in result["interactions"]}
        assert scores[("B", "A")] == pytest.approx((10.0 * 10.0) / (5.5 * 5.5))  # ligand high in B, receptor high in A
        assert scores[("A", "B")] == pytest.approx((1.0 * 1.0) / (5.5 * 5.5))


class TestGeneNamesAreRequired:
    def test_a_plain_matrix_without_names_raises(self) -> None:
        expression, labels = _planted_expression()

        with pytest.raises(ValueError, match="gene_names is required"):
            compute_ligand_receptor_interactions(expression, labels, DATABASE)

    def test_the_built_in_database_without_names_raises_too(self) -> None:
        expression, labels = _planted_expression()

        with pytest.raises(ValueError, match="gene_names is required"):
            compute_ligand_receptor_interactions(expression, labels)

    def test_fewer_names_than_columns_raise(self) -> None:
        expression, labels = _planted_expression()

        with pytest.raises(ValueError, match="gene_names has 7 entries but the expression matrix has 8 columns"):
            compute_ligand_receptor_interactions(expression, labels, DATABASE, gene_names=GENES[:-1])

    def test_more_names_than_columns_raise(self) -> None:
        expression, labels = _planted_expression()

        with pytest.raises(ValueError, match="gene_names has 9 entries but the expression matrix has 8 columns"):
            compute_ligand_receptor_interactions(expression, labels, DATABASE, gene_names=[*GENES, "EXTRA"])

    def test_a_label_for_every_cell_is_required(self) -> None:
        expression, labels = _planted_expression()

        with pytest.raises(ValueError, match="cell_types has 89 entries but the expression matrix has 90 cells"):
            compute_ligand_receptor_interactions(expression, labels[:-1], DATABASE, gene_names=GENES)

    def test_a_vector_is_not_an_expression_matrix(self) -> None:
        with pytest.raises(ValueError, match="2-D matrix"):
            compute_ligand_receptor_interactions(np.ones(8), ["A"] * 8, DATABASE, gene_names=GENES)


class TestNoPairIsPresent:
    def test_no_database_gene_in_the_data_raises(self) -> None:
        expression, labels = _planted_expression()
        unrelated = [f"gene{i}" for i in range(len(GENES))]

        with pytest.raises(ValueError, match="None of the 2 ligand-receptor pairs has both of its genes"):
            compute_ligand_receptor_interactions(expression, labels, DATABASE, gene_names=unrelated)

    def test_a_pair_with_only_one_of_its_genes_present_is_not_enough(self) -> None:
        expression, labels = _planted_expression()
        names = [("FZD_RENAMED" if name == "FZD5" else name) for name in GENES]
        only_planted = {"pairs": [DATABASE["pairs"][0]]}

        with pytest.raises(ValueError, match="None of the 1 ligand-receptor pairs"):
            compute_ligand_receptor_interactions(expression, labels, only_planted, gene_names=names)

    def test_matching_is_exact_and_case_sensitive(self) -> None:
        expression, labels = _planted_expression()
        lower_case = {"pairs": [{"ligand": "wnt5a", "receptor": "fzd5"}]}

        with pytest.raises(ValueError, match="None of the 1 ligand-receptor pairs"):
            compute_ligand_receptor_interactions(expression, labels, lower_case, gene_names=GENES)

    def test_an_empty_database_raises(self) -> None:
        expression, labels = _planted_expression()

        with pytest.raises(ValueError, match="None of the 0 ligand-receptor pairs"):
            compute_ligand_receptor_interactions(expression, labels, {"pairs": []}, gene_names=GENES)

    def test_the_human_database_does_not_match_mouse_style_names(self) -> None:
        expression, labels = _planted_expression()
        mouse_style = ["Actb", "Tgfb1", "Gapdh", "Wnt5a", "Tgfbr1", "Fzd5", "B2m", "Malat1"]

        with pytest.raises(ValueError, match="same gene-name convention"):
            compute_ligand_receptor_interactions(expression, labels, gene_names=mouse_style)


class TestPartialOverlap:
    def test_only_pairs_with_both_genes_present_are_scored_and_the_rest_are_counted(self) -> None:
        expression, labels = _planted_expression()
        database = {
            "pairs": [
                {"ligand": "WNT5A", "receptor": "FZD5"},
                {"ligand": "TGFB1", "receptor": "TGFBR1"},
                {"ligand": "WNT5A", "receptor": "NOT_IN_THE_DATA"},
                {"ligand": "NOT_EITHER", "receptor": "FZD5"},
            ]
        }

        result = compute_ligand_receptor_interactions(expression, labels, database, gene_names=GENES, seed=1)

        assert {(row["ligand"], row["receptor"]) for row in result["interactions"]} == {
            ("WNT5A", "FZD5"),
            ("TGFB1", "TGFBR1"),
        }
        summary = result["summary"]
        assert summary["n_pairs_in_database"] == 4
        assert summary["n_pairs_scored"] == 2
        assert summary["n_pairs_skipped"] == 2
        assert summary["total_pairs_tested"] == 2 * 3 * 3
        assert summary["unique_ligands"] == 2
        assert summary["unique_receptors"] == 2
        assert summary["unique_cell_types"] == 3

    def test_the_built_in_database_scores_only_what_the_data_contain(self) -> None:
        expression, labels = _planted_expression()
        names = ["CXCL12", "CXCR4", "WNT5A", "FZD5", "EGF", "EGFR", "ACTB", "GAPDH"]

        result = compute_ligand_receptor_interactions(expression, labels, gene_names=names, seed=1)

        present = {(pair["ligand"], pair["receptor"]) for pair in default_lr_database()["pairs"]} & {
            (ligand, receptor) for ligand in names for receptor in names
        }
        assert {(row["ligand"], row["receptor"]) for row in result["interactions"]} == present
        assert ("CXCL12", "CXCR4") in present and ("EGF", "EGFR") in present
        assert result["summary"]["n_pairs_scored"] == len(present)
        assert result["summary"]["n_pairs_in_database"] == len(default_lr_database()["pairs"])


class TestInputKinds:
    def test_a_dataset_supplies_its_own_gene_names(self) -> None:
        expression, labels = _planted_expression()
        dataset = SpatialDataset(
            expression=expression,
            coordinates=np.zeros((90, 2)),
            barcodes=[f"cell{i}" for i in range(90)],
            gene_names=list(GENES),
        )

        from_dataset = compute_ligand_receptor_interactions(dataset, labels, DATABASE, seed=7)
        from_matrix = compute_ligand_receptor_interactions(expression, labels, DATABASE, gene_names=GENES, seed=7)

        assert from_dataset["interactions"] == from_matrix["interactions"]

    def test_explicit_names_win_over_the_names_a_dataset_carries(self) -> None:
        expression, labels = _planted_expression()
        dataset = SpatialDataset(
            expression=expression,
            coordinates=np.zeros((90, 2)),
            barcodes=[f"cell{i}" for i in range(90)],
            gene_names=[f"wrong{i}" for i in range(len(GENES))],
        )

        result = compute_ligand_receptor_interactions(dataset, labels, DATABASE, gene_names=GENES, seed=7)

        assert max(result["interactions"], key=lambda row: row["score"])["ligand"] == "WNT5A"

    def test_a_dataframe_supplies_its_columns_as_names(self) -> None:
        pandas = pytest.importorskip("pandas")
        expression, labels = _planted_expression()
        frame = pandas.DataFrame(expression, columns=GENES)

        from_frame = compute_ligand_receptor_interactions(frame, labels, DATABASE, seed=7)
        from_matrix = compute_ligand_receptor_interactions(expression, labels, DATABASE, gene_names=GENES, seed=7)

        assert from_frame["interactions"] == from_matrix["interactions"]

    def test_sparse_and_list_inputs_give_the_dense_result(self) -> None:
        expression, labels = _planted_expression()
        dense = compute_ligand_receptor_interactions(expression, labels, DATABASE, gene_names=GENES, seed=7)

        from_sparse = compute_ligand_receptor_interactions(
            sparse.csr_matrix(expression), labels, DATABASE, gene_names=GENES, seed=7
        )
        from_lists = compute_ligand_receptor_interactions(
            expression.tolist(), labels, DATABASE, gene_names=GENES, seed=7
        )

        assert from_sparse["interactions"] == dense["interactions"]
        assert from_lists["interactions"] == dense["interactions"]


class TestDuplicateNames:
    def test_a_scored_gene_that_labels_two_columns_is_ambiguous(self) -> None:
        expression, labels = _planted_expression()
        names = list(GENES)
        names[0] = "WNT5A"  # now both column 0 and column 3 are called WNT5A

        with pytest.raises(ValueError, match="more than one expression column"):
            compute_ligand_receptor_interactions(expression, labels, DATABASE, gene_names=names)

    def test_a_repeated_name_that_is_not_scored_does_no_harm(self) -> None:
        expression, labels = _planted_expression()
        names = list(GENES)
        names[0] = "B2M"  # B2M now labels columns 0 and 6, but no database pair uses it

        result = compute_ligand_receptor_interactions(expression, labels, DATABASE, gene_names=names, seed=7)
        baseline = compute_ligand_receptor_interactions(expression, labels, DATABASE, gene_names=GENES, seed=7)

        assert result["interactions"] == baseline["interactions"]
