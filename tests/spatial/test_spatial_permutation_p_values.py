"""Permutation p-values are (b + 1) / (N + 1), never 0, and a z-score without spread is NaN.

``neighborhood_enrichment`` and ``compute_ligand_receptor_interactions`` reported ``b / N``: when no permutation
reached the observed value the p-value was exactly 0, claiming that the value is impossible under the null although
only N permutations were tried. ``neighborhood_enrichment`` also divided by a made-up 1 wherever the null had no
variance. Real data and real calls only.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.spatial.analysis.neighborhood import neighborhood_enrichment
from metainformant.spatial.communication.cell_communication import compute_ligand_receptor_interactions

N_PERMUTATIONS = 50


def _two_halves() -> tuple[np.ndarray, np.ndarray]:
    """Type A on the left half of a line of 40 cells and type B on the right: perfectly segregated."""
    coordinates = np.column_stack([np.arange(40, dtype=float), np.zeros(40)])
    cell_types = np.array(["A"] * 20 + ["B"] * 20)
    return cell_types, coordinates


class TestNeighborhoodEnrichmentPValues:
    def test_a_count_no_permutation_reaches_has_the_smallest_possible_p_not_zero(self) -> None:
        cell_types, coordinates = _two_halves()

        result = neighborhood_enrichment(cell_types, coordinates, n_neighbors=4, n_permutations=N_PERMUTATIONS, seed=1)

        smallest = 1.0 / (N_PERMUTATIONS + 1)
        names = result.cell_type_names
        a, b = names.index("A"), names.index("B")
        # Same-type neighbors are far more common than in any shuffled labelling
        assert result.enrichment_matrix[a, a] > 3.0 and result.enrichment_matrix[b, b] > 3.0
        assert result.p_values[a, a] == pytest.approx(smallest, rel=1e-12)
        assert result.p_values[b, b] == pytest.approx(smallest, rel=1e-12)
        assert (result.p_values > 0.0).all()

    def test_a_depleted_pair_has_p_one_by_the_same_formula(self) -> None:
        cell_types, coordinates = _two_halves()

        result = neighborhood_enrichment(cell_types, coordinates, n_neighbors=4, n_permutations=N_PERMUTATIONS, seed=1)

        a, b = result.cell_type_names.index("A"), result.cell_type_names.index("B")
        assert result.enrichment_matrix[a, b] < -3.0
        assert result.p_values[a, b] == pytest.approx(1.0)  # every permutation has at least as many A-B neighbors

    def test_the_p_value_counts_permutations_at_least_as_extreme(self) -> None:
        """Recompute (b + 1) / (N + 1) from an independent permutation loop with the same seed."""
        cell_types, coordinates = _two_halves()
        rng = np.random.RandomState(7)
        codes = np.array([0] * 20 + [1] * 20)
        neighbors = {i: [j for j in (i - 2, i - 1, i + 1, i + 2) if 0 <= j < 40] for i in range(40)}

        def counts(labels: np.ndarray) -> np.ndarray:
            matrix = np.zeros((2, 2))
            for i, listed in neighbors.items():
                for j in listed:
                    matrix[labels[i], labels[j]] += 1
            return matrix

        result = neighborhood_enrichment(cell_types, coordinates, radius=2.0, n_permutations=30, seed=7)

        observed = counts(codes)
        null = np.array([counts(rng.permutation(codes)) for _ in range(30)])
        expected = ((null >= observed).sum(axis=0) + 1.0) / 31.0
        np.testing.assert_allclose(result.count_matrix, observed)
        np.testing.assert_allclose(result.p_values, expected, rtol=1e-12)

    def test_a_null_without_spread_has_a_nan_enrichment_not_a_zero(self) -> None:
        coordinates = np.column_stack([np.arange(12, dtype=float), np.zeros(12)])

        result = neighborhood_enrichment(["A"] * 12, coordinates, n_neighbors=3, n_permutations=10, seed=0)

        assert np.isnan(result.enrichment_matrix).all()  # one cell type: every shuffle gives the same counts
        assert result.p_values.shape == (1, 1) and result.p_values[0, 0] == pytest.approx(1.0)

    def test_enrichment_is_finite_where_the_null_has_spread(self) -> None:
        cell_types, coordinates = _two_halves()

        result = neighborhood_enrichment(cell_types, coordinates, n_neighbors=4, n_permutations=N_PERMUTATIONS, seed=1)

        assert np.isfinite(result.enrichment_matrix).all()

    @pytest.mark.parametrize("n_permutations", [0, -3])
    def test_fewer_than_one_permutation_is_refused(self, n_permutations: int) -> None:
        cell_types, coordinates = _two_halves()

        with pytest.raises(ValueError, match="n_permutations"):
            neighborhood_enrichment(cell_types, coordinates, n_permutations=n_permutations)


class TestLigandReceptorPValues:
    @staticmethod
    def _expression() -> tuple[np.ndarray, list[str], list[str]]:
        """Ligand L only in type A and receptor R only in type B: A -> B is far above any shuffled labelling."""
        rng = np.random.default_rng(0)
        expression = rng.uniform(0.0, 0.1, (40, 3))
        expression[:20, 0] += 5.0  # L in the A cells
        expression[20:, 1] += 5.0  # R in the B cells
        return expression, ["A"] * 20 + ["B"] * 20, ["L", "R", "other"]

    def test_an_interaction_no_permutation_reaches_has_p_one_over_101_not_zero(self) -> None:
        expression, cell_types, gene_names = self._expression()

        result = compute_ligand_receptor_interactions(
            expression,
            cell_types,
            lr_database={"pairs": [{"ligand": "L", "receptor": "R"}]},
            gene_names=gene_names,
            seed=3,
        )

        by_pair = {(ix["source_type"], ix["target_type"]): ix for ix in result["interactions"]}
        assert by_pair[("A", "B")]["score"] > by_pair[("B", "A")]["score"]
        assert by_pair[("A", "B")]["p_value"] == pytest.approx(1.0 / 101.0, rel=1e-12)
        assert min(ix["p_value"] for ix in result["interactions"]) > 0.0
        assert result["n_significant"] >= 1

    def test_the_p_value_is_b_plus_one_over_n_plus_one(self) -> None:
        expression, cell_types, gene_names = self._expression()
        result = compute_ligand_receptor_interactions(
            expression,
            cell_types,
            lr_database={"pairs": [{"ligand": "L", "receptor": "R"}]},
            gene_names=gene_names,
            seed=3,
        )

        for ix in result["interactions"]:
            numerator = ix["p_value"] * 101.0
            assert numerator == pytest.approx(round(numerator), abs=1e-9)  # an integer count b + 1
            assert 1 <= round(numerator) <= 101
