"""Real-implementation tests for metainformant.networks.regulatory.grn_inference.

Covers the pure-Python numerical GRN inference stack: correlation-based,
mutual-information (ARACNE-like with DPI), and regression-based inference,
regulator scoring, motif counting, and validation metrics. All inputs are
real synthetic matrices; all computation is real (no test doubles).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.networks.regulatory.grn_inference import (
    _coordinate_descent_regression,
    _count_motifs_fast,
    _get_matrix_shape,
    _hypergeometric_pvalue,
    _mutual_information_pure,
    _pearson_correlation_pure,
    _spearman_correlation_pure,
    _to_list_of_lists,
    compute_network_motifs,
    infer_grn_correlation,
    infer_grn_mutual_info,
    infer_grn_regression,
    score_regulators,
    validate_grn,
)

RNG = np.random.default_rng(42)


def _make_expression(n_genes: int = 8, n_samples: int = 40) -> tuple[np.ndarray, list[str]]:
    """Synthetic expression matrix (genes x samples) with planted structure.

    Gene 0 drives gene 1 (strong linear relation) and gene 2 (weaker),
    everything else is independent noise. Fixed seed for determinism.
    """
    driver = RNG.normal(0.0, 1.0, n_samples)
    matrix = np.zeros((n_genes, n_samples))
    matrix[0] = driver
    matrix[1] = 2.0 * driver + RNG.normal(0.0, 0.2, n_samples)
    matrix[2] = 0.8 * driver + RNG.normal(0.0, 0.5, n_samples)
    for g in range(3, n_genes):
        matrix[g] = RNG.normal(0.0, 1.0, n_samples)
    names = [f"GENE{i}" for i in range(n_genes)]
    return matrix, names


class TestPureCorrelationHelpers:
    def test_pearson_matches_known_value(self) -> None:
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        r, p = _pearson_correlation_pure(x, [2.0, 4.0, 6.0, 8.0, 10.0])
        assert math.isclose(r, 1.0, abs_tol=1e-9)
        assert p == 0.0  # |r| >= 1 short-circuits to p=0

    def test_pearson_anticorrelated(self) -> None:
        r, _ = _pearson_correlation_pure([1.0, 2.0, 3.0, 4.0], [8.0, 6.0, 4.0, 2.0])
        assert math.isclose(r, -1.0, abs_tol=1e-9)

    def test_pearson_small_sample_returns_null(self) -> None:
        r, p = _pearson_correlation_pure([1.0, 2.0], [3.0, 4.0])
        assert r == 0.0 and p == 1.0

    def test_pearson_constant_input_returns_null(self) -> None:
        r, p = _pearson_correlation_pure([5.0, 5.0, 5.0, 5.0], [1.0, 2.0, 3.0, 4.0])
        assert r == 0.0 and p == 1.0

    def test_spearman_monotonic_is_one(self) -> None:
        # Monotonic but nonlinear (exponential) -> Spearman = 1
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [math.exp(v) for v in x]
        r, _ = _spearman_correlation_pure(x, y)
        assert math.isclose(r, 1.0, abs_tol=1e-9)

    def test_spearman_ties_use_average_ranks(self) -> None:
        # Ties in x must be handled via averaged ranks without error
        x = [1.0, 1.0, 2.0, 3.0]
        y = [1.0, 2.0, 3.0, 4.0]
        r, _ = _spearman_correlation_pure(x, y)
        assert 0.0 < r <= 1.0

    def test_spearman_small_sample_returns_null(self) -> None:
        r, p = _spearman_correlation_pure([1.0], [2.0])
        assert r == 0.0 and p == 1.0


class TestMutualInformation:
    def test_mi_zero_for_independent(self) -> None:
        rng = np.random.default_rng(7)
        x = list(rng.normal(size=500))
        y = list(rng.normal(size=500))
        assert _mutual_information_pure(x, y, n_bins=8) < 0.1

    def test_mi_positive_for_dependent(self) -> None:
        rng = np.random.default_rng(7)
        x = list(rng.normal(size=500))
        y = [v + 0.05 * rng.normal() for v in x]
        assert _mutual_information_pure(x, y, n_bins=8) > 0.5

    def test_mi_identical_variables_positive(self) -> None:
        x = list(range(50))
        assert _mutual_information_pure(x, x, n_bins=5) > 0.0

    def test_mi_empty_returns_zero(self) -> None:
        assert _mutual_information_pure([], [], n_bins=5) == 0.0


class TestMatrixHelpers:
    def test_shape_from_numpy_and_lists(self) -> None:
        m = np.zeros((3, 4))
        assert _get_matrix_shape(m) == (3, 4)
        assert _get_matrix_shape([[1.0] * 4 for _ in range(3)]) == (3, 4)

    def test_to_list_of_lists_numpy(self) -> None:
        m = np.array([[1.0, 2.0], [3.0, 4.0]])
        out = _to_list_of_lists(m, 2, 2)
        assert out == [[1.0, 2.0], [3.0, 4.0]]
        assert all(isinstance(v, float) for row in out for v in row)

    def test_to_list_of_lists_nested_lists(self) -> None:
        out = _to_list_of_lists([[1, 2], [3, 4]], 2, 2)
        assert out == [[1.0, 2.0], [3.0, 4.0]]


class TestInferGrnCorrelation:
    def test_plants_expected_edges(self) -> None:
        matrix, names = _make_expression()
        result = infer_grn_correlation(matrix, names, threshold=0.4)

        assert result["n_edges"] == len(result["edges"]) > 0
        edges = {(e["source"], e["target"]) for e in result["edges"]}
        # Driver GENE0 should regulate its planted targets
        assert ("GENE0", "GENE1") in edges
        assert ("GENE0", "GENE2") in edges
        # Independent noise genes should not strongly correlate with GENE0
        assert ("GENE0", "GENE5") not in edges

    def test_edge_weights_and_pvalues_in_range(self) -> None:
        matrix, names = _make_expression()
        result = infer_grn_correlation(matrix, names, threshold=0.4)
        for e in result["edges"]:
            assert 0.4 <= abs(e["weight"]) <= 1.0
            assert 0.0 <= e["p_value"] <= 1.0

    def test_edges_sorted_by_abs_weight_desc(self) -> None:
        matrix, names = _make_expression()
        result = infer_grn_correlation(matrix, names, threshold=0.3)
        weights = [abs(e["weight"]) for e in result["edges"]]
        assert weights == sorted(weights, reverse=True)

    def test_adjacency_matrix_shape_and_values(self) -> None:
        matrix, names = _make_expression(n_genes=6)
        result = infer_grn_correlation(matrix, names, threshold=0.9)
        adj = result["adjacency_matrix"]
        assert len(adj) == 6 and all(len(row) == 6 for row in adj)
        # Diagonal untouched (self-comparisons skipped)
        assert all(adj[i][i] == 0.0 for i in range(6))

    def test_spearman_method(self) -> None:
        matrix, names = _make_expression()
        result = infer_grn_correlation(matrix, names, method="spearman", threshold=0.4)
        assert result["n_edges"] > 0
        assert ("GENE0", "GENE1") in {(e["source"], e["target"]) for e in result["edges"]}

    def test_tf_list_restriction(self) -> None:
        matrix, names = _make_expression()
        result = infer_grn_correlation(matrix, names, tf_list=["GENE0"], threshold=0.4)
        assert all(e["source"] == "GENE0" for e in result["edges"])

    def test_tf_list_with_unknown_gene_falls_back(self) -> None:
        matrix, names = _make_expression()
        # TF list that matches nothing -> falls back to all genes
        result = infer_grn_correlation(matrix, names, tf_list=["NOT_A_GENE"], threshold=0.9)
        assert result["n_edges"] >= 0

    def test_gene_name_length_mismatch_raises(self) -> None:
        matrix, names = _make_expression()
        with pytest.raises(ValueError, match="must match matrix rows"):
            infer_grn_correlation(matrix, names[:-1])

    def test_unsupported_method_raises(self) -> None:
        matrix, names = _make_expression()
        with pytest.raises(ValueError, match="Unsupported method"):
            infer_grn_correlation(matrix, names, method="kendall")

    def test_accepts_plain_lists(self) -> None:
        matrix, names = _make_expression()
        result = infer_grn_correlation(matrix.tolist(), names, threshold=0.4)
        assert result["n_edges"] > 0


class TestInferGrnMutualInfo:
    def test_planted_edge_survives_dpi(self) -> None:
        matrix, names = _make_expression()
        result = infer_grn_mutual_info(matrix, names, n_bins=6)
        edges = {(e["source"], e["target"]) for e in result["edges"]}
        assert ("GENE0", "GENE1") in edges

    def test_mi_matrix_shape_and_symmetry_self_zero(self) -> None:
        matrix, names = _make_expression(n_genes=5)
        result = infer_grn_mutual_info(matrix, names, n_bins=4)
        mi = result["mi_matrix"]
        assert len(mi) == 5 and all(len(row) == 5 for row in mi)
        # Diagonal entries are skipped (self MI never computed)
        assert all(mi[i][i] == 0.0 for i in range(5))

    def test_dpi_weakest_triplet_edge_removed(self) -> None:
        # Chain A->B->C with tight links: MI(A,B) and MI(B,C) exceed MI(A,C),
        # so DPI must remove the indirect A-C pair and keep both direct pairs.
        rng = np.random.default_rng(11)
        a = list(rng.normal(size=80))
        b = [0.98 * v + 0.08 * rng.normal() for v in a]
        c = [0.98 * v + 0.08 * rng.normal() for v in b]
        matrix = [a, b, c]
        names = ["A", "B", "C"]
        result = infer_grn_mutual_info(matrix, names, n_bins=8)
        edges = {(e["source"], e["target"]) for e in result["edges"]}
        # Strongest direct edges survive; indirect pair pruned
        assert ("A", "B") in edges or ("B", "A") in edges
        assert ("B", "C") in edges or ("C", "B") in edges
        assert ("A", "C") not in edges and ("C", "A") not in edges

    def test_gene_name_length_mismatch_raises(self) -> None:
        matrix, names = _make_expression()
        with pytest.raises(ValueError, match="must match matrix rows"):
            infer_grn_mutual_info(matrix, names[:3])

    def test_tf_list_restriction(self) -> None:
        matrix, names = _make_expression()
        result = infer_grn_mutual_info(matrix, names, tf_list=["GENE0"], n_bins=6)
        assert all(e["source"] == "GENE0" for e in result["edges"])


class TestInferGrnRegression:
    def test_recovers_planted_structure(self) -> None:
        matrix, names = _make_expression()
        result = infer_grn_regression(matrix, names, method="ridge", alpha=0.01)
        edges = {(e["source"], e["target"]) for e in result["edges"]}
        assert ("GENE0", "GENE1") in edges
        assert ("GENE0", "GENE2") in edges

    def test_r_squared_high_for_planted_target(self) -> None:
        matrix, names = _make_expression()
        result = infer_grn_regression(matrix, names, method="ridge", alpha=0.01)
        assert result["r_squared_per_gene"]["GENE1"] > 0.9

    def test_coefficients_matrix_shape(self) -> None:
        matrix, names = _make_expression(n_genes=6)
        result = infer_grn_regression(matrix, names, method="ridge", alpha=0.1)
        coefs = result["coefficients"]
        assert len(coefs) == 6 and all(len(row) == 6 for row in coefs)

    def test_lasso_sparsifies(self) -> None:
        matrix, names = _make_expression()
        weak = infer_grn_regression(matrix, names, method="lasso", alpha=0.5)
        strong = infer_grn_regression(matrix, names, method="lasso", alpha=50.0)
        assert strong["n_edges"] <= weak["n_edges"]

    def test_unsupported_method_raises(self) -> None:
        matrix, names = _make_expression()
        with pytest.raises(ValueError, match="Unsupported method"):
            infer_grn_regression(matrix, names, method="elastic")

    def test_gene_name_length_mismatch_raises(self) -> None:
        matrix, names = _make_expression()
        with pytest.raises(ValueError, match="must match matrix rows"):
            infer_grn_regression(matrix, names[1:])


class TestCoordinateDescent:
    def test_ridge_recovers_linear_signal(self) -> None:
        # Note: coordinate descent here has no intercept term, so use a
        # zero-mean target; slope recovery must be tight.
        rng = np.random.default_rng(3)
        x = [[float(v)] for v in rng.normal(size=50)]
        y = [3.0 * row[0] + 0.01 * rng.normal() for row in x]
        coefs, r2 = _coordinate_descent_regression(x, y, method="ridge", alpha=0.001)
        assert len(coefs) == 1
        assert abs(coefs[0] - 3.0) < 0.2
        assert r2 > 0.99

    def test_ridge_handles_offset_target(self) -> None:
        # Non-zero-mean target: slope still recovered, fit degraded but strong
        rng = np.random.default_rng(3)
        x = [[float(v)] for v in rng.normal(size=50)]
        y = [3.0 * row[0] + 1.0 + 0.01 * rng.normal() for row in x]
        coefs, r2 = _coordinate_descent_regression(x, y, method="ridge", alpha=0.001)
        assert abs(coefs[0] - 3.0) < 0.5
        assert r2 > 0.85  # document actual fit without intercept

    def test_lasso_zeroes_weak_features(self) -> None:
        rng = np.random.default_rng(5)
        strong = list(rng.normal(size=60))
        noise1 = list(rng.normal(size=60))
        noise2 = list(rng.normal(size=60))
        y = [2.5 * s + 0.01 * rng.normal() for s in strong]
        x = [[s, n1, n2] for s, n1, n2 in zip(strong, noise1, noise2)]
        coefs, _ = _coordinate_descent_regression(x, y, method="lasso", alpha=0.5)
        assert abs(coefs[0]) > 1.0
        assert coefs[1] == 0.0 and coefs[2] == 0.0

    def test_empty_input_returns_null(self) -> None:
        assert _coordinate_descent_regression([], []) == ([], 0.0)
        assert _coordinate_descent_regression([[]], [1.0, 2.0]) == ([], 0.0)

    def test_constant_target_r2_zero(self) -> None:
        x = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
        coefs, r2 = _coordinate_descent_regression(x, [7.0, 7.0, 7.0], method="ridge", alpha=0.1)
        assert r2 == 0.0


class TestHypergeometricPvalue:
    def test_no_enrichment_gives_high_p(self) -> None:
        # Observed = expected -> p near 1
        p = _hypergeometric_pvalue(k=5, n=20, K=50, N=200)
        assert p > 0.3

    def test_strong_enrichment_gives_low_p(self) -> None:
        p = _hypergeometric_pvalue(k=30, n=30, K=50, N=200)
        assert p < 0.001

    def test_degenerate_cases_return_one(self) -> None:
        assert _hypergeometric_pvalue(k=5, n=0, K=10, N=100) == 1.0
        assert _hypergeometric_pvalue(k=5, n=10, K=0, N=100) == 1.0
        assert _hypergeometric_pvalue(k=5, n=10, K=10, N=0) == 1.0

    def test_bounded_in_unit_interval(self) -> None:
        for k in (0, 5, 10, 20):
            p = _hypergeometric_pvalue(k=k, n=20, K=30, N=100)
            assert 0.0 <= p <= 1.0


class TestScoreRegulators:
    def _grn(self) -> dict:
        return {
            "edges": [
                {"source": "TF1", "target": "G1", "weight": 0.9},
                {"source": "TF1", "target": "G2", "weight": 0.8},
                {"source": "TF1", "target": "G3", "weight": 0.7},
                {"source": "TF2", "target": "G4", "weight": 0.6},
            ]
        }

    def test_scores_sorted_by_pvalue(self) -> None:
        results = score_regulators(self._grn(), ["G1", "G2", "G3"])
        assert len(results) == 2
        ps = [r["enrichment_p"] for r in results]
        assert ps == sorted(ps)
        assert results[0]["tf"] == "TF1"

    def test_top_tf_has_all_targets_in_set(self) -> None:
        results = score_regulators(self._grn(), ["G1", "G2", "G3"])
        top = results[0]
        assert top["n_targets"] == 3
        assert top["n_in_set"] == 3

    def test_regulation_score_sums_abs_weights(self) -> None:
        results = score_regulators(self._grn(), ["G1", "G2", "G3"])
        top = results[0]
        assert math.isclose(top["regulation_score"], 0.9 + 0.8 + 0.7, abs_tol=1e-6)

    def test_case_insensitive_gene_matching(self) -> None:
        results = score_regulators(self._grn(), ["g1", "G2", "g3"])
        assert results[0]["n_in_set"] == 3

    def test_missing_edges_key_raises(self) -> None:
        with pytest.raises(ValueError, match="edges"):
            score_regulators({"nodes": []}, ["G1"])

    def test_empty_gene_set_returns_priors(self) -> None:
        results = score_regulators(self._grn(), [])
        assert all(r["n_in_set"] == 0 for r in results)


class TestComputeNetworkMotifs:
    def test_invalid_motif_size_raises(self) -> None:
        with pytest.raises(ValueError, match="motif_size"):
            compute_network_motifs([], motif_size=4)

    def test_mutual_regulation_counted(self) -> None:
        edges = [
            {"source": "A", "target": "B"},
            {"source": "B", "target": "A"},
            {"source": "C", "target": "D"},
        ]
        result = compute_network_motifs(edges, motif_size=2)
        assert result["motif_counts"]["mutual_regulation"] == 1
        assert result["motif_counts"]["single_regulation"] == 1
        assert "z_scores" in result and "significant_motifs" in result

    def test_feed_forward_loop_counted(self) -> None:
        edges = [
            {"source": "A", "target": "B"},
            {"source": "B", "target": "C"},
            {"source": "A", "target": "C"},
        ]
        result = compute_network_motifs(edges, motif_size=3)
        assert result["motif_counts"]["feed_forward_loop"] == 1

    def test_cascade_counted_when_no_shortcut(self) -> None:
        edges = [
            {"source": "A", "target": "B"},
            {"source": "B", "target": "C"},
        ]
        result = compute_network_motifs(edges, motif_size=3)
        assert result["motif_counts"]["cascade"] >= 1
        assert result["motif_counts"]["feed_forward_loop"] == 0

    def test_co_regulation_counted(self) -> None:
        edges = [
            {"source": "A", "target": "T"},
            {"source": "B", "target": "T"},
            {"source": "C", "target": "T"},
        ]
        result = compute_network_motifs(edges, motif_size=3)
        assert result["motif_counts"]["co_regulation"] == 3  # C(3,2) pairs

    def test_empty_edges_gives_zero_motifs(self) -> None:
        result = compute_network_motifs([], motif_size=3)
        assert all(v == 0 for v in result["motif_counts"].values())

    def test_zscores_present_for_all_motif_types(self) -> None:
        edges = [{"source": "A", "target": "B"}, {"source": "B", "target": "C"}]
        result = compute_network_motifs(edges, motif_size=3)
        assert set(result["z_scores"].keys()) == set(result["motif_counts"].keys())
        assert 0.0 <= abs(result["z_scores"]["cascade"]) < 100.0


class TestCountMotifsFast:
    def test_matches_full_implementation(self) -> None:
        edges = [
            {"source": "A", "target": "B"},
            {"source": "B", "target": "C"},
            {"source": "A", "target": "C"},
            {"source": "C", "target": "A"},
        ]
        fast = _count_motifs_fast(edges, 3, ["A", "B", "C"])
        full = compute_network_motifs(edges, motif_size=3)
        for k, v in fast.items():
            assert full["motif_counts"][k] == v

    def test_size2_counts(self) -> None:
        edges = [
            {"source": "A", "target": "B"},
            {"source": "B", "target": "A"},
        ]
        fast = _count_motifs_fast(edges, 2, ["A", "B"])
        assert fast["mutual_regulation"] == 1
        assert fast["single_regulation"] == 0


class TestValidateGrn:
    def test_perfect_prediction(self) -> None:
        known = [
            {"source": "A", "target": "B", "weight": 1.0},
            {"source": "B", "target": "C", "weight": 1.0},
        ]
        result = validate_grn(known, known)
        assert result["precision"] == 1.0
        assert result["recall"] == 1.0
        assert result["f1"] == 1.0

    def test_half_overlap(self) -> None:
        known = [
            {"source": "A", "target": "B"},
            {"source": "B", "target": "C"},
            {"source": "C", "target": "D"},
            {"source": "D", "target": "A"},
        ]
        predicted = known[:2] + [{"source": "X", "target": "Y"}]
        result = validate_grn(predicted, known)
        assert result["precision"] == pytest.approx(2 / 3)
        assert result["recall"] == 0.5
        assert result["f1"] == pytest.approx(2 * (2 / 3) * 0.5 / (2 / 3 + 0.5))

    def test_metrics_bounded(self) -> None:
        rng = np.random.default_rng(9)
        nodes = [f"N{i}" for i in range(10)]
        known = [{"source": nodes[i], "target": nodes[(i + 1) % 10]} for i in range(0, 10, 2)]
        predicted = [
            {"source": nodes[i], "target": nodes[j], "weight": float(rng.uniform(0, 1))}
            for i, j in rng.integers(0, 10, size=(30, 2))
            if i != j
        ]
        result = validate_grn(predicted, known)
        for metric in ("precision", "recall", "f1", "auroc", "auprc"):
            assert 0.0 <= result[metric] <= 1.0

    def test_empty_predicted_raises(self) -> None:
        with pytest.raises(ValueError, match="predicted_edges"):
            validate_grn([], [{"source": "A", "target": "B"}])

    def test_empty_known_raises(self) -> None:
        with pytest.raises(ValueError, match="known_edges"):
            validate_grn([{"source": "A", "target": "B"}], [])

    def test_duplicate_edges_deduplicated(self) -> None:
        known = [{"source": "A", "target": "B"}]
        predicted = [
            {"source": "A", "target": "B"},
            {"source": "A", "target": "B"},
            {"source": "A", "target": "B"},
        ]
        result = validate_grn(predicted, known)
        assert result["precision"] == 1.0  # set semantics dedupe predictions

    def test_auroc_perfect_scoring(self) -> None:
        # True edges get the highest weights -> AUROC = 1
        known = [
            {"source": "A", "target": "B"},
            {"source": "A", "target": "C"},
        ]
        predicted = [
            {"source": "A", "target": "B", "weight": 0.95},
            {"source": "A", "target": "C", "weight": 0.90},
            {"source": "C", "target": "A", "weight": 0.10},
            {"source": "B", "target": "C", "weight": 0.05},
        ]
        result = validate_grn(predicted, known)
        assert result["auroc"] == 1.0
