"""Real-implementation tests for metainformant.rna.deconvolution.bulk_deconvolution.

Depth coverage: NNLS deconvolution (scipy + pure-Python projected gradient),
SVR (CIBERSORT-style) with NNLS fallback, signature matrix construction, marker
gene selection, validation metrics, and batch processing. All inputs are real
synthetic numeric data; no test doubles. Descriptive-only statistics.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.rna.deconvolution.bulk_deconvolution import (
    _compute_r_squared,
    _compute_rmse,
    _get_expression_list,
    _get_mean_expression,
    _matrix_vector_multiply,
    _nnls_projected_gradient_numpy,
    _nnls_projected_gradient_pure,
    _pearson_corr_numpy,
    _simple_t_test,
    _to_flat_list,
    _to_nested_list,
    _validate_and_convert,
    batch_deconvolve,
    build_signature_matrix,
    deconvolve_nnls,
    deconvolve_svr,
    select_marker_genes,
    validate_deconvolution,
)

# Well-conditioned 2-cell-type system: 6 genes x 2 types.
SIGNATURE = [
    [8.0, 2.0],
    [2.5, 1.0],
    [7.0, 3.0],
    [1.0, 0.5],
    [5.0, 4.0],
    [6.5, 1.5],
]
TRUE_PROPS = [0.6, 0.4]
MIXTURE = [sum(s * p for s, p in zip(row, TRUE_PROPS)) for row in SIGNATURE]


class TestDeconvolveNnls:
    def test_recovers_known_proportions(self) -> None:
        result = deconvolve_nnls(MIXTURE, SIGNATURE)
        props = result["proportions"]
        assert len(props) == 2
        assert props[0] == pytest.approx(0.6, abs=1e-3)
        assert props[1] == pytest.approx(0.4, abs=1e-3)
        assert result["r_squared"] > 0.999
        assert result["rmse"] < 1e-3
        assert result["method"] in ("scipy_nnls", "projected_gradient")
        assert result["n_genes"] == 6 and result["n_cell_types"] == 2

    def test_unnormalized_sums_to_unconstrained_scale(self) -> None:
        # Without normalization, an exact 0.6/0.4 mixture is recovered as
        # 0.6/0.4 (solution norm 1.0) but with a 2x mixture the raw solution
        # scales to 2.0 total.
        result = deconvolve_nnls([2 * v for v in MIXTURE], SIGNATURE, normalize=False)
        total = sum(result["proportions"].values())
        assert total == pytest.approx(2.0, abs=1e-3)

    def test_normalized_proportions_sum_to_one(self) -> None:
        result = deconvolve_nnls([2 * v for v in MIXTURE], SIGNATURE, normalize=True)
        assert sum(result["proportions"].values()) == pytest.approx(1.0, abs=1e-6)

    def test_accepts_numpy_arrays(self) -> None:
        result = deconvolve_nnls(np.array(MIXTURE), np.array(SIGNATURE))
        assert result["proportions"][0] == pytest.approx(0.6, abs=1e-3)

    def test_mixture_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="doesn't match"):
            deconvolve_nnls(MIXTURE[:-1], SIGNATURE)

    def test_non_2d_signature_raises(self) -> None:
        with pytest.raises(ValueError, match="2D"):
            deconvolve_nnls(MIXTURE, [1.0, 2.0, 3.0])

    def test_result_metrics_bounded(self) -> None:
        result = deconvolve_nnls(MIXTURE, SIGNATURE)
        assert result["residual"] >= 0.0
        assert 0.0 <= result["rmse"]
        assert 0.0 <= result["r_squared"] <= 1.0 + 1e-9


class TestPureNnlsSolver:
    def test_pure_solver_matches_scipy(self) -> None:
        x_pure, res_pure = _nnls_projected_gradient_pure(SIGNATURE, MIXTURE)
        # Recovered proportions must be near the planted 0.6/0.4
        total = sum(x_pure)
        scaled = [v / total for v in x_pure]
        assert scaled[0] == pytest.approx(0.6, abs=5e-3)
        assert scaled[1] == pytest.approx(0.4, abs=5e-3)
        assert res_pure < 1e-3

    def test_numpy_projected_gradient_solver(self) -> None:
        A = np.array(SIGNATURE)
        b = np.array(MIXTURE)
        x, residual = _nnls_projected_gradient_numpy(A, b)
        total = x.sum()
        scaled = (x / total).tolist()
        assert scaled[0] == pytest.approx(0.6, abs=5e-3)
        assert (x >= 0).all()
        assert residual < 1e-3

    def test_nonnegativity_constraint(self) -> None:
        # A deliberately over-fit system still returns non-negative solution
        x, _ = _nnls_projected_gradient_pure([[1.0, 0.0], [0.0, 1.0]], [-1.0, 1.0])
        assert all(v >= 0.0 for v in x)


class TestValidationAndHelpers:
    def test_validate_and_convert_lists(self) -> None:
        mix, sig, n_genes, n_types = _validate_and_convert(MIXTURE, SIGNATURE)
        assert n_genes == 6 and n_types == 2
        assert len(mix) == 6

    def test_validate_and_convert_numpy(self) -> None:
        mix, sig, n_genes, n_types = _validate_and_convert(np.array(MIXTURE), np.array(SIGNATURE))
        assert n_genes == 6 and n_types == 2
        assert sig.shape == (6, 2)

    def test_matrix_vector_multiply(self) -> None:
        result = _matrix_vector_multiply([[1.0, 2.0], [3.0, 4.0]], [1.0, 1.0])
        assert result == [3.0, 7.0]

    def test_to_flat_and_nested(self) -> None:
        assert _to_flat_list([1, 2, 3], 3) == [1.0, 2.0, 3.0]
        assert _to_nested_list([[1, 2], [3, 4]], 2, 2) == [[1.0, 2.0], [3.0, 4.0]]

    def test_r_squared_and_rmse(self) -> None:
        observed = [1.0, 2.0, 3.0, 4.0]
        assert _compute_r_squared(observed, observed) == 1.0
        assert _compute_rmse(observed, observed) == 0.0
        predicted = [2.0, 3.0, 4.0, 5.0]  # off by 1 everywhere
        assert _compute_rmse(observed, predicted) == pytest.approx(1.0)
        # Perfect linear fit: R^2 = 1; huge misfit: negative
        assert _compute_r_squared(observed, [1.0, 1.0, 1.0, 1.0]) < 0.0

    def test_r_squared_empty(self) -> None:
        assert _compute_r_squared([], []) == 0.0
        assert _compute_rmse([], []) == 0.0

    def test_pearson_corr(self) -> None:
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        assert _pearson_corr_numpy(np.array(x), np.array(x)) == pytest.approx(1.0)
        assert _pearson_corr_numpy(np.array(x), np.array(x[::-1])) == pytest.approx(-1.0)
        assert _pearson_corr_numpy(np.array([1.0]), np.array([1.0])) == 0.0  # too short

    @pytest.mark.filterwarnings("ignore:invalid value encountered")
    def test_pearson_corr_constant_returns_zero(self) -> None:
        # Constant input -> NaN correlation -> 0.0 (numpy emits a RuntimeWarning)
        assert _pearson_corr_numpy(np.array([2.0, 2.0, 2.0]), np.array([1.0, 2.0, 3.0])) == 0.0


class TestDeconvolveSvr:
    def test_recovers_approximate_proportions(self) -> None:
        result = deconvolve_svr(MIXTURE, SIGNATURE)
        assert result["method"] == "svr_cibersort"
        props = result["proportions"]
        assert props[0] == pytest.approx(0.6, abs=0.2)
        assert sum(props.values()) == pytest.approx(1.0, abs=1e-6)
        assert result["correlation"] > 0.9

    def test_nu_validation(self) -> None:
        for bad_nu in (0.0, -0.5, 1.5):
            with pytest.raises(ValueError, match="nu must be"):
                deconvolve_svr(MIXTURE, SIGNATURE, nu=bad_nu)

    def test_small_gene_set_skips_permutation(self) -> None:
        result = deconvolve_svr(MIXTURE, SIGNATURE)  # 6 genes <= 20
        assert result["p_value"] == 0.0

    def test_dimension_mismatch_raises(self) -> None:
        with pytest.raises(ValueError):
            deconvolve_svr(MIXTURE[:-1], SIGNATURE)


class TestSignatureMatrix:
    def _profiles(self) -> dict:
        return {
            "T_cell": {"CD3": [100, 104], "CD4": [80, 76], "CD19": [2, 3], "HBA1": [1, 1]},
            "B_cell": {"CD3": [5, 4], "CD4": [2, 2], "CD19": [95, 90], "HBA1": [1, 2]},
        }

    def test_build_basic(self) -> None:
        result = build_signature_matrix(self._profiles(), ["T_cell", "B_cell"], n_markers=2)
        assert result["n_cell_types"] == 2
        assert result["cell_types"] == ["T_cell", "B_cell"]
        assert len(result["matrix"]) == result["n_genes"]
        assert all(len(row) == 2 for row in result["matrix"])
        assert set(result["markers_per_type"]) <= {"T_cell", "B_cell"}

    def test_markers_are_type_specific(self) -> None:
        result = build_signature_matrix(self._profiles(), ["T_cell", "B_cell"], n_markers=2)
        assert "CD3" in result["markers_per_type"]["T_cell"]
        assert "CD4" in result["markers_per_type"]["T_cell"]
        assert "CD19" in result["markers_per_type"]["B_cell"]

    def test_empty_profiles_raises(self) -> None:
        with pytest.raises(ValueError, match="empty"):
            build_signature_matrix({}, ["T_cell"])

    def test_unknown_cell_type_raises(self) -> None:
        with pytest.raises(ValueError, match="not found"):
            build_signature_matrix(self._profiles(), ["T_cell", "NK_cell"])


class TestSelectMarkerGenes:
    def _expr(self) -> dict:
        return {
            "type_a": {"g1": 100, "g2": 5, "g3": 50, "shared": 10},
            "type_b": {"g1": 3, "g2": 90, "g3": 45, "shared": 10},
        }

    def test_fold_change_selects_specific_genes(self) -> None:
        markers = select_marker_genes(self._expr(), ["type_a", "type_b"], n_genes=1)
        assert "g1" in markers  # top marker for type_a
        assert "g2" in markers  # top marker for type_b

    def test_specificity_method(self) -> None:
        markers = select_marker_genes(self._expr(), ["type_a", "type_b"], n_genes=1, method="specificity")
        assert "g1" in markers and "g2" in markers

    def test_t_test_method_with_sample_lists(self) -> None:
        expr = {
            "type_a": {"g1": [100.0, 102.0, 98.0], "g2": [1.0, 1.5, 0.5]},
            "type_b": {"g1": [2.0, 3.0, 2.5], "g2": [90.0, 92.0, 88.0]},
        }
        markers = select_marker_genes(expr, ["type_a", "type_b"], n_genes=1, method="t_test")
        assert "g1" in markers and "g2" in markers

    def test_t_test_falls_back_without_replicates(self) -> None:
        # Scalar values -> no t-test possible -> fold-change fallback
        markers = select_marker_genes(self._expr(), ["type_a", "type_b"], n_genes=1, method="t_test")
        assert len(markers) > 0

    def test_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown method"):
            select_marker_genes(self._expr(), ["type_a", "type_b"], method="auc")

    def test_unknown_label_raises(self) -> None:
        with pytest.raises(ValueError, match="not found"):
            select_marker_genes(self._expr(), ["type_z"])

    def test_empty_expression_raises(self) -> None:
        with pytest.raises(ValueError, match="empty"):
            select_marker_genes({}, ["type_a"])

    def test_dedup_preserves_order(self) -> None:
        # Same top gene for both types is deduplicated
        expr = {
            "a": {"g": 100, "h": 1},
            "b": {"g": 90, "h": 2},
        }
        markers = select_marker_genes(expr, ["a", "b"], n_genes=2)
        assert len(markers) == len(set(markers))


class TestHelperExtractors:
    def test_get_mean_expression_scalar_and_list(self) -> None:
        assert _get_mean_expression({"g": 5}, "g") == 5.0
        assert _get_mean_expression({"g": [2, 4, 6]}, "g") == 4.0
        assert _get_mean_expression({"g": []}, "g") == 0.0
        assert _get_mean_expression({}, "missing") == 0.0

    def test_get_expression_list(self) -> None:
        assert _get_expression_list({"g": [1, 2]}, "g") == [1.0, 2.0]
        assert _get_expression_list({"g": 7}, "g") == [7.0]
        assert _get_expression_list({}, "g") == [0.0]

    def test_simple_t_test_direction(self) -> None:
        t = _simple_t_test([10.0, 11.0, 12.0], [1.0, 2.0, 3.0])
        assert t > 0
        assert _simple_t_test([1.0, 2.0, 3.0], [10.0, 11.0, 12.0]) < 0
        assert _simple_t_test([1.0], [2.0]) == 0.0  # too few samples
        assert _simple_t_test([5.0, 5.0], [5.0, 5.0]) == 0.0  # zero variance


class TestValidateDeconvolution:
    def test_internal_consistency(self) -> None:
        result = validate_deconvolution({0: 0.6, 1: 0.3, 2: 0.1})
        assert result["valid"] is True
        assert result["sum_proportions"] == pytest.approx(1.0)
        assert result["n_cell_types"] == 3
        assert result["n_nonzero"] == 3
        assert result["max_proportion"] == 0.6
        assert result["max_proportion_type"] == "0"
        assert result["min_proportion"] == pytest.approx(0.1)
        # Entropy of uniform-ish distribution is high but below log2(3)
        assert 0.0 < result["entropy"] <= math.log2(3) + 1e-9

    def test_invalid_negative_proportions(self) -> None:
        result = validate_deconvolution({0: 1.2, 1: -0.2})
        assert result["valid"] is False

    def test_sum_range_validity(self) -> None:
        # Sum 0.6 lies inside the [0.5, 1.5] tolerance window -> valid
        result = validate_deconvolution({0: 0.3, 1: 0.3})
        assert result["valid"] is True
        # Sum 4.0 far out of range -> invalid
        result2 = validate_deconvolution({0: 2.0, 1: 2.0})
        assert result2["valid"] is False

    def test_empty_estimated(self) -> None:
        result = validate_deconvolution({})
        assert result["valid"] is False
        assert result["n_cell_types"] == 0

    def test_ground_truth_metrics(self) -> None:
        result = validate_deconvolution({0: 0.6, 1: 0.4}, ground_truth={0: 0.6, 1: 0.4})
        assert result["correlation"] == pytest.approx(1.0)
        assert result["rmse"] == 0.0
        assert result["mae"] == 0.0
        assert result["max_error"] == 0.0

    def test_ground_truth_imperfect(self) -> None:
        result = validate_deconvolution({0: 0.8, 1: 0.2}, ground_truth={0: 0.6, 1: 0.4})
        assert result["rmse"] > 0.0
        assert result["max_error"] == pytest.approx(0.2)
        assert result["max_error_type"] in ("0", "1")
        assert -1.0 <= result["correlation"] <= 1.0

    def test_ground_truth_extra_keys(self) -> None:
        result = validate_deconvolution({0: 1.0}, ground_truth={0: 0.5, 1: 0.5})
        assert result["rmse"] > 0.0  # missing key treated as 0.0


class TestBatchDeconvolve:
    def test_batch_nnls(self) -> None:
        samples = {
            "s1": MIXTURE,
            "s2": [0.5 * v for v in MIXTURE],
        }
        sig = {"matrix": SIGNATURE, "gene_ids": [f"g{i}" for i in range(6)], "cell_types": ["T_cell", "B_cell"]}
        result = batch_deconvolve(samples, sig, method="nnls")
        assert result["n_samples"] == 2
        assert result["n_cell_types"] == 2
        assert result["failed_samples"] == []
        assert set(result["proportions"]["s1"]) == {"T_cell", "B_cell"}
        assert result["proportions"]["s1"]["T_cell"] == pytest.approx(0.6, abs=1e-3)
        # Summary statistics present per cell type
        assert "T_cell" in result["summary"]
        assert set(result["summary"]["T_cell"]) == {"mean", "std", "min", "max"}

    def test_batch_with_dict_expression(self) -> None:
        gene_ids = [f"g{i}" for i in range(6)]
        sig = {"matrix": SIGNATURE, "gene_ids": gene_ids, "cell_types": ["A", "B"]}
        sample_dict = {gid: val for gid, val in zip(gene_ids, MIXTURE)}
        result = batch_deconvolve({"s1": sample_dict}, sig, method="nnls")
        assert result["proportions"]["s1"]["A"] == pytest.approx(0.6, abs=1e-3)

    def test_batch_tracks_failures(self) -> None:
        samples = {"good": MIXTURE, "bad_dim": [1.0, 2.0]}  # wrong length
        result = batch_deconvolve(samples, SIGNATURE, method="nnls")
        assert result["failed_samples"] == ["bad_dim"]
        assert result["n_samples"] == 1

    def test_batch_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown method"):
            batch_deconvolve({"s": MIXTURE}, SIGNATURE, method="wizardry")

    def test_batch_empty_samples_raises(self) -> None:
        with pytest.raises(ValueError, match="No samples"):
            batch_deconvolve({}, SIGNATURE)

    def test_batch_svr_method(self) -> None:
        samples = {"s1": MIXTURE}
        sig = {"matrix": SIGNATURE, "gene_ids": [f"g{i}" for i in range(6)], "cell_types": ["T", "B"]}
        result = batch_deconvolve(samples, sig, method="svr")
        assert result["method"] == "svr"
        assert result["n_samples"] == 1
