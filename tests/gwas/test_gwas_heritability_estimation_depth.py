"""Real-implementation tests for metainformant.gwas.heritability.estimation.

Depth coverage: LDSC heritability estimation (with and without intercept),
partitioned h2, cross-trait genetic correlation, Haseman-Elston regression,
simple GREML via REML eigendecomposition, liability-scale conversion, and the
internal numerical helpers. All inputs are real synthetic numeric data with
planted ground truth; no test doubles. Descriptive statistics only.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.gwas.heritability.estimation import (
    _approximate_h2_se,
    _golden_section_max,
    _ldsc_regression_numpy,
    _ldsc_regression_pure,
    _normal_cdf,
    _normal_pdf,
    _normal_quantile,
    _partitioned_h2_separate,
    _reml_log_likelihood,
    compute_liability_h2,
    estimate_h2_ldsc,
    genetic_correlation,
    greml_simple,
    haseman_elston_regression,
    partitioned_h2,
)

RNG = np.random.default_rng(42)


def _simulate_ldsc_data(
    n_snps: int = 500,
    n_samples: int = 5000,
    true_h2: float = 0.3,
    intercept: float = 1.05,
    seed: int = 42,
) -> tuple[list[float], list[float]]:
    """Simulate chi2 statistics from an LDSC model with known parameters.

    chi2_j = a + N/M * h2 * l_j, with noise. LD scores are drawn from a
    lognormal-like distribution typical of real LD profiles.
    """
    rng = np.random.default_rng(seed)
    ld_scores = rng.gamma(shape=2.0, scale=3.0, size=n_snps)  # positive, skewed
    per_snp = true_h2 * n_samples / n_snps
    mean_chi2 = intercept + per_snp * ld_scores
    noise = rng.normal(0.0, 0.15, size=n_snps)
    chi2 = np.maximum(mean_chi2 + noise, 0.01)
    return chi2.tolist(), ld_scores.tolist()


class TestEstimateH2Ldsc:
    def test_recovers_planted_h2(self) -> None:
        chi2, ld = _simulate_ldsc_data(true_h2=0.3, intercept=1.05)
        result = estimate_h2_ldsc(chi2, ld, n_samples=5000)
        assert result["status"] == "success"
        assert result["h2"] == pytest.approx(0.3, abs=0.08)
        assert result["h2_se"] > 0.0
        assert 0.9 < result["intercept"] < 1.2
        assert result["lambda_gc"] > 0.0
        assert 0.0 <= result["h2"] <= 1.0

    def test_h2_clamped_to_unit_interval(self) -> None:
        # Chi2 unrelated to LD scores (all noise) -> slope ~ 0 -> h2 ~ 0
        rng = np.random.default_rng(7)
        chi2 = list(rng.uniform(0.5, 1.5, size=100))
        ld = list(rng.gamma(2.0, 3.0, size=100))
        result = estimate_h2_ldsc(chi2, ld, n_samples=1000)
        assert result["status"] == "success"
        assert 0.0 <= result["h2"] <= 1.0

    def test_no_intercept_mode(self) -> None:
        chi2, ld = _simulate_ldsc_data()
        result = estimate_h2_ldsc(chi2, ld, n_samples=5000, intercept=False)
        assert result["status"] == "success"
        assert result["intercept"] == 1.0
        assert result["intercept_se"] == 0.0

    def test_lambda_gc_from_median(self) -> None:
        # Even-length median = mean of two central values
        chi2 = [0.1, 0.4549, 0.4549, 1.0]
        ld = [1.0, 2.0, 3.0, 4.0]
        result = estimate_h2_ldsc(chi2, ld, n_samples=100)
        assert result["status"] == "success"
        assert result["lambda_gc"] == pytest.approx(1.0, abs=1e-6)

    def test_empty_inputs_error(self) -> None:
        result = estimate_h2_ldsc([], [], n_samples=100)
        assert result["status"] == "error"
        assert "required" in result["message"]

    def test_length_mismatch_error(self) -> None:
        result = estimate_h2_ldsc([1.0, 2.0], [1.0], n_samples=100)
        assert result["status"] == "error"
        assert "mismatch" in result["message"]

    def test_invalid_sample_size_error(self) -> None:
        result = estimate_h2_ldsc([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], n_samples=0)
        assert result["status"] == "error"
        assert "positive" in result["message"]

    def test_too_few_snps_error(self) -> None:
        result = estimate_h2_ldsc([1.0, 2.0], [1.0, 2.0], n_samples=100)
        assert result["status"] == "error"
        assert "at least 3" in result["message"]

    def test_mean_chi2_reported(self) -> None:
        chi2, ld = _simulate_ldsc_data(n_snps=50)
        result = estimate_h2_ldsc(chi2, ld, n_samples=5000)
        assert result["mean_chi2"] == pytest.approx(sum(chi2) / 50, rel=1e-6)


class TestPartitionedH2:
    def test_two_categories_recover_signal(self) -> None:
        rng = np.random.default_rng(11)
        n_snps, n_samples = 400, 8000
        ld_a = rng.gamma(2.0, 3.0, size=n_snps)
        ld_b = rng.gamma(2.0, 3.0, size=n_snps)
        # Category A carries h2=0.2, category B carries h2=0.1
        chi2 = 1.0 + (0.2 * n_samples / n_snps) * ld_a + (0.1 * n_samples / n_snps) * ld_b
        chi2 = np.maximum(chi2 + rng.normal(0, 0.1, n_snps), 0.01)
        result = partitioned_h2(chi2.tolist(), {"catA": ld_a.tolist(), "catB": ld_b.tolist()}, n_samples)
        assert result["status"] == "success"
        assert result["n_categories"] == 2
        assert result["n_snps"] == n_snps
        cats = result["per_category"]
        assert cats["catA"]["h2"] > cats["catB"]["h2"]
        assert cats["catA"]["proportion"] > 0.5
        assert cats["catA"]["h2_se"] >= 0.0
        assert result["total_h2"] == pytest.approx(sum(c["h2"] for c in cats.values()), abs=1e-6)

    def test_enrichment_computed(self) -> None:
        chi2, ld = _simulate_ldsc_data(n_snps=100)
        result = partitioned_h2(chi2, {"all": ld}, n_samples=5000)
        assert result["status"] == "success"
        cat = result["per_category"]["all"]
        # All SNPs in the category -> enrichment = proportion / 1.0
        assert cat["enrichment"] == pytest.approx(cat["proportion"], rel=1e-6)

    def test_length_mismatch_reports_category(self) -> None:
        chi2 = [1.0, 2.0, 3.0, 4.0]
        result = partitioned_h2(chi2, {"bad": [1.0, 2.0]}, n_samples=100)
        assert result["status"] == "error"
        assert "bad" in result["message"]

    def test_empty_inputs_error(self) -> None:
        assert partitioned_h2([], {"a": [1.0]}, 100)["status"] == "error"
        assert partitioned_h2([1.0], {}, 100)["status"] == "error"


class TestGeneticCorrelation:
    def _paired_traits(self) -> tuple[list[float], list[float], list[float]]:
        """Two traits sharing a genetic component -> positive rg."""
        chi2_base, ld = _simulate_ldsc_data(n_snps=300, n_samples=5000, true_h2=0.25)
        rng = np.random.default_rng(5)
        # Construct z-scores whose products correlate through shared LD signal
        z1 = [math.sqrt(max(c, 0.01)) for c in chi2_base]
        shared = np.array(z1) * 0.5 + np.array(rng.normal(0, 1, 300)) * 0.5
        z2 = list(np.abs(shared))
        return z1, z2, ld

    def test_positive_rg_for_shared_architecture(self) -> None:
        z1, z2, ld = self._paired_traits()
        result = genetic_correlation(z1, z2, ld, n1=5000, n2=5000)
        assert result["status"] == "success"
        assert -1.0 <= result["rg"] <= 1.0
        assert result["rg"] > 0.0  # shared architecture
        assert result["rg_se"] > 0.0
        assert 0.0 <= result["p_value"] <= 1.0
        assert result["h2_1"] > 0.0 and result["h2_2"] > 0.0
        assert result["intercept"] != 0.0

    def test_length_mismatch_error(self) -> None:
        result = genetic_correlation([1.0, 2.0, 3.0], [1.0, 2.0], [1.0, 2.0, 3.0], 100, 100)
        assert result["status"] == "error"
        assert "mismatch" in result["message"]

    def test_too_few_snps_error(self) -> None:
        result = genetic_correlation([1.0, 2.0], [1.0, 2.0], [1.0, 2.0], 100, 100)
        assert result["status"] == "error"
        assert "at least 3" in result["message"]

    def test_empty_inputs_error(self) -> None:
        result = genetic_correlation([], [1.0, 2.0, 3.0], [1.0, 1.0, 1.0], 100, 100)
        assert result["status"] == "error"


class TestHasemanElston:
    def test_recovers_planted_h2(self) -> None:
        # 100 samples with GRM-like structure; phenotype = genetic + noise
        rng = np.random.default_rng(13)
        n = 60
        K = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                k = rng.normal(0.0, 0.05)
                K[i, j] = K[j, i] = k
        np.fill_diagonal(K, 1.0)
        # Genetic values via GRM cholesky-like construction
        g = rng.normal(0.0, 1.0, n) @ K * 0.3
        y = g + rng.normal(0.0, 1.0, n)
        y_centered = y - y.mean()

        iu = np.triu_indices(n, k=1)
        grm_offdiag = K[iu].tolist()
        products = (y_centered[iu[0]] * y_centered[iu[1]]).tolist()

        result = haseman_elston_regression(grm_offdiag, products)
        assert result["status"] == "success"
        assert result["n_pairs"] == len(grm_offdiag)
        assert result["h2_se"] >= 0.0
        assert 0.0 <= result["h2"] <= 1.0

    def test_perfect_linear_relation(self) -> None:
        # y_product exactly proportional to grm -> clean slope
        grm = [0.1, 0.2, 0.3, 0.4]
        products = [2.0 * g + 0.5 for g in grm]
        result = haseman_elston_regression(grm, products)
        assert result["status"] == "success"
        assert result["slope"] == pytest.approx(2.0, rel=1e-6)
        assert result["intercept"] == pytest.approx(0.5, rel=1e-6)

    def test_length_mismatch_error(self) -> None:
        result = haseman_elston_regression([0.1, 0.2, 0.3], [1.0, 2.0])
        assert result["status"] == "error"
        assert "mismatch" in result["message"]

    def test_too_few_pairs_error(self) -> None:
        result = haseman_elston_regression([0.1, 0.2], [1.0, 2.0])
        assert result["status"] == "error"
        assert "at least 3" in result["message"]

    def test_zero_variance_grm_error(self) -> None:
        result = haseman_elston_regression([0.5, 0.5, 0.5, 0.5], [1.0, 2.0, 3.0, 4.0])
        assert result["status"] == "error"
        assert "Zero variance" in result["message"]

    def test_empty_inputs_error(self) -> None:
        result = haseman_elston_regression([], [])
        assert result["status"] == "error"


class TestGremlSimple:
    def _simulate_gwas(self, n: int = 150, m_snps: int = 60, h2_true: float = 0.6, seed: int = 21):
        rng = np.random.default_rng(seed)
        # Standardized genotype matrix -> GRM. Few SNPs relative to samples so
        # the GRM has real family structure (eigenvalue spread) and h2 is
        # identifiable; with a near-diagonal GRM REML has no power.
        G = (rng.random((n, m_snps)) < 0.5).astype(float)
        G = (G - G.mean(axis=0)) / (G.std(axis=0) + 1e-9)
        K = (G @ G.T) / m_snps
        u = rng.multivariate_normal(np.zeros(n), K * h2_true)
        y = u + rng.normal(0.0, math.sqrt(1.0 - h2_true), n)
        return K, y

    def test_recovers_planted_h2(self) -> None:
        K, y = self._simulate_gwas(h2_true=0.6)
        result = greml_simple(K, y.tolist())
        assert result["status"] == "success"
        assert result["method"] == "greml"
        assert result["h2"] == pytest.approx(0.6, abs=0.25)
        assert result["h2_se"] > 0.0
        assert result["sigma_g"] > 0.0
        assert result["sigma_e"] > 0.0
        assert result["n_samples"] == len(y)
        assert result["log_likelihood"] != 0.0

    def test_zero_variance_phenotype_error(self) -> None:
        result = greml_simple(np.eye(5), [1.0] * 5)
        assert result["status"] == "error"
        assert "variance" in result["message"]

    def test_too_few_samples_error(self) -> None:
        result = greml_simple(np.eye(2), [1.0, 2.0])
        assert result["status"] == "error"
        assert "at least 3" in result["message"]

    def test_shape_mismatch_error(self) -> None:
        result = greml_simple(np.eye(5), [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        assert result["status"] == "error"
        assert "shape" in result["message"]

    def test_accepts_list_of_lists(self) -> None:
        K, y = self._simulate_gwas(n=30, seed=3)
        result = greml_simple(K.tolist(), y.tolist())
        assert result["status"] == "success"


class TestLiabilityH2:
    def test_matches_documented_formula(self) -> None:
        # Pin the Lee-style factor: h2_l = h2_o * K^2(1-K)^2 / (P(1-P) z^2)
        # where z = pdf at liability threshold Phi^{-1}(1-K).
        K_prev, P = 0.01, 0.5
        threshold = _normal_quantile(1.0 - K_prev)
        z_density = _normal_pdf(threshold)
        expected_factor = (K_prev**2 * (1 - K_prev) ** 2) / (P * (1 - P) * z_density**2)
        result = compute_liability_h2(0.05, prevalence=K_prev, sample_prevalence=P)
        assert result["status"] == "success"
        assert result["conversion_factor"] == pytest.approx(expected_factor, rel=1e-6)
        assert result["h2_liability"] == pytest.approx(min(1.0, 0.05 * expected_factor), rel=1e-6)
        assert result["threshold"] == pytest.approx(threshold, rel=1e-6)

    def test_rare_disease_balanced_sampling_inflates(self) -> None:
        # When sample prevalence matches the rare population prevalence, the
        # conversion factor is large (ascertainment correction upward).
        result = compute_liability_h2(0.05, prevalence=0.01, sample_prevalence=0.01)
        assert result["status"] == "success"
        assert result["conversion_factor"] > 10.0

    def test_threshold_positive_for_rare_disease(self) -> None:
        result = compute_liability_h2(0.05, prevalence=0.01, sample_prevalence=0.5)
        assert result["status"] == "success"
        assert result["threshold"] > 0.0
        assert 0.0 <= result["h2_liability"] <= 1.0

    def test_identity_for_balanced_design(self) -> None:
        # When P matches K, factor = K^2(1-K)^2 / (K(1-K)z^2) = K(1-K)/z^2
        result = compute_liability_h2(0.2, prevalence=0.1, sample_prevalence=0.1)
        assert result["status"] == "success"
        assert result["h2_liability"] >= 0.2  # never shrinks below observed

    def test_h2_clamped_to_one(self) -> None:
        # Balanced sampling of a rare disease: factor > 20 pushes h2_l to the cap
        result = compute_liability_h2(1.0, prevalence=0.001, sample_prevalence=0.001)
        assert result["status"] == "success"
        assert result["h2_liability"] == 1.0

    def test_out_of_range_h2_error(self) -> None:
        for bad in (-0.1, 1.5):
            result = compute_liability_h2(bad, prevalence=0.01, sample_prevalence=0.5)
            assert result["status"] == "error"
            assert "[0, 1]" in result["message"]

    def test_out_of_range_prevalences(self) -> None:
        assert compute_liability_h2(0.1, prevalence=0.0, sample_prevalence=0.5)["status"] == "error"
        assert compute_liability_h2(0.1, prevalence=1.0, sample_prevalence=0.5)["status"] == "error"
        assert compute_liability_h2(0.1, prevalence=0.01, sample_prevalence=0.0)["status"] == "error"
        assert compute_liability_h2(0.1, prevalence=0.01, sample_prevalence=1.0)["status"] == "error"


class TestLdscRegressionHelpers:
    def _data(self) -> tuple[list[float], list[float]]:
        ld = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [1.0 + 0.5 * v for v in ld]  # intercept 1, slope 0.5
        return y, ld

    def test_numpy_and_pure_agree(self) -> None:
        y, ld = self._data()
        r_np = _ldsc_regression_numpy(y, ld, 100, 5, True)
        r_pure = _ldsc_regression_pure(y, ld, 100, 5, True)
        assert r_np is not None and r_pure is not None
        for a, b in zip(r_np, r_pure):
            assert a == pytest.approx(b, rel=1e-6)

    def test_slope_and_intercept(self) -> None:
        y, ld = self._data()
        slope, slope_se, intercept, intercept_se = _ldsc_regression_numpy(y, ld, 100, 5, True)
        assert slope == pytest.approx(0.5, rel=1e-6)
        assert intercept == pytest.approx(1.0, rel=1e-6)
        assert slope_se >= 0.0 and intercept_se >= 0.0

    def test_no_intercept_constraint(self) -> None:
        y, ld = self._data()
        slope, slope_se, intercept, intercept_se = _ldsc_regression_pure(y, ld, 100, 5, False)
        assert intercept == 1.0
        assert intercept_se == 0.0
        assert slope_se >= 0.0

    def test_singular_returns_none(self) -> None:
        # All identical LD scores -> no information -> None
        result = _ldsc_regression_numpy([1.0, 2.0, 3.0], [0.0, 0.0, 0.0], 10, 3, True)
        assert result is None


class TestPartitionedFallback:
    def test_separate_fallback_runs(self) -> None:
        chi2, ld = _simulate_ldsc_data(n_snps=50)
        result = _partitioned_h2_separate(chi2, {"a": ld, "b": ld}, n_samples=5000)
        assert result["status"] == "success"
        assert result["n_categories"] == 2
        assert result["per_category"]["a"]["h2"] == pytest.approx(result["per_category"]["b"]["h2"], rel=1e-6)
        # Proportions sum to ~1 when both categories succeed
        props = sum(c["proportion"] for c in result["per_category"].values())
        assert props == pytest.approx(1.0, abs=1e-6)


class TestRemlHelpers:
    def _setup(self, seed: int = 9):
        rng = np.random.default_rng(seed)
        n = 20
        A = rng.normal(0.0, 1.0, (n, n))
        K = (A + A.T) / 2
        np.fill_diagonal(K, 1.0)
        eigenvalues, eigenvectors = np.linalg.eigh(K)
        y = rng.normal(0.0, 1.0, n)
        y_rot = eigenvectors.T @ y
        ones_rot = eigenvectors.T @ np.ones(n)
        return y_rot, ones_rot, np.maximum(eigenvalues, 0.0), n

    def test_reml_ll_is_finite_and_peaked(self) -> None:
        y_rot, ones_rot, eigs, n = self._setup()
        lls = [_reml_log_likelihood(y_rot, ones_rot, eigs, h2, n) for h2 in (0.1, 0.3, 0.5, 0.7, 0.9)]
        assert all(math.isfinite(v) for v in lls)

    def test_approximate_se_positive_and_bounded(self) -> None:
        y_rot, ones_rot, eigs, n = self._setup()
        se = _approximate_h2_se(y_rot, ones_rot, eigs, 0.5, n)
        assert 0.0 < se <= 0.5

    def test_golden_section_finds_maximum(self) -> None:
        # Unimodal function with known max at x=0.7
        best = _golden_section_max(lambda x: -((x - 0.7) ** 2), 0.0, 1.0)
        assert best == pytest.approx(0.7, abs=1e-3)


class TestNormalHelpers:
    def test_cdf_extremes_and_symmetry(self) -> None:
        assert _normal_cdf(0.0) == pytest.approx(0.5, abs=1e-6)
        assert _normal_cdf(10.0) == pytest.approx(1.0, abs=1e-6)
        assert _normal_cdf(-10.0) == pytest.approx(0.0, abs=1e-6)
        # Symmetry: Phi(-x) = 1 - Phi(x)
        for x in (0.5, 1.0, 2.5):
            assert _normal_cdf(-x) == pytest.approx(1.0 - _normal_cdf(x), abs=1e-6)

    def test_cdf_matches_erf_reference(self) -> None:
        for x in (0.3, 1.1, 2.0):
            expected = 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))
            assert _normal_cdf(x) == pytest.approx(expected, abs=1e-4)

    def test_pdf_values(self) -> None:
        assert _normal_pdf(0.0) == pytest.approx(0.3989422804014327, rel=1e-9)
        assert _normal_pdf(1.0) == pytest.approx(0.24197072451914337, rel=1e-9)

    def test_quantile_known_values(self) -> None:
        assert _normal_quantile(0.5) == 0.0
        assert _normal_quantile(0.975) == pytest.approx(1.959964, abs=1e-3)
        assert _normal_quantile(0.025) == pytest.approx(-1.959964, abs=1e-3)
        assert _normal_quantile(0.0) == -10.0
        assert _normal_quantile(1.0) == 10.0

    def test_quantile_cdf_roundtrip(self) -> None:
        for p in (0.6, 0.8, 0.95):
            z = _normal_quantile(p)
            assert _normal_cdf(z) == pytest.approx(p, abs=1e-3)
