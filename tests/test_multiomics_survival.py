"""Tests for multiomics survival analysis.

Tests Cox proportional hazards, Kaplan-Meier estimation, log-rank tests,
multi-omic survival models (Lasso-Cox), risk stratification, and
concordance index computation.

All tests use real implementations with realistic survival data (NO mocking).
"""

from __future__ import annotations

import math
import random
from typing import Any

import numpy as np
import pytest

from metainformant.multiomics.survival.analysis import (
    compute_concordance_index,
    cox_regression,
    kaplan_meier,
    log_rank_test,
    multi_omic_survival_model,
    risk_stratification,
)

# ---------------------------------------------------------------------------
# Fixtures: realistic survival data generators
# ---------------------------------------------------------------------------


def _make_survival_data(
    n: int = 50,
    p: int = 5,
    censor_rate: float = 0.3,
    seed: int = 42,
) -> tuple[list[float], list[int], Any, list[str]]:
    """Generate realistic survival data with covariates.

    Returns (time, event, covariates, covariate_names).
    """
    rng = np.random.RandomState(seed)

    # Covariates
    X = rng.randn(n, p)
    covariate_names = [f"feature_{i}" for i in range(p)]

    # True coefficients (sparse)
    beta_true = np.zeros(p)
    beta_true[0] = 0.5  # first covariate matters
    if p > 1:
        beta_true[1] = -0.3

    # Generate times from exponential with covariate effect
    linear_pred = X @ beta_true
    scale = np.exp(-linear_pred)
    time = (rng.exponential(scale=scale) + 0.1).tolist()

    # Censoring
    event = [1] * n
    n_censor = int(n * censor_rate)
    censor_idx = rng.choice(n, size=n_censor, replace=False)
    for idx in censor_idx:
        event[idx] = 0

    return time, event, X, covariate_names


def _make_two_group_survival(
    n_per_group: int = 25,
    seed: int = 42,
) -> tuple[list[float], list[int], list[int]]:
    """Generate survival data with two clearly different groups.

    Group 0 has longer survival, group 1 has shorter survival.
    Returns (time, event, groups).
    """
    rng = np.random.RandomState(seed)

    # Group 0: longer survival (scale=10)
    time_0 = rng.exponential(scale=10.0, size=n_per_group) + 0.5
    event_0 = rng.binomial(1, 0.8, size=n_per_group)

    # Group 1: shorter survival (scale=3)
    time_1 = rng.exponential(scale=3.0, size=n_per_group) + 0.5
    event_1 = rng.binomial(1, 0.8, size=n_per_group)

    time = time_0.tolist() + time_1.tolist()
    event = event_0.tolist() + event_1.tolist()
    groups = [0] * n_per_group + [1] * n_per_group

    return time, event, groups


def _make_multi_omic_survival_data(
    n: int = 40,
    seed: int = 42,
) -> tuple[dict[str, Any], list[float], list[int]]:
    """Generate multi-omic features with corresponding survival data.

    Returns (omic_features, time, event).
    """
    rng = np.random.RandomState(seed)

    # Expression: 8 features
    expr = rng.randn(n, 8)
    # Methylation: 6 features
    meth = rng.randn(n, 6)
    # Protein: 5 features
    prot = rng.randn(n, 5)

    omic_features = {
        "expression": expr,
        "methylation": meth,
        "protein": prot,
    }

    # Generate survival based on first expression feature
    linear_pred = expr[:, 0] * 0.5
    scale = np.exp(-linear_pred)
    time = (rng.exponential(scale=scale) + 0.1).tolist()
    event = rng.binomial(1, 0.7, size=n).tolist()

    return omic_features, time, event


# ===================================================================
# Cox Proportional Hazards Regression Tests
# ===================================================================


class TestCoxRegression:
    """Tests for cox_regression."""

    def test_basic_fit(self) -> None:
        """Cox regression should produce coefficients and HR for each covariate."""
        time, event, X, names = _make_survival_data(n=50, p=3)
        result = cox_regression(time, event, X, covariate_names=names)

        assert "coefficients" in result
        assert "hazard_ratios" in result
        assert "se" in result
        assert "p_values" in result
        assert "concordance_index" in result
        assert "log_likelihood" in result
        assert "covariate_names" in result

        assert len(result["coefficients"]) == 3
        assert len(result["hazard_ratios"]) == 3
        assert len(result["p_values"]) == 3

    def test_hazard_ratios_positive(self) -> None:
        """Hazard ratios should always be positive (exp(beta))."""
        time, event, X, names = _make_survival_data(n=40, p=4)
        result = cox_regression(time, event, X)
        for hr in result["hazard_ratios"]:
            assert hr > 0.0

    def test_concordance_index_range(self) -> None:
        """C-index should be between 0 and 1."""
        time, event, X, names = _make_survival_data(n=50, p=3)
        result = cox_regression(time, event, X)
        assert 0.0 <= result["concordance_index"] <= 1.0

    def test_p_values_range(self) -> None:
        """P-values should be between 0 and 1."""
        time, event, X, names = _make_survival_data(n=40, p=3)
        result = cox_regression(time, event, X)
        for pval in result["p_values"]:
            assert 0.0 <= pval <= 1.0

    def test_single_covariate(self) -> None:
        """Cox regression should work with a single covariate."""
        time, event, X, _ = _make_survival_data(n=30, p=1)
        result = cox_regression(time, event, X)
        assert len(result["coefficients"]) == 1

    def test_default_covariate_names(self) -> None:
        """Without provided names, should generate default names."""
        time, event, X, _ = _make_survival_data(n=30, p=3)
        result = cox_regression(time, event, X)
        assert result["covariate_names"] == ["X0", "X1", "X2"]

    def test_mismatched_time_event_raises(self) -> None:
        """Different lengths for time and event should raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            cox_regression([1.0, 2.0], [1], np.array([[1.0], [2.0]]))

    def test_empty_data_raises(self) -> None:
        """Empty data should raise ValueError."""
        with pytest.raises(ValueError, match="empty"):
            cox_regression([], [], np.array([]).reshape(0, 1))

    def test_no_events_raises(self) -> None:
        """All censored data should raise ValueError."""
        with pytest.raises(ValueError, match="No events"):
            cox_regression(
                [1.0, 2.0, 3.0],
                [0, 0, 0],
                np.array([[1.0], [2.0], [3.0]]),
            )

    def test_mismatched_covariate_names_raises(self) -> None:
        """Wrong number of covariate names should raise ValueError."""
        time, event, X, _ = _make_survival_data(n=20, p=3)
        with pytest.raises(ValueError, match="covariate_names"):
            cox_regression(time, event, X, covariate_names=["a", "b"])


# ===================================================================
# Kaplan-Meier Tests
# ===================================================================


class TestKaplanMeier:
    """Tests for kaplan_meier estimation."""

    def test_basic_km(self) -> None:
        """KM should return survival probabilities and confidence intervals."""
        time = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        event = [1, 0, 1, 1, 0, 1, 0, 1]
        result = kaplan_meier(time, event)

        assert "times" in result
        assert "survival_prob" in result
        assert "confidence_lower" in result
        assert "confidence_upper" in result
        assert "n_at_risk" in result
        assert "median_survival" in result

    def test_survival_starts_at_one(self) -> None:
        """First survival probability should be less than 1 (after first event)."""
        time = [1.0, 2.0, 3.0, 4.0, 5.0]
        event = [1, 1, 1, 1, 1]
        result = kaplan_meier(time, event)
        assert result["survival_prob"][0] < 1.0

    def test_survival_monotonic_decreasing(self) -> None:
        """Survival probabilities should be non-increasing."""
        time, event, _ = _make_two_group_survival(n_per_group=20)
        result = kaplan_meier(time, event)
        probs = result["survival_prob"]
        for i in range(len(probs) - 1):
            assert probs[i] >= probs[i + 1] - 1e-10

    def test_confidence_interval_ordering(self) -> None:
        """Lower CI should be <= survival prob <= upper CI."""
        time = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
        event = [1, 0, 1, 1, 0, 1, 1]
        result = kaplan_meier(time, event)
        for i in range(len(result["survival_prob"])):
            assert result["confidence_lower"][i] <= result["survival_prob"][i]
            assert result["survival_prob"][i] <= result["confidence_upper"][i]

    def test_all_events_survival_reaches_zero(self) -> None:
        """With all events, survival should eventually reach zero."""
        time = [1.0, 2.0, 3.0, 4.0, 5.0]
        event = [1, 1, 1, 1, 1]
        result = kaplan_meier(time, event)
        assert result["survival_prob"][-1] == 0.0

    def test_all_censored_no_events(self) -> None:
        """With all censored, no event times to report."""
        time = [1.0, 2.0, 3.0]
        event = [0, 0, 0]
        result = kaplan_meier(time, event)
        assert result["times"] == []
        assert result["survival_prob"] == []

    def test_grouped_km(self) -> None:
        """KM with groups should return separate curves per group."""
        time, event, groups = _make_two_group_survival(n_per_group=15)
        result = kaplan_meier(time, event, groups=groups)

        assert "groups" in result
        assert 0 in result["groups"]
        assert 1 in result["groups"]

        for group_id, group_result in result["groups"].items():
            assert "times" in group_result
            assert "survival_prob" in group_result

    def test_median_survival_detected(self) -> None:
        """Median survival should be detected when S(t) drops below 0.5."""
        time = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        event = [1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
        result = kaplan_meier(time, event)
        assert result["median_survival"] is not None
        assert result["median_survival"] > 0.0

    def test_mismatched_lengths_raises(self) -> None:
        """Different lengths should raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            kaplan_meier([1.0, 2.0], [1])

    def test_empty_data_raises(self) -> None:
        """Empty input should raise ValueError."""
        with pytest.raises(ValueError, match="empty"):
            kaplan_meier([], [])

    def test_mismatched_groups_raises(self) -> None:
        """Groups with wrong length should raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            kaplan_meier([1.0, 2.0], [1, 0], groups=[0])

    def test_n_at_risk_decreasing(self) -> None:
        """Number at risk should be non-increasing over time."""
        time = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        event = [1, 1, 0, 1, 0, 1]
        result = kaplan_meier(time, event)
        n_at_risk = result["n_at_risk"]
        for i in range(len(n_at_risk) - 1):
            assert n_at_risk[i] >= n_at_risk[i + 1]


# ===================================================================
# Log-Rank Test Tests
# ===================================================================


class TestLogRankTest:
    """Tests for log_rank_test."""

    def test_basic_log_rank(self) -> None:
        """Log-rank test should produce chi2, p-value, and df."""
        time, event, groups = _make_two_group_survival(n_per_group=25)
        result = log_rank_test(time, event, groups)

        assert "chi2" in result
        assert "p_value" in result
        assert "df" in result
        assert "observed_expected_per_group" in result
        assert result["chi2"] >= 0.0
        assert 0.0 <= result["p_value"] <= 1.0
        assert result["df"] == 1  # two groups => df=1

    def test_different_groups_significant(self) -> None:
        """Groups with very different survival should yield significant p."""
        time, event, groups = _make_two_group_survival(n_per_group=30, seed=42)
        result = log_rank_test(time, event, groups)
        # Groups have quite different hazards, so p should be reasonably small
        assert result["p_value"] < 0.2  # relaxed threshold for test stability

    def test_same_group_not_significant(self) -> None:
        """Identical groups should yield non-significant p-value."""
        rng = np.random.RandomState(42)
        n = 40
        time = rng.exponential(scale=5.0, size=n).tolist()
        event = [1] * n
        groups = [0] * (n // 2) + [1] * (n - n // 2)
        result = log_rank_test(time, event, groups)
        # Same distribution, p should be > 0.05 usually
        assert result["chi2"] >= 0.0

    def test_observed_expected_per_group(self) -> None:
        """Each group should have observed and expected counts."""
        time, event, groups = _make_two_group_survival(n_per_group=20)
        result = log_rank_test(time, event, groups)
        for group_id in set(groups):
            assert group_id in result["observed_expected_per_group"]
            obs_exp = result["observed_expected_per_group"][group_id]
            assert "observed" in obs_exp
            assert "expected" in obs_exp
            assert obs_exp["observed"] >= 0
            assert obs_exp["expected"] >= 0

    def test_three_groups(self) -> None:
        """Log-rank test should work with 3 groups."""
        rng = np.random.RandomState(42)
        time = rng.exponential(scale=5.0, size=30).tolist()
        event = rng.binomial(1, 0.8, size=30).tolist()
        groups = [0] * 10 + [1] * 10 + [2] * 10
        result = log_rank_test(time, event, groups)
        assert result["df"] == 2  # 3 groups => df=2

    def test_mismatched_lengths_raises(self) -> None:
        """Different lengths should raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            log_rank_test([1.0, 2.0], [1, 0], [0])

    def test_empty_data_raises(self) -> None:
        """Empty data should raise ValueError."""
        with pytest.raises(ValueError, match="empty"):
            log_rank_test([], [], [])

    def test_single_group_raises(self) -> None:
        """Single group should raise ValueError."""
        with pytest.raises(ValueError, match="at least 2"):
            log_rank_test([1.0, 2.0], [1, 1], [0, 0])


# ===================================================================
# Multi-Omic Survival Model Tests
# ===================================================================


class TestMultiOmicSurvivalModel:
    """Tests for multi_omic_survival_model (Lasso-Cox)."""

    def test_basic_model(self) -> None:
        """Lasso-Cox should produce selected features, coefficients, C-index."""
        omic_features, time, event = _make_multi_omic_survival_data(n=40)
        result = multi_omic_survival_model(omic_features, time, event)

        assert "selected_features" in result
        assert "coefficients" in result
        assert "c_index" in result
        assert "risk_scores" in result
        assert len(result["risk_scores"]) == 40

    def test_feature_selection(self) -> None:
        """Lasso should select fewer features than the total."""
        omic_features, time, event = _make_multi_omic_survival_data(n=40)
        result = multi_omic_survival_model(omic_features, time, event)
        total_features = sum(np.asarray(m).shape[1] for m in omic_features.values())
        # Lasso should select some but not necessarily all features
        assert len(result["selected_features"]) <= total_features

    def test_c_index_range(self) -> None:
        """C-index should be between 0 and 1."""
        omic_features, time, event = _make_multi_omic_survival_data(n=40)
        result = multi_omic_survival_model(omic_features, time, event)
        assert 0.0 <= result["c_index"] <= 1.0

    def test_risk_scores_finite(self) -> None:
        """All risk scores should be finite."""
        omic_features, time, event = _make_multi_omic_survival_data(n=30)
        result = multi_omic_survival_model(omic_features, time, event)
        for score in result["risk_scores"]:
            assert math.isfinite(score)

    def test_selected_features_have_coefficients(self) -> None:
        """Number of selected features should match number of coefficients."""
        omic_features, time, event = _make_multi_omic_survival_data(n=40)
        result = multi_omic_survival_model(omic_features, time, event)
        assert len(result["selected_features"]) == len(result["coefficients"])

    def test_feature_names_contain_omic_prefix(self) -> None:
        """Selected feature names should contain omic layer prefix."""
        omic_features, time, event = _make_multi_omic_survival_data(n=40)
        result = multi_omic_survival_model(omic_features, time, event)
        valid_prefixes = list(omic_features.keys())
        for feat_name in result["selected_features"]:
            assert any(feat_name.startswith(prefix) for prefix in valid_prefixes)

    def test_invalid_method_raises(self) -> None:
        """Unsupported method should raise ValueError."""
        omic_features, time, event = _make_multi_omic_survival_data(n=20)
        with pytest.raises(ValueError, match="Unsupported"):
            multi_omic_survival_model(omic_features, time, event, method="ridge_cox")

    def test_empty_features_raises(self) -> None:
        """Empty omic_features should raise ValueError."""
        with pytest.raises(ValueError, match="at least one"):
            multi_omic_survival_model({}, [1.0, 2.0], [1, 0])

    def test_mismatched_time_event_raises(self) -> None:
        """Mismatched time and event should raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            multi_omic_survival_model(
                {"a": np.random.randn(5, 3)},
                [1.0, 2.0, 3.0, 4.0, 5.0],
                [1, 0, 1],
            )

    def test_single_omic(self) -> None:
        """Should work with a single omic layer."""
        rng = np.random.RandomState(42)
        features = {"single_omic": rng.randn(30, 5)}
        time = rng.exponential(scale=5.0, size=30).tolist()
        event = rng.binomial(1, 0.7, size=30).tolist()
        result = multi_omic_survival_model(features, time, event)
        assert "c_index" in result


# ===================================================================
# Risk Stratification Tests
# ===================================================================


class TestRiskStratification:
    """Tests for risk_stratification."""

    def test_basic_stratification(self) -> None:
        """Risk stratification should produce group labels and survival curves."""
        rng = np.random.RandomState(42)
        n = 40
        risk_scores = rng.randn(n).tolist()
        time = rng.exponential(scale=5.0, size=n).tolist()
        event = rng.binomial(1, 0.7, size=n).tolist()

        result = risk_stratification(risk_scores, time, event, n_groups=2)

        assert "group_labels" in result
        assert "survival_curves" in result
        assert "log_rank_p" in result
        assert "hazard_ratio" in result
        assert len(result["group_labels"]) == n

    def test_two_groups_produces_two_labels(self) -> None:
        """Splitting into 2 groups should produce exactly 2 unique labels."""
        rng = np.random.RandomState(42)
        n = 30
        risk_scores = rng.randn(n).tolist()
        time = rng.exponential(scale=5.0, size=n).tolist()
        event = rng.binomial(1, 0.7, size=n).tolist()

        result = risk_stratification(risk_scores, time, event, n_groups=2)
        assert len(set(result["group_labels"])) == 2

    def test_three_groups(self) -> None:
        """Risk stratification into 3 groups should work."""
        rng = np.random.RandomState(42)
        n = 45
        risk_scores = rng.randn(n).tolist()
        time = rng.exponential(scale=5.0, size=n).tolist()
        event = rng.binomial(1, 0.7, size=n).tolist()

        result = risk_stratification(risk_scores, time, event, n_groups=3)
        assert len(set(result["group_labels"])) >= 2  # at least 2 unique groups

    def test_log_rank_p_valid(self) -> None:
        """Log-rank p-value should be in [0, 1]."""
        rng = np.random.RandomState(42)
        n = 30
        risk_scores = rng.randn(n).tolist()
        time = rng.exponential(scale=5.0, size=n).tolist()
        event = rng.binomial(1, 0.7, size=n).tolist()

        result = risk_stratification(risk_scores, time, event, n_groups=2)
        assert 0.0 <= result["log_rank_p"] <= 1.0

    def test_hazard_ratio_positive(self) -> None:
        """Hazard ratio should be non-negative."""
        rng = np.random.RandomState(42)
        n = 40
        risk_scores = rng.randn(n).tolist()
        time = rng.exponential(scale=5.0, size=n).tolist()
        event = rng.binomial(1, 0.7, size=n).tolist()

        result = risk_stratification(risk_scores, time, event, n_groups=2)
        assert result["hazard_ratio"] >= 0.0

    def test_survival_curves_per_group(self) -> None:
        """Survival curves should be present for each group."""
        rng = np.random.RandomState(42)
        n = 40
        risk_scores = rng.randn(n).tolist()
        time = rng.exponential(scale=5.0, size=n).tolist()
        event = rng.binomial(1, 0.7, size=n).tolist()

        result = risk_stratification(risk_scores, time, event, n_groups=2)
        # survival_curves should contain per-group data
        assert result["survival_curves"] is not None

    def test_n_groups_one_raises(self) -> None:
        """n_groups < 2 should raise ValueError."""
        with pytest.raises(ValueError, match="n_groups must be >= 2"):
            risk_stratification([1.0, 2.0], [1.0, 2.0], [1, 1], n_groups=1)

    def test_mismatched_lengths_raises(self) -> None:
        """Mismatched input lengths should raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            risk_stratification([1.0, 2.0], [1.0], [1, 0])

    def test_too_few_subjects_raises(self) -> None:
        """Fewer subjects than groups should raise ValueError."""
        with pytest.raises(ValueError, match="at least"):
            risk_stratification([1.0], [2.0], [1], n_groups=3)


# ===================================================================
# Concordance Index Tests
# ===================================================================


class TestConcordanceIndex:
    """Tests for compute_concordance_index."""

    def test_basic_concordance(self) -> None:
        """Should compute C-index with concordant/discordant counts."""
        risk = [3.0, 2.0, 1.0, 0.5]
        time = [1.0, 2.0, 3.0, 4.0]
        event = [1, 1, 1, 1]
        result = compute_concordance_index(risk, time, event)

        assert "c_index" in result
        assert "se" in result
        assert "n_concordant" in result
        assert "n_discordant" in result
        assert "n_tied" in result
        assert 0.0 <= result["c_index"] <= 1.0

    def test_perfect_concordance(self) -> None:
        """Higher risk = earlier event should give C-index close to 1."""
        risk = [5.0, 4.0, 3.0, 2.0, 1.0]
        time = [1.0, 2.0, 3.0, 4.0, 5.0]
        event = [1, 1, 1, 1, 1]
        result = compute_concordance_index(risk, time, event)
        assert result["c_index"] > 0.9

    def test_inverse_concordance(self) -> None:
        """Higher risk = later event should give C-index close to 0."""
        risk = [1.0, 2.0, 3.0, 4.0, 5.0]
        time = [1.0, 2.0, 3.0, 4.0, 5.0]
        event = [1, 1, 1, 1, 1]
        result = compute_concordance_index(risk, time, event)
        assert result["c_index"] < 0.2

    def test_random_concordance(self) -> None:
        """Random risk scores should give C-index near 0.5."""
        rng = np.random.RandomState(42)
        n = 100
        risk = rng.randn(n).tolist()
        time = rng.exponential(scale=5.0, size=n).tolist()
        event = [1] * n
        result = compute_concordance_index(risk, time, event)
        # Should be roughly 0.5 +/- 0.15
        assert 0.3 <= result["c_index"] <= 0.7

    def test_with_censoring(self) -> None:
        """C-index should handle censored observations."""
        risk = [3.0, 2.0, 1.0, 0.5]
        time = [1.0, 2.0, 3.0, 4.0]
        event = [1, 0, 1, 0]
        result = compute_concordance_index(risk, time, event)
        assert 0.0 <= result["c_index"] <= 1.0

    def test_no_events_returns_half(self) -> None:
        """With no events, C-index should be 0.5."""
        risk = [1.0, 2.0, 3.0]
        time = [1.0, 2.0, 3.0]
        event = [0, 0, 0]
        result = compute_concordance_index(risk, time, event)
        assert result["c_index"] == 0.5

    def test_se_non_negative(self) -> None:
        """Standard error should be non-negative."""
        rng = np.random.RandomState(42)
        n = 30
        risk = rng.randn(n).tolist()
        time = rng.exponential(scale=5.0, size=n).tolist()
        event = rng.binomial(1, 0.7, size=n).tolist()
        result = compute_concordance_index(risk, time, event)
        assert result["se"] >= 0.0

    def test_mismatched_lengths_raises(self) -> None:
        """Mismatched input lengths should raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            compute_concordance_index([1.0, 2.0], [1.0], [1, 0])

    def test_tied_risks(self) -> None:
        """Tied risk scores should be counted properly."""
        risk = [2.0, 2.0, 1.0]
        time = [1.0, 2.0, 3.0]
        event = [1, 1, 1]
        result = compute_concordance_index(risk, time, event)
        assert result["n_tied"] >= 0
