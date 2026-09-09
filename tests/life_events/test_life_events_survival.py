"""Tests for life_events survival and time-to-event analysis."""

from __future__ import annotations

import math

import pytest

from metainformant.life_events.survival.time_to_event import (
    competing_risks,
    cox_ph_model,
    kaplan_meier_estimator,
    recurrent_events,
    time_varying_covariates,
)


class TestKaplanMeier:
    def test_all_events(self):
        """Survival drops multiplicatively when every subject has an event."""
        result = kaplan_meier_estimator(times=[5.0, 10.0, 15.0], events=[1, 1, 1])

        assert result["time_points"] == [5.0, 10.0, 15.0]
        assert result["survival_probability"] == pytest.approx([2 / 3, 1 / 3, 0.0])
        assert result["n_events"] == 3
        assert result["n_censored"] == 0
        assert result["median_survival"] == 10.0
        # Survival never exceeds 1 and CIs stay within [0, 1]
        for prob, (low, high) in zip(result["survival_probability"], result["confidence_interval"]):
            assert 0.0 <= prob <= 1.0
            assert 0.0 <= low <= high <= 1.0

    def test_with_censoring(self):
        """Censored observations are removed from the risk set without dropping survival."""
        result = kaplan_meier_estimator(times=[5.0, 8.0, 10.0], events=[1, 0, 1])

        assert result["survival_probability"] == pytest.approx([2 / 3, 0.0])
        assert result["n_events"] == 2
        assert result["n_censored"] == 1
        assert result["n_at_risk"] == [3, 1]

    def test_no_events_survival_stays_one(self):
        """With no events the curve is not evaluated (no event times)."""
        result = kaplan_meier_estimator(times=[5.0, 10.0], events=[0, 0])

        assert result["time_points"] == []
        assert result["survival_probability"] == []
        assert result["n_events"] == 0
        assert result["n_censored"] == 2
        assert result["median_survival"] is None

    def test_length_mismatch_raises(self):
        """Mismatched times/events lengths raise ValueError."""
        with pytest.raises(ValueError, match="must match"):
            kaplan_meier_estimator(times=[1.0, 2.0], events=[1])

    def test_stratified_groups(self):
        """Stratified analysis returns one curve per group plus aggregated totals."""
        result = kaplan_meier_estimator(
            times=[5.0, 10.0, 5.0, 20.0],
            events=[1, 1, 1, 0],
            groups=["a", "a", "b", "b"],
        )
        assert set(result["groups"]) == {"a", "b"}
        assert result["n_events"] == 3
        assert result["n_censored"] == 1
        assert result["groups"]["a"]["median_survival"] == 5.0  # survival hits 0.5 at t=5
        assert result["time_points"] == [5.0, 10.0]


class TestCoxPH:
    def test_positive_effect_fitted(self):
        """A covariate that raises hazard gets a positive beta and hazard ratio > 1."""
        # X=4 has the shortest observed time -> higher X, higher hazard
        times = [10.0, 20.0, 30.0, 40.0]
        events = [1, 1, 1, 0]
        covariates = [[4.0], [3.0], [2.0], [1.0]]

        result = cox_ph_model(times, events, covariates, covariate_names=["exposure"])

        assert result["n_subjects"] == 4
        assert result["n_events"] == 3
        assert result["coefficients"][0] > 0
        assert result["hazard_ratios"][0] > 1.0
        # Hazard ratio is exp(beta) and consistent with the coefficient
        assert result["hazard_ratios"][0] == pytest.approx(math.exp(result["coefficients"][0]))
        assert 0.0 <= result["concordance"] <= 1.0
        assert 0.0 <= result["p_values"][0] <= 1.0
        assert result["covariate_names"] == ["exposure"]

    def test_default_covariate_names(self):
        """Without names, covariates are labelled X0..Xp-1."""
        result = cox_ph_model([10.0, 20.0], [1, 0], [[1.0, 2.0], [2.0, 1.0]])

        assert result["covariate_names"] == ["X0", "X1"]
        assert len(result["coefficients"]) == 2
        assert len(result["se"]) == 2

    def test_length_mismatch_raises(self):
        """Mismatched input lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            cox_ph_model([1.0, 2.0], [1], [[1.0], [2.0]])

    def test_no_covariates_raises(self):
        """An empty covariate vector raises ValueError."""
        with pytest.raises(ValueError, match="At least one covariate"):
            cox_ph_model([1.0, 2.0], [1, 0], [[], []])


class TestCompetingRisks:
    def test_two_event_types(self):
        """Cumulative incidence accumulates per cause and overall survival drops."""
        result = competing_risks(times=[5.0, 5.0, 10.0], events=[1, 1, 0], event_types=[1, 2, 0])

        ci = result["cumulative_incidence_per_type"]
        assert set(ci) == {1, 2}
        assert ci[1] == pytest.approx([(5.0, 1 / 3)])
        assert ci[2] == pytest.approx([(5.0, 1 / 3)])
        hazards = result["cause_specific_hazards"]
        assert hazards[1] == pytest.approx([(5.0, 1 / 3)])
        assert hazards[2] == pytest.approx([(5.0, 1 / 3)])
        # Survival drops to 1/3 at t=5; the t=10 observation is censored-only,
        # so it is recorded with unchanged survival
        overall = result["overall_survival"]
        assert [t for t, _ in overall] == [5.0, 10.0]
        assert all(abs(s - 1 / 3) < 1e-9 for _, s in overall)

    def test_length_mismatch_raises(self):
        """Mismatched input lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            competing_risks(times=[1.0], events=[1, 0], event_types=[1, 1])


class TestRecurrentEvents:
    def test_mean_cumulative_function_and_rates(self):
        """MCF accumulates per event time and rates are events over observation time."""
        result = recurrent_events(
            times=[[1.0, 3.0, 5.0, 7.0], [2.0, 4.0, 6.0]],
            events=[[1, 0, 1, 0], [1, 1, 1]],
            subject_ids=["a", "b"],
        )

        assert result["n_subjects"] == 2
        assert result["total_events"] == 5
        assert result["rate_per_subject"]["a"] == pytest.approx(2 / 7)
        assert result["rate_per_subject"]["b"] == pytest.approx(0.5)

        # At each unique event time: cumulative events / subjects still at risk
        mcf = result["mean_cumulative_function"]
        assert [t for t, _ in mcf] == [1.0, 2.0, 4.0, 5.0, 6.0]
        assert mcf[0] == pytest.approx((1.0, 0.5))
        assert mcf[-1] == pytest.approx((6.0, 2.5))
        # MCF is non-decreasing
        values = [v for _, v in mcf]
        assert values == sorted(values)

        gaps = result["gap_time_distribution"]
        assert gaps["n_gaps"] == 3
        assert gaps["median"] == pytest.approx(2.0)
        assert gaps["mean"] == pytest.approx(8 / 3)

    def test_subject_with_no_events(self):
        """Subjects without events contribute no events but do not crash."""
        result = recurrent_events(times=[[1.0]], events=[[0]], subject_ids=["lonely"])

        assert result["total_events"] == 0
        assert result["mean_cumulative_function"] == []
        assert result["gap_time_distribution"] == {}
        assert result["rate_per_subject"]["lonely"] == 0.0

    def test_length_mismatch_raises(self):
        """Mismatched input lengths raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            recurrent_events(times=[[1.0]], events=[[1], [0]], subject_ids=["a"])


class TestTimeVaryingCovariates:
    def test_expansion_and_sorting(self):
        """Interval data is expanded to counting-process rows sorted by subject and start."""
        intervals = [
            {"subject_id": "s2", "start": 0.0, "stop": 15.0, "event": 0, "covariates": {"treat": 1, "age": 50}},
            {"subject_id": "s1", "start": 10.0, "stop": 20.0, "event": 1, "covariates": {"treat": 0}},
            {"subject_id": "s1", "start": 0.0, "stop": 10.0, "event": 0, "covariates": {"treat": 1}},
        ]

        result = time_varying_covariates(intervals)

        assert result["n_subjects"] == 2
        assert result["n_intervals"] == 3
        assert result["interval_counts"] == {"s2": 1, "s1": 2}
        assert result["covariate_names"] == ["age", "treat"]

        rows = result["expanded_data"]
        assert [(r["subject_id"], r["start"]) for r in rows] == [("s1", 0.0), ("s1", 10.0), ("s2", 0.0)]
        assert rows[2]["age"] == 50
        assert rows[0]["event"] == 0 and rows[1]["event"] == 1
