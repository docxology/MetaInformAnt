"""Independent arithmetic controls for acquisition cost and time scenarios."""

from __future__ import annotations
import math
import pytest
from metainformant.rna.engine.acquisition_estimates import (
    AcquisitionEstimateError,
    ThroughputEvidence,
    LaneCosts,
    estimate_lane,
    combine_estimates,
)


def test_linear_capacity_preserves_work_cost_without_overheads() -> None:
    evidence = ThroughputEvidence(2, 20, 40, "measured two-worker window")
    baseline = estimate_lane(200, 2, evidence, LaneCosts(0.5))
    larger = estimate_lane(200, 4, evidence, LaneCosts(0.5))
    assert baseline.hours_low == 5
    assert baseline.hours_high == 10
    assert larger.hours_high == 5
    assert larger.cost_high_usd == baseline.cost_high_usd == 10
    assert larger.extrapolated and not baseline.extrapolated


def test_saturation_and_setup_break_the_linear_cost_assumption() -> None:
    evidence = ThroughputEvidence(2, 40, 40, "measured", fleet_rate_cap=40)
    baseline = estimate_lane(80, 2, evidence, LaneCosts(1, setup_hours=0.5))
    larger = estimate_lane(80, 4, evidence, LaneCosts(1, setup_hours=0.5))
    assert baseline.hours_high == larger.hours_high == 2.5
    assert larger.cost_high_usd == 2 * baseline.cost_high_usd


def test_hybrid_duration_is_max_cost_is_sum_and_credits_never_enter() -> None:
    a = estimate_lane(80, 2, ThroughputEvidence(2, 20, 40, "local"), LaneCosts(0))
    b = estimate_lane(
        120, 3, ThroughputEvidence(3, 30, 60, "aws"), LaneCosts(1, fixed_usd=2)
    )
    result = combine_estimates((a, b), spent_usd=10, reserved_usd=3, ceiling_usd=26)
    assert result.hours_low == 2
    assert result.hours_high == 4
    assert result.gross_high_usd == 27
    assert result.fits_conservative_ceiling is False


def test_retry_scenario_and_no_pending_work() -> None:
    evidence = ThroughputEvidence(1, 10, 10, "observed")
    retried = estimate_lane(
        10, 1, evidence, LaneCosts(2, fixed_usd=3, retries_per_task=0.5)
    )
    assert retried.hours_high == 1.5 and retried.cost_high_usd == 6
    empty = estimate_lane(0, 1, evidence, LaneCosts(2, fixed_usd=3, setup_hours=1))
    assert empty.hours_high == empty.cost_high_usd == 0


@pytest.mark.parametrize("value", [0, -1, math.nan, math.inf, True])
def test_invalid_rate_evidence_never_produces_eta(value: float) -> None:
    with pytest.raises(AcquisitionEstimateError):
        ThroughputEvidence(1, value, 10, "observed")


@pytest.mark.parametrize("value", [-1, math.nan, math.inf, True])
def test_invalid_cost_or_budget_fails_closed(value: float) -> None:
    with pytest.raises(AcquisitionEstimateError):
        LaneCosts(value)
    lane = estimate_lane(1, 1, ThroughputEvidence(1, 1, 1, "source"), LaneCosts(1))
    with pytest.raises(AcquisitionEstimateError):
        combine_estimates((lane,), spent_usd=value)


def test_missing_source_and_reversed_bounds_refused() -> None:
    with pytest.raises(AcquisitionEstimateError):
        ThroughputEvidence(1, 10, 1, "source")
    with pytest.raises(AcquisitionEstimateError):
        ThroughputEvidence(1, 1, 10, "")
