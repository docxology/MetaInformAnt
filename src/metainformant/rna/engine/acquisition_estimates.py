"""Transparent cost and elapsed-time models for configurable acquisition lanes."""
from __future__ import annotations

import math
from dataclasses import dataclass


class AcquisitionEstimateError(ValueError):
    """Missing or inconsistent evidence prevents a meaningful estimate."""

    def __init__(self, field: str, reason: str) -> None:
        self.field, self.reason = field, reason
        super().__init__(f"{field}: {reason}")


def _nonnegative(field: str, value: float) -> None:
    if isinstance(value, bool) or not math.isfinite(value) or value < 0:
        raise AcquisitionEstimateError(field, "must be finite and nonnegative")


@dataclass(frozen=True, slots=True)
class ThroughputEvidence:
    """Observed whole-lane rates at a declared capacity, not a scaling guarantee.

    Units are VM instances for AWS or sample-worker slots for a local process.
    Rate bounds are operator-supplied observed scenarios, not confidence limits.
    """

    observed_units: int
    samples_per_hour_low: float
    samples_per_hour_high: float
    source: str
    fleet_rate_cap: float | None = None

    def __post_init__(self) -> None:
        if type(self.observed_units) is not int or self.observed_units <= 0:
            raise AcquisitionEstimateError("observed_units", "must be a positive integer")
        for field, value in (("rate_low", self.samples_per_hour_low), ("rate_high", self.samples_per_hour_high)):
            _nonnegative(field, value)
            if value == 0:
                raise AcquisitionEstimateError(field, "zero throughput cannot support an ETA")
        if self.samples_per_hour_high < self.samples_per_hour_low or not self.source.strip():
            raise AcquisitionEstimateError("evidence", "ordered rate bounds and a source are required")
        if self.fleet_rate_cap is not None:
            _nonnegative("fleet_rate_cap", self.fleet_rate_cap)
            if self.fleet_rate_cap == 0:
                raise AcquisitionEstimateError("fleet_rate_cap", "must be positive")


@dataclass(frozen=True, slots=True)
class LaneCosts:
    """Gross costs; credits are deliberately excluded.

    AWS hourly quotes should include compute, attached EBS, public IP and margin.
    Storage, requests and network charges outside that quote belong in fixed_usd.
    A zero local rate means an explicit marginal-cost assumption, not free compute.
    """

    hourly_usd_per_unit: float
    fixed_usd: float = 0
    setup_hours: float = 0
    retries_per_task: float = 0

    def __post_init__(self) -> None:
        for name in ("hourly_usd_per_unit", "fixed_usd", "setup_hours", "retries_per_task"):
            _nonnegative(name, getattr(self, name))


@dataclass(frozen=True, slots=True)
class LaneEstimate:
    samples: int
    units: int
    rate_low: float
    rate_high: float
    hours_low: float
    hours_high: float
    cost_low_usd: float
    cost_high_usd: float
    extrapolated: bool
    evidence_source: str
    assumptions: tuple[str, ...]


def estimate_lane(samples: int, units: int, evidence: ThroughputEvidence, costs: LaneCosts) -> LaneEstimate:
    """Bound aggregate completion scenarios with an optional saturation cap."""
    if type(samples) is not int or samples < 0 or type(units) is not int or units <= 0:
        raise AcquisitionEstimateError("capacity", "samples must be nonnegative and units positive integers")
    factor = units / evidence.observed_units
    low, high = evidence.samples_per_hour_low * factor, evidence.samples_per_hour_high * factor
    if evidence.fleet_rate_cap is not None:
        low, high = min(low, evidence.fleet_rate_cap), min(high, evidence.fleet_rate_cap)
    tasks = samples * (1 + costs.retries_per_task)
    hours_low = costs.setup_hours + tasks / high if samples else 0.0
    hours_high = costs.setup_hours + tasks / low if samples else 0.0
    hourly = units * costs.hourly_usd_per_unit
    fixed = costs.fixed_usd if samples else 0.0
    return LaneEstimate(samples, units, low, high, hours_low, hours_high,
                        fixed + hourly * hours_low, fixed + hourly * hours_high,
                        units != evidence.observed_units, evidence.source,
                        ("Rates assume a comparable sample-size and source mix.",
                         "Scaling to a different capacity is extrapolation; benchmark before admission.",
                         "Setup, retries and a shared rate cap can increase cost when parallelism rises.",
                         "All configured units are charged for the modeled lane duration; idle shutdown can reduce cost."))


@dataclass(frozen=True, slots=True)
class CampaignEstimate:
    hours_low: float
    hours_high: float
    gross_low_usd: float
    gross_high_usd: float
    ceiling_usd: float | None
    fits_conservative_ceiling: bool | None


def combine_estimates(lanes: tuple[LaneEstimate, ...], *, spent_usd: float = 0,
                      reserved_usd: float = 0, ceiling_usd: float | None = None) -> CampaignEstimate:
    """Concurrent lanes finish at their slowest lane; gross cost is additive.

    reserved_usd is committed work outside the supplied pending-lane estimates.
    Do not include the same active worker runtime in both terms.
    """
    _nonnegative("spent_usd", spent_usd)
    _nonnegative("reserved_usd", reserved_usd)
    if not lanes:
        raise AcquisitionEstimateError("lanes", "at least one modeled lane is required")
    if ceiling_usd is not None:
        _nonnegative("ceiling_usd", ceiling_usd)
    low = spent_usd + reserved_usd + math.fsum(lane.cost_low_usd for lane in lanes)
    high = spent_usd + reserved_usd + math.fsum(lane.cost_high_usd for lane in lanes)
    return CampaignEstimate(max(lane.hours_low for lane in lanes), max(lane.hours_high for lane in lanes), low, high,
                            ceiling_usd, high <= ceiling_usd if ceiling_usd is not None else None)
