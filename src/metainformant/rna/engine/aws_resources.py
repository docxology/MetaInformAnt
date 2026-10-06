"""Conservative worker resource prices and elapsed-runtime accounting."""

from __future__ import annotations
from dataclasses import dataclass
import json
import math
from typing import Final, Sequence, TypedDict, NotRequired

MIN_MONTH_HOURS: Final = 28 * 24
PUBLIC_IPV4_HOURLY: Final = 0.005
HOURLY_MARGIN: Final = 0.05


class WorkerPricingError(ValueError):
    """Resource pricing or persisted runtime cannot produce a safe bound."""

    def __init__(self, field: str, reason: str) -> None:
        self.field = field
        self.reason = reason
        super().__init__(f"{field}: {reason}")


def catalog_unit_price(products: Sequence[str], unit: str) -> float:
    """Extract exactly one positive finite rate in the requested catalog unit."""
    rates = []
    for encoded in products:
        product = json.loads(encoded)
        for term in product["terms"]["OnDemand"].values():
            for dimension in term["priceDimensions"].values():
                if dimension["unit"] == unit:
                    value = dimension["pricePerUnit"]["USD"]
                    if isinstance(value, bool):
                        raise WorkerPricingError(unit, "invalid boolean rate")
                    rates.append(float(value))
    if len(rates) != 1 or not math.isfinite(rates[0]) or rates[0] <= 0:
        raise WorkerPricingError(unit, "catalog must supply one positive finite rate")
    return rates[0]


@dataclass(frozen=True, slots=True)
class WorkerPrices:
    """Current compute/storage quotes and the operator's hourly floor."""

    compute_hourly: float
    gp3_gib_month: float
    hourly_floor: float

    def __post_init__(self) -> None:
        for name, value in [
            ("compute_hourly", self.compute_hourly),
            ("gp3_gib_month", self.gp3_gib_month),
            ("hourly_floor", self.hourly_floor),
        ]:
            if isinstance(value, bool) or not math.isfinite(value) or value <= 0:
                raise WorkerPricingError(name, "must be finite and positive")

    def hourly_bound(self, disk_gib: int) -> float:
        """Include baseline gp3, one public IPv4 address and an operating margin."""
        if type(disk_gib) is not int or disk_gib <= 0:
            raise WorkerPricingError("disk_gib", "must be a positive integer")
        return max(
            self.hourly_floor,
            self.compute_hourly
            + self.gp3_gib_month * disk_gib / MIN_MONTH_HOURS
            + PUBLIC_IPV4_HOURLY
            + HOURLY_MARGIN,
        )


def runtime_charge(started_at: float, ended_at: float, hourly_bound: float) -> float:
    """Charge an observed interval using that job's immutable admitted bound."""
    for name, value in [
        ("started_at", started_at),
        ("ended_at", ended_at),
        ("hourly_bound", hourly_bound),
    ]:
        if isinstance(value, bool) or not math.isfinite(value) or value < 0:
            raise WorkerPricingError(name, "must be finite and nonnegative")
    if ended_at < started_at:
        raise WorkerPricingError("ended_at", "precedes start")
    return (ended_at - started_at) / 3600 * hourly_bound


class AccountedJob(TypedDict):
    """Billing fields retained in an admitted job checkpoint."""

    started_at: float
    finished_at: NotRequired[float]
    hourly_upper_bound: NotRequired[float]


class CampaignBilling(TypedDict):
    """Historical ledger rate remains the fallback for legacy jobs only."""

    historical_gross: float
    hourly_upper_bound: float
    jobs: list[AccountedJob]


def campaign_charge(state: CampaignBilling, now: float) -> float:
    """Sum observed job charges without repricing earlier admission decisions."""
    historical = state["historical_gross"]
    if not math.isfinite(historical) or historical < 0:
        raise WorkerPricingError("historical_gross", "must be finite and nonnegative")
    charges = [
        runtime_charge(
            job["started_at"],
            job.get("finished_at", max(now, job["started_at"])),
            job.get("hourly_upper_bound", state["hourly_upper_bound"]),
        )
        for job in state["jobs"]
    ]
    return historical + math.fsum(charges)
