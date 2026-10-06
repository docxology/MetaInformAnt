"""Independent cost arithmetic and historical checkpoint rate controls."""

import json
import math
import pytest
from metainformant.rna.engine.aws_resources import (
    WorkerPrices,
    WorkerPricingError,
    catalog_unit_price,
    runtime_charge,
)
from metainformant.rna.engine.aws_completion import budget_allows


def test_large_disk_reservation_rejects_old_underestimate() -> None:
    profile = WorkerPrices(0.357, 0.08, 0.55)
    rate = profile.hourly_bound(2000)
    expected = 0.357 + 2000 * 0.08 / (28 * 24) + 0.005 + 0.05
    assert rate == pytest.approx(expected)
    assert budget_allows(733.3, 750, 43200, 0.55)
    assert not budget_allows(733.3, 750, 43200, rate)


def test_current_disk_keeps_the_configured_floor() -> None:
    assert WorkerPrices(0.357, 0.08, 0.55).hourly_bound(600) == 0.55


def test_each_runtime_keeps_its_original_rate() -> None:
    assert runtime_charge(0, 3600, 0.55) + runtime_charge(
        3600, 10800, 0.8
    ) == pytest.approx(2.15)
    assert runtime_charge(0, 3660, 0.55) > runtime_charge(0, 3600, 0.55)


@pytest.mark.parametrize("bad", [math.nan, math.inf, -1, True])
def test_invalid_resource_prices_fail_closed(bad: float) -> None:
    with pytest.raises(WorkerPricingError):
        WorkerPrices(0.357, bad, 0.55)


@pytest.mark.parametrize("disk", [0, -1, True, 600.5])
def test_invalid_disk_profiles_fail_closed(disk: int) -> None:
    with pytest.raises(WorkerPricingError):
        WorkerPrices(0.357, 0.08, 0.55).hourly_bound(disk)


def test_catalog_selection_respects_units_and_cardinality() -> None:
    product = json.dumps(
        {
            "terms": {
                "OnDemand": {
                    "term": {
                        "priceDimensions": {
                            "volume": {
                                "unit": "GB-Mo",
                                "pricePerUnit": {"USD": "0.08"},
                            },
                            "compute": {
                                "unit": "Hrs",
                                "pricePerUnit": {"USD": "0.357"},
                            },
                        }
                    }
                }
            }
        }
    )
    assert catalog_unit_price([product], "GB-Mo") == 0.08
    assert catalog_unit_price([product], "Hrs") == 0.357
    with pytest.raises(WorkerPricingError):
        catalog_unit_price([product, product], "Hrs")
    with pytest.raises(WorkerPricingError):
        catalog_unit_price([product], "unknown")


def test_mixed_checkpoint_preserves_legacy_and_new_rates() -> None:
    from metainformant.rna.engine.aws_resources import campaign_charge

    state = {
        "historical_gross": 67.92,
        "hourly_upper_bound": 0.55,
        "jobs": [
            {"started_at": 0.0, "finished_at": 3600.0},
            {"started_at": 3600.0, "hourly_upper_bound": 0.8},
        ],
    }
    assert campaign_charge(state, 10800.0) == pytest.approx(70.07)
    assert campaign_charge(state, 14400.0) == pytest.approx(70.87)


def test_terminating_checkpoint_keeps_charging_until_observed_finish() -> None:
    from metainformant.rna.engine.aws_resources import campaign_charge

    state = {
        "historical_gross": 0.0,
        "hourly_upper_bound": 0.55,
        "jobs": [{"started_at": 0.0}],
    }
    assert campaign_charge(state, 3660.0) == pytest.approx(0.55 * 61 / 60)
    state["jobs"][0]["finished_at"] = 3660.0
    assert campaign_charge(state, 7200.0) == pytest.approx(0.55 * 61 / 60)
