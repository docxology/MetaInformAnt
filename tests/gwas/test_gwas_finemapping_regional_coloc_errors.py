"""Regression tests: ``regional_coloc`` must not report success when no window could be analysed.

Defect: every window of a region is analysed by ``multi_trait_coloc``, which returns ``{"status": "error", ...}``
for invalid per-window input (a non-positive standard error, a prior variance named for an unknown trait, a
non-finite z-score). ``regional_coloc`` kept those error dictionaries as ordinary window results, picked the
"best" window among them (an error window has no ``multi_trait_pp_shared``, which was read as 0.0 and beat the
initial -1.0) and returned top-level ``status == "success"`` with an error dictionary as ``best_window``.

The invariants checked here do not depend on how the PP(shared) values are computed: an error status when every
window failed, a ``best_window`` that is a successful window with the largest PP(shared) among the windows that
succeeded, and a failure count that matches the windows carrying an error status.
"""

from __future__ import annotations

import pytest

from metainformant.gwas.finemapping.colocalization import regional_coloc

REGION = {"chrom": "1", "start": 0, "end": 5000}

# Six variants: one window (a region of at most 100 variants is analysed as a single window).
Z_SMALL = {
    "height": [0.5, 1.0, 5.0, 0.2, 0.1, 0.3],
    "bmi": [0.1, 0.9, 4.0, 0.3, 0.2, 0.1],
}

N_LARGE = 150  # more than 100 variants: windows of 100 variants advancing by 50 -> (0, 100) and (50, 150)


def _large_z() -> dict[str, list[float]]:
    """150 variants with a shared signal at index 120 and a second, trait-specific one at index 30."""
    height = [0.3 * ((i * 7) % 5 - 2) for i in range(N_LARGE)]
    bmi = [0.3 * ((i * 3) % 5 - 2) for i in range(N_LARGE)]
    height[120], bmi[120] = 9.0, 8.0
    height[30] = 6.0
    return {"height": height, "bmi": bmi}


@pytest.mark.parametrize(
    ("kwargs", "expected_reason"),
    [
        ({"prior_variances": {"no_such_trait": 0.1}}, "unknown traits"),
        ({"standard_errors": {"height": -1.0, "bmi": 0.1}}, "standard_errors must be positive"),
        ({"standard_errors": {"height": 0.0, "bmi": 0.1}}, "standard_errors must be positive"),
    ],
)
def test_every_window_failing_is_an_error_not_a_success(kwargs: dict, expected_reason: str) -> None:
    result = regional_coloc(Z_SMALL, REGION, **kwargs)

    assert result["status"] == "error"
    assert expected_reason in result["message"]
    assert "all 1 windows" in result["message"]
    assert "best_window" not in result
    # The per-window errors stay available for diagnosis.
    assert result["n_windows"] == 1
    assert result["n_failed_windows"] == 1
    assert [window["status"] for window in result["all_windows"]] == ["error"]
    assert result["region"]["n_variants"] == 6


def test_a_non_finite_z_score_is_reported_as_an_error() -> None:
    z = {"height": [0.5, float("nan"), 5.0], "bmi": [0.1, 0.9, 4.0]}

    result = regional_coloc(z, REGION)

    assert result["status"] == "error"
    assert "finite" in result["message"]


def test_failures_in_every_window_of_a_multi_window_region_are_reported_once_per_reason() -> None:
    result = regional_coloc(_large_z(), REGION, prior_variances={"no_such_trait": 0.1})

    assert result["status"] == "error"
    assert result["n_windows"] == 2
    assert result["n_failed_windows"] == 2
    assert result["message"].count("unknown traits") == 1  # two windows, one distinct reason
    assert "all 2 windows" in result["message"]


def test_a_failed_window_is_never_the_best_window() -> None:
    """A bad standard error at variant 10 breaks window (0, 100) only; window (50, 150) still succeeds."""
    standard_errors = [0.1] * N_LARGE
    standard_errors[10] = -1.0

    result = regional_coloc(_large_z(), REGION, standard_errors={"height": standard_errors, "bmi": 0.1})

    assert result["status"] == "success"
    assert result["n_windows"] == 2
    assert result["n_failed_windows"] == 1
    statuses = {window["window_index"]: window["status"] for window in result["all_windows"]}
    assert statuses == {0: "error", 1: "success"}
    assert result["best_window"]["status"] == "success"
    assert result["best_window"]["window_index"] == 1
    # The failed window keeps its own message instead of being dropped.
    failed = next(window for window in result["all_windows"] if window["status"] == "error")
    assert "standard_errors must be positive" in failed["message"]


def test_best_window_is_the_successful_window_with_the_largest_pp_shared() -> None:
    result = regional_coloc(_large_z(), REGION)

    assert result["status"] == "success"
    assert result["n_failed_windows"] == 0
    assert result["n_windows"] == 2
    pp_values = [window["multi_trait_pp_shared"] for window in result["all_windows"]]
    assert result["best_window"]["multi_trait_pp_shared"] == max(pp_values)
    # The shared signal at variant 120 lies in both windows; the trait-specific one at 30 only in the first.
    # Each window therefore sees a different mixture of evidence, so the two values must differ.
    assert pp_values[0] != pp_values[1]


def test_a_single_window_region_keeps_its_result() -> None:
    result = regional_coloc(Z_SMALL, REGION)

    assert result["status"] == "success"
    assert result["n_windows"] == 1
    assert result["n_failed_windows"] == 0
    assert result["best_window"] is result["all_windows"][0]
    assert 0.0 <= result["best_window"]["multi_trait_pp_shared"] <= 1.0
    assert result["region"] == {"chrom": "1", "start": 0, "end": 5000, "n_variants": 6}


def test_input_validation_errors_are_unchanged() -> None:
    assert regional_coloc({}, REGION)["status"] == "error"
    assert regional_coloc({"only_one": [0.1, 0.2]}, REGION)["status"] == "error"
    assert regional_coloc({"a": [0.1, 0.2], "b": [0.1]}, REGION)["status"] == "error"
    assert regional_coloc({"a": [], "b": []}, REGION)["status"] == "error"
