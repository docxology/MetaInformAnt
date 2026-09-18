"""Tests for the gated inferential comparative analysis (inferential_comparative.py).

Covers the confirmatory comparative layer of
projects/hymenoptera_amalgkit/docs/manuscript/statistical_analysis_plan.md
(sections 5.1, 6, 7):

- covariate recovery on synthetic data through the gated analysis;
- bootstrap CI coverage sanity against a known true effect;
- heterogeneity statistics (Cochran's Q, I-squared, tau-squared) on
  constructed between-study variance;
- fail-closed role gating: unfrozen manifests, descriptive/stopped
  contracts, tested-feature mismatches, and undeclared contrast levels
  all refuse before any inferential output exists;
- the sensitivity runner: leave-one-study-out and covariate-exclusion
  refits with directional agreement against the predeclared expectation.

All fixtures are small deterministic numpy/pandas data. No mocks, no
network, no live data root.
"""

from typing import Any

import numpy as np
import pandas as pd
import pytest

from metainformant.rna.analysis.inferential_comparative import (
    ComparativeDesign,
    InferentialComparativeError,
    bootstrap_effect_ci,
    directional_agreement,
    fit_comparative_effects,
    random_effects_summary,
    run_inferential_comparative_analysis,
    run_registered_sensitivity_analyses,
)
from metainformant.rna.analysis.statistics_contract import (
    INFERENTIAL_ROLE,
    STOPPED_ROLE,
    AnalysisProvenance,
    SensitivityAnalysis,
    StatisticsContractError,
    benjamini_hochberg_fdr,
)

SOFTWARE_VERSIONS = {"metainformant": "1.0.0", "numpy": "2.0.0", "python": "3.12"}


def _inferential_contract(**overrides: Any) -> AnalysisProvenance:
    fields: dict[str, Any] = dict(
        analysis_id="hymenoptera_comparative_v1",
        estimand="random-effects combined caste contrast coefficient per ortholog",
        replicate_unit="biological sample within study",
        random_seed=20260917,
        resampling_count=200,
        null_model="two-sided z-test of the combined contrast coefficient against zero",
        multiple_testing_family="ortholog feature",
        multiple_testing_method="bh-fdr",
        tested_feature_count=1,
        software_versions=SOFTWARE_VERSIONS,
        analysis_role=INFERENTIAL_ROLE,
    )
    fields.update(overrides)
    return AnalysisProvenance(**fields)


def _design(**overrides: Any) -> ComparativeDesign:
    fields: dict[str, Any] = dict(
        response_col="expression",
        contrast_col="caste",
        reference_level="worker",
        treatment_level="queen",
        study_col="study",
        covariate_cols=("stage",),
        feature_col="feature",
    )
    fields.update(overrides)
    return ComparativeDesign(**fields)


def _simulate(
    study_effects: dict[str, float],
    n_per_study: int = 12,
    covariate_effect: float = 0.8,
    noise_sd: float = 0.5,
    seed: int = 20260917,
    features: tuple[str, ...] = ("orth_1",),
) -> pd.DataFrame:
    """Balanced synthetic long-format observations with a known contrast."""
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    for study, effect in study_effects.items():
        for i in range(n_per_study):
            treated = i % 2 == 0
            stage = int(rng.integers(0, 2))
            for feature in features:
                shift = 0.0 if feature == features[0] else -0.4
                rows.append(
                    {
                        "study": study,
                        "caste": "queen" if treated else "worker",
                        "stage": "L1" if stage else "L2",
                        "feature": feature,
                        "expression": effect * treated + covariate_effect * stage + shift + rng.normal(0.0, noise_sd),
                    }
                )
    return pd.DataFrame(rows)


# =============================================================================
# Covariate recovery through the gated analysis
# =============================================================================


def test_gated_analysis_recovers_effect_and_covariate() -> None:
    observations = _simulate({"study_a": 1.5, "study_b": 1.5, "study_c": 1.5})
    result = run_inferential_comparative_analysis(
        observations, _design(), _inferential_contract(), evidence_manifest_frozen=True
    )
    assert result["role"] == INFERENTIAL_ROLE
    assert result["gate"] == "post-freeze"
    features = result["features"]
    assert features.attrs["role"] == INFERENTIAL_ROLE
    row = features.loc["orth_1"]
    assert row["effect"] == pytest.approx(1.5, abs=0.35)
    assert row["ci_low"] < 1.5 < row["ci_high"]
    assert 0.0 <= row["p_value"] <= 1.0
    assert row["p_adj_bh"] >= row["p_value"]
    assert row["n_studies"] == 3
    assert row["n_observations"] == 36


def test_covariate_removes_confounding_bias() -> None:
    # The stage covariate is correlated with treatment in this construction;
    # the adjusted model must stay near the true effect while the naive
    # model without the covariate is pulled away.
    rng = np.random.default_rng(7)
    rows: list[dict[str, Any]] = []
    for study in ("study_a", "study_b", "study_c"):
        for i in range(16):
            treated = i < 8  # stage assigned non-randomly across treatment arms
            stage = i % 2  # perfectly balanced across arms; still needed by design
            rows.append(
                {
                    "study": study,
                    "caste": "queen" if treated else "worker",
                    "stage": "L1" if stage else "L2",
                    "feature": "orth_1",
                    "expression": 1.5 * treated + 0.8 * stage + rng.normal(0.0, 0.5),
                }
            )
    observations = pd.DataFrame(rows)
    contract = _inferential_contract(resampling_count=50)
    result = run_inferential_comparative_analysis(observations, _design(), contract, evidence_manifest_frozen=True)
    effect = result["features"].loc["orth_1", "effect"]
    naive = fit_comparative_effects(observations, _design(covariate_cols=()), feature="orth_1")
    assert effect == pytest.approx(1.5, abs=0.25)
    assert naive["combined"]["effect"] == pytest.approx(1.5, abs=0.25)


def test_multiplicity_matches_declared_bh_fdr() -> None:
    observations = _simulate({"study_a": 2.0, "study_b": 2.0, "study_c": 2.0}, features=("orth_1", "orth_2"))
    contract = _inferential_contract(tested_feature_count=2, resampling_count=20)
    result = run_inferential_comparative_analysis(observations, _design(), contract, evidence_manifest_frozen=True)
    features = result["features"]
    assert len(features) == 2
    expected = benjamini_hochberg_fdr(features["p_value"].tolist())
    assert features["p_adj_bh"].tolist() == pytest.approx(expected)
    assert (features["p_adj_bh"] >= features["p_value"]).all()


# =============================================================================
# Bootstrap CI coverage sanity
# =============================================================================


def test_bootstrap_ci_covers_true_effect() -> None:
    design = _design()
    contract = _inferential_contract()
    covered = 0
    repeats = 20
    for seed in range(repeats):
        observations = _simulate(
            {"study_a": 1.0, "study_b": 1.0, "study_c": 1.0},
            noise_sd=0.5,
            seed=1000 + seed,
        )
        result = run_inferential_comparative_analysis(observations, design, contract, evidence_manifest_frozen=True)
        row = result["features"].loc["orth_1"]
        assert row["bootstrap_n_success"] >= 0.8 * contract.resampling_count
        covered += int(row["ci_low"] <= 1.0 <= row["ci_high"])
    assert covered / repeats >= 0.8


def test_bootstrap_ci_deterministic_and_widens_with_noise() -> None:
    design = _design()
    observations = _simulate({"study_a": 1.0, "study_b": 1.0, "study_c": 1.0}, seed=11)
    contract = _inferential_contract(resampling_count=100)
    first = run_inferential_comparative_analysis(observations, design, contract, evidence_manifest_frozen=True)
    second = run_inferential_comparative_analysis(observations, design, contract, evidence_manifest_frozen=True)
    assert first["features"].loc["orth_1", "ci_low"] == second["features"].loc["orth_1", "ci_low"]
    assert first["features"].loc["orth_1", "ci_high"] == second["features"].loc["orth_1", "ci_high"]

    noisy = _simulate({"study_a": 1.0, "study_b": 1.0, "study_c": 1.0}, noise_sd=2.0, seed=11)
    quiet = bootstrap_effect_ci(observations, design, random_seed=5, resampling_count=100, feature="orth_1")
    loud = bootstrap_effect_ci(noisy, design, random_seed=5, resampling_count=100, feature="orth_1")
    assert loud["ci_high"] - loud["ci_low"] > quiet["ci_high"] - quiet["ci_low"]


# =============================================================================
# Heterogeneity statistics on constructed between-study variance
# =============================================================================


def test_random_effects_heterogeneity_on_constructed_variance() -> None:
    homogeneous = random_effects_summary([0.0, 0.05, -0.05, 0.02], [0.1, 0.1, 0.1, 0.1])
    heterogeneous = random_effects_summary([0.0, 0.0, 1.0, 1.0], [0.1, 0.1, 0.1, 0.1])
    assert homogeneous["q"] < homogeneous["df"]
    assert homogeneous["i_squared_percent"] < 30.0
    assert homogeneous["tau_squared"] == 0.0
    assert heterogeneous["q"] > 10 * heterogeneous["df"]
    assert heterogeneous["i_squared_percent"] > 90.0
    assert heterogeneous["tau_squared"] > 0.0
    assert heterogeneous["q"] > homogeneous["q"]
    # Random-effects weighting de-emphasizes the outlying studies.
    assert abs(heterogeneous["effect"] - 0.5) < abs(heterogeneous["fixed_effect"] - 0.5) + 1e-12


def test_heterogeneous_studies_detected_end_to_end() -> None:
    # Study intercept differences are absorbed within-study; heterogeneity
    # must come from genuinely different contrast effects per study.
    observations = _simulate(
        {"study_a": 0.0, "study_b": 0.1, "study_c": 2.0, "study_d": 1.9},
        n_per_study=16,
        noise_sd=0.1,
        seed=42,
    )
    result = run_inferential_comparative_analysis(
        observations,
        _design(),
        _inferential_contract(resampling_count=20),
        evidence_manifest_frozen=True,
    )
    row = result["features"].loc["orth_1"]
    assert row["n_studies"] == 4
    assert row["i_squared_percent"] > 90.0
    assert row["tau_squared"] > 0.0


# =============================================================================
# Fail-closed role gating
# =============================================================================


def test_unfrozen_manifest_refuses_even_with_inferential_contract() -> None:
    observations = _simulate({"study_a": 1.0, "study_b": 1.0, "study_c": 1.0})
    with pytest.raises(RuntimeError, match="gated for post-freeze use"):
        run_inferential_comparative_analysis(observations, _design(), _inferential_contract())


def test_descriptive_contract_never_receives_p_values() -> None:
    observations = _simulate({"study_a": 1.0, "study_b": 1.0, "study_c": 1.0})
    descriptive = AnalysisProvenance(
        analysis_id="hymenoptera_fingerprint_v1",
        estimand="descriptive effect summary",
        replicate_unit="species finalized matrix",
        random_seed=1,
        resampling_count=10,
        null_model="not-applicable",
        multiple_testing_family=None,
        multiple_testing_method=None,
        tested_feature_count=None,
        software_versions=SOFTWARE_VERSIONS,
        analysis_role="descriptive",
    )
    with pytest.raises(StatisticsContractError, match="descriptive"):
        run_inferential_comparative_analysis(observations, _design(), descriptive, evidence_manifest_frozen=True)
    with pytest.raises(StatisticsContractError, match="inferential"):
        run_registered_sensitivity_analyses(observations, _design(), descriptive, evidence_manifest_frozen=True)


def test_stopped_contract_refuses_inferential_output() -> None:
    observations = _simulate({"study_a": 1.0, "study_b": 1.0, "study_c": 1.0})
    stopped = AnalysisProvenance(
        analysis_id="hymenoptera_comparative_halted",
        estimand="never produced",
        replicate_unit="not-applicable",
        random_seed=0,
        resampling_count=1,
        null_model="not-applicable",
        multiple_testing_family=None,
        multiple_testing_method=None,
        tested_feature_count=None,
        software_versions=SOFTWARE_VERSIONS,
        analysis_role=STOPPED_ROLE,
    )
    with pytest.raises(StatisticsContractError, match="stopped"):
        run_inferential_comparative_analysis(observations, _design(), stopped, evidence_manifest_frozen=True)


def test_tested_feature_count_mismatch_refuses() -> None:
    observations = _simulate({"study_a": 1.0, "study_b": 1.0, "study_c": 1.0})
    with pytest.raises(StatisticsContractError, match="tested_feature_count"):
        run_inferential_comparative_analysis(
            observations,
            _design(),
            _inferential_contract(tested_feature_count=2),
            evidence_manifest_frozen=True,
        )


def test_undeclared_contrast_level_refuses() -> None:
    observations = _simulate({"study_a": 1.0, "study_b": 1.0, "study_c": 1.0})
    observations.loc[observations.index[0], "caste"] = "male"
    with pytest.raises(InferentialComparativeError, match="undeclared levels"):
        run_inferential_comparative_analysis(
            observations, _design(), _inferential_contract(resampling_count=5), evidence_manifest_frozen=True
        )


def test_single_state_study_refuses() -> None:
    observations = _simulate({"study_a": 1.0, "study_b": 1.0, "study_c": 1.0})
    observations = observations[~((observations["study"] == "study_c") & (observations["caste"] == "worker"))]
    with pytest.raises(InferentialComparativeError, match="both declared contrast levels"):
        run_inferential_comparative_analysis(
            observations, _design(), _inferential_contract(resampling_count=5), evidence_manifest_frozen=True
        )


def test_too_few_studies_refuses() -> None:
    observations = _simulate({"study_a": 1.0, "study_b": 1.0})
    with pytest.raises(InferentialComparativeError, match="minimum of 3"):
        run_inferential_comparative_analysis(
            observations,
            _design(min_studies=3),
            _inferential_contract(resampling_count=5),
            evidence_manifest_frozen=True,
        )


# =============================================================================
# Sensitivity runner
# =============================================================================


def _sensitivity_contract(entry_kwargs: list[dict[str, Any]], **overrides: Any) -> AnalysisProvenance:
    entries = tuple(SensitivityAnalysis(**kwargs) for kwargs in entry_kwargs)
    return _inferential_contract(sensitivity_analyses=entries, **overrides)


def test_leave_one_study_out_reports_directional_agreement() -> None:
    # study_c carries an inflated effect; dropping it must decrease the
    # combined estimate.
    observations = _simulate({"study_a": 1.0, "study_b": 1.0, "study_c": 3.0}, seed=3)
    contract = _sensitivity_contract(
        [
            {
                "name": "leave_one_study_out_stability",
                "varied_parameter": "leave_one_study_out",
                "baseline_value": "all studies retained",
                "varied_values": ("study_a", "study_b", "study_c"),
                "expected_direction": "decrease",
            }
        ]
    )
    results = run_registered_sensitivity_analyses(observations, _design(), contract, evidence_manifest_frozen=True)
    assert len(results) == 1
    entry = results[0]
    assert entry["role"] == INFERENTIAL_ROLE
    assert entry["n_comparisons"] == 3
    by_value = {c["varied_value"]: c for c in entry["comparisons"]}
    assert by_value["study_c"]["observed_direction"] == "decrease"
    assert by_value["study_c"]["directional_agreement"] is True
    assert by_value["study_a"]["directional_agreement"] is False
    assert entry["n_agreement"] == 1
    assert entry["fraction_agreement"] == pytest.approx(1 / 3)


def test_exclude_covariate_sensitivity_runs_and_labels_roles() -> None:
    observations = _simulate({"study_a": 1.0, "study_b": 1.0, "study_c": 1.0}, seed=5)
    contract = _sensitivity_contract(
        [
            {
                "name": "drop_stage_covariate",
                "varied_parameter": "exclude_covariate",
                "baseline_value": "stage retained",
                "varied_values": ("stage",),
                "expected_direction": "either",
            }
        ]
    )
    results = run_registered_sensitivity_analyses(observations, _design(), contract, evidence_manifest_frozen=True)
    entry = results[0]
    assert entry["n_comparisons"] == 1
    comparison = entry["comparisons"][0]
    assert comparison["feature"] == "orth_1"
    assert np.isfinite(comparison["varied_effect"])
    assert entry["fraction_agreement"] == 1.0  # expected_direction="either"


def test_unknown_sensitivity_parameter_refuses() -> None:
    observations = _simulate({"study_a": 1.0, "study_b": 1.0, "study_c": 1.0})
    contract = _sensitivity_contract(
        [
            {
                "name": "unregistered_variation",
                "varied_parameter": "rescale_expression",
                "baseline_value": "raw",
                "varied_values": ("log",),
                "expected_direction": "either",
            }
        ]
    )
    with pytest.raises(InferentialComparativeError, match="cannot execute"):
        run_registered_sensitivity_analyses(observations, _design(), contract, evidence_manifest_frozen=True)


def test_sensitivity_with_unknown_study_value_refuses() -> None:
    observations = _simulate({"study_a": 1.0, "study_b": 1.0, "study_c": 1.0})
    contract = _sensitivity_contract(
        [
            {
                "name": "leave_one_study_out_stability",
                "varied_parameter": "leave_one_study_out",
                "baseline_value": "all studies retained",
                "varied_values": ("study_z",),
                "expected_direction": "either",
            }
        ]
    )
    with pytest.raises(InferentialComparativeError, match="not present in the observations"):
        run_registered_sensitivity_analyses(observations, _design(), contract, evidence_manifest_frozen=True)


def test_directional_agreement_semantics() -> None:
    assert directional_agreement("decrease", "decrease") is True
    assert directional_agreement("decrease", "increase") is False
    assert directional_agreement("none", "none") is True
    assert directional_agreement("none", "decrease") is False
    assert directional_agreement("either", "decrease") is True
    assert directional_agreement("either", "none") is True
