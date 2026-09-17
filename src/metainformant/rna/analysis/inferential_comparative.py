"""Inferential comparative analysis over study-aware expression responses.

Implements the confirmatory comparative layer (layer 3) of
``projects/hymenoptera_amalgkit/docs/manuscript/statistical_analysis_plan.md``
(sections 5.1, 6, and 7):

- Study-aware comparative model: for each feature (ortholog/pathway member)
  the contrast effect is estimated *within* every study by closed-form
  ordinary least squares with the declared biological contrast as the focal
  term and tissue/caste/sex/stage-style covariates retained as available.
  Libraries are never stacked across studies as one experiment.
- Study-level random-effect approximation: the study-specific effects are
  combined with the DerSimonian-Laird method-of-moments random-effects
  estimator (no likelihood machinery; deterministic and closed-form).
- Effect sizes with uncertainty: standard errors from the within-study
  fits and the between-study combination, plus percentile bootstrap
  confidence intervals from a seeded ``numpy`` generator (resampling
  biological observations within each study).
- Heterogeneity statistics: Cochran's Q, its degrees of freedom, I-squared
  (percent), and the between-study variance tau-squared.
- Multiplicity control: the raw per-feature p-values of the combined
  z-test are adjusted only through
  :func:`metainformant.rna.analysis.statistics_contract.declared_inferential_bh_fdr`,
  which re-validates the contract and its tested-feature family.
- Sensitivity runner: executes the :class:`SensitivityAnalysis` entries
  registered on the contract (leave-one-study-out, covariate exclusion)
  and reports directional agreement against the primary estimand.

Fail-closed boundary: the entire inferential path is GATED for post-freeze
use. :func:`run_inferential_comparative_analysis` and
:func:`run_registered_sensitivity_analyses` refuse to run without an
explicit ``evidence_manifest_frozen=True`` affirmation and a validated
:class:`metainformant.rna.analysis.statistics_contract.AnalysisProvenance`
declared ``analysis_role="inferential"`` with a BH-FDR procedure. A
descriptive, stopped, or unavailable contract never produces p-values from
this module: the call raises before any output exists. Every result object
carries an explicit ``role`` label; no helper in this module emits p-values.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from metainformant.rna.analysis.statistics_contract import (
    INFERENTIAL_ROLE,
    AnalysisProvenance,
    SensitivityAnalysis,
    StatisticsContractError,
    declared_inferential_bh_fdr,
    validate_analysis_provenance,
    validate_sensitivity_analysis,
)

__all__ = [
    "ComparativeDesign",
    "InferentialComparativeError",
    "bootstrap_effect_ci",
    "directional_agreement",
    "fit_comparative_effects",
    "fit_study_effects",
    "random_effects_summary",
    "require_inferential_contract",
    "run_inferential_comparative_analysis",
    "run_registered_sensitivity_analyses",
]


class InferentialComparativeError(StatisticsContractError):
    """A comparative model cannot be estimated as declared (fail-closed)."""


@dataclass(frozen=True)
class ComparativeDesign:
    """Predeclared design of one comparative contrast.

    The design is frozen like the provenance record: it names the response
    column, the focal biological contrast (reference vs treatment level),
    the study column, optional covariates (tissue, caste, sex, stage, ...),
    and the optional feature column that splits the observations into the
    tested family (orthologs or pathways).
    """

    response_col: str
    contrast_col: str
    reference_level: str
    treatment_level: str
    study_col: str = "study"
    covariate_cols: tuple[str, ...] = ()
    feature_col: str | None = None
    min_studies: int = 2
    min_observations_per_study: int = 4


# =============================================================================
# Fail-closed gating (mirrors declared_inferential_bh_fdr)
# =============================================================================


def require_inferential_contract(
    contract: AnalysisProvenance,
    evidence_manifest_frozen: bool = False,
) -> None:
    """GATE: refuse any inferential output unless the contract allows it.

    Raises:
        RuntimeError: If ``evidence_manifest_frozen`` is not True.
        ProvenanceError: If the contract fails
            :func:`validate_analysis_provenance`.
        StatisticsContractError: If the contract is not declared
            ``analysis_role="inferential"`` with a BH-FDR procedure.
    """
    if not evidence_manifest_frozen:
        raise RuntimeError(
            "the inferential comparative analysis is gated for post-freeze use: "
            "refusing to run while the evidence manifest is unfrozen"
        )
    validate_analysis_provenance(contract)
    if contract.analysis_role != INFERENTIAL_ROLE:
        raise StatisticsContractError(
            "the inferential comparative analysis requires a contract declared "
            f"analysis_role={INFERENTIAL_ROLE!r}, got {contract.analysis_role!r}; "
            "a descriptive contract never receives p-values from this module"
        )
    method = contract.multiple_testing_method
    if method is None or str(method).strip().lower() not in {"bh-fdr", "benjamini-hochberg"}:
        raise StatisticsContractError(
            "the inferential comparative analysis requires a declared BH-FDR "
            f"procedure, got {contract.multiple_testing_method!r}"
        )


# =============================================================================
# Validation helpers
# =============================================================================


def _validate_design(design: ComparativeDesign) -> None:
    if not isinstance(design, ComparativeDesign):
        raise TypeError("design must be a ComparativeDesign record")
    for field in ("response_col", "contrast_col", "reference_level", "treatment_level", "study_col"):
        value = getattr(design, field)
        if not isinstance(value, str) or not value:
            raise ValueError(f"design field {field!r} must be a non-empty string, got {value!r}")
    if design.reference_level == design.treatment_level:
        raise ValueError("design reference_level and treatment_level must differ")
    if not isinstance(design.covariate_cols, tuple) or any(
        not isinstance(name, str) or not name for name in design.covariate_cols
    ):
        raise ValueError("design covariate_cols must be a tuple of non-empty column names")
    if design.feature_col is not None and (not isinstance(design.feature_col, str) or not design.feature_col):
        raise ValueError("design feature_col must be None or a non-empty column name")
    if not isinstance(design.min_studies, int) or isinstance(design.min_studies, bool) or design.min_studies < 2:
        raise ValueError("design min_studies must be an integer >= 2 (heterogeneity needs >= 2 studies)")
    if (
        not isinstance(design.min_observations_per_study, int)
        or isinstance(design.min_observations_per_study, bool)
        or design.min_observations_per_study < 1
    ):
        raise ValueError("design min_observations_per_study must be a positive integer")


def _validate_observations(observations: pd.DataFrame, design: ComparativeDesign) -> None:
    if not isinstance(observations, pd.DataFrame):
        raise TypeError("observations must be a pandas DataFrame in long format")
    required = [design.response_col, design.contrast_col, design.study_col, *design.covariate_cols]
    if design.feature_col is not None:
        required.append(design.feature_col)
    missing = [column for column in required if column not in observations.columns]
    if missing:
        raise ValueError(f"observations is missing required columns: {missing}")
    used = observations[required]
    if not pd.api.types.is_numeric_dtype(observations[design.response_col]):
        raise ValueError(f"response column {design.response_col!r} must be numeric")
    if used[required].isna().any().any():
        raise ValueError(
            "observations contain missing values in the columns used by the design "
            "(fail-closed: missingness must be resolved by the caller, not silently dropped)"
        )
    response = observations[design.response_col].to_numpy(dtype=float)
    if not np.isfinite(response).all():
        raise ValueError(f"response column {design.response_col!r} contains non-finite values")


def _feature_labels(observations: pd.DataFrame, design: ComparativeDesign) -> list[Any]:
    if design.feature_col is None:
        return [None]
    return sorted(observations[design.feature_col].unique(), key=str)


def _feature_frame(observations: pd.DataFrame, design: ComparativeDesign, feature: Any) -> pd.DataFrame:
    frame = observations
    if design.feature_col is not None:
        if feature is None:
            raise ValueError("feature must be provided when design.feature_col is set")
        frame = frame[frame[design.feature_col] == feature]
        if frame.empty:
            raise InferentialComparativeError(f"feature {feature!r} has no observations")
    return frame


def _study_arrays(
    sub: pd.DataFrame,
    design: ComparativeDesign,
    exclude_covariates: frozenset[str],
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Build the within-study design matrix and response as numpy arrays.

    Columns: intercept, treatment indicator, then covariates (numeric kept
    as-is; categorical expanded to drop-first indicators over sorted
    levels). No rank check here; :func:`_ols_effect` fails closed on
    rank-deficient matrices so bootstrap resamples share the same check.

    Raises:
        InferentialComparativeError: On a contrast level outside the
            declared pair or a missing contrast level.
    """
    contrast = sub[design.contrast_col]
    unexpected = set(contrast.unique()) - {design.reference_level, design.treatment_level}
    if unexpected:
        raise InferentialComparativeError(
            f"contrast column {design.contrast_col!r} carries undeclared levels "
            f"{sorted(unexpected, key=str)}; ambiguous labels must be excluded by the caller"
        )
    if design.reference_level not in set(contrast.unique()) or design.treatment_level not in set(contrast.unique()):
        raise InferentialComparativeError(
            f"study does not contain both declared contrast levels "
            f"({design.reference_level!r}, {design.treatment_level!r}); a one-state "
            "study cannot support the contrast"
        )

    columns = [np.ones(len(sub), dtype=float), (contrast == design.treatment_level).to_numpy(dtype=float)]
    names = ["intercept", "treatment"]
    for covariate in design.covariate_cols:
        if covariate in exclude_covariates:
            continue
        series = sub[covariate]
        if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series):
            columns.append(series.to_numpy(dtype=float))
            names.append(covariate)
        else:
            labels = series.astype(str).to_numpy()
            for category in sorted(set(labels))[1:]:
                columns.append((labels == category).astype(float))
                names.append(f"{covariate}[{category}]")
    return np.column_stack(columns), sub[design.response_col].to_numpy(dtype=float), names


def _ols_effect(matrix: np.ndarray, response: np.ndarray, names: list[str]) -> dict[str, Any]:
    """Closed-form OLS contrast fit on prebuilt arrays (fail-closed rank).

    Returns the treatment coefficient, its standard error, and fit sizes.
    No p-value is produced here; p-values exist only behind the gate.
    """
    n_obs, n_par = matrix.shape
    if n_obs <= n_par:
        raise InferentialComparativeError(
            f"fit has {n_obs} observations for {n_par} design columns; the model is not estimable"
        )
    if np.linalg.matrix_rank(matrix) < n_par:
        raise InferentialComparativeError(
            f"design matrix is rank-deficient (columns: {names}); a covariate is "
            "perfectly confounded with the contrast or constant within this study "
            "(plan section 8 stopping rule)"
        )
    beta, *_ = np.linalg.lstsq(matrix, response, rcond=None)
    residual = response - matrix @ beta
    sigma_squared = float(residual @ residual) / (n_obs - n_par)
    xtx_inverse = np.linalg.pinv(matrix.T @ matrix)
    standard_error = float(np.sqrt(max(sigma_squared * xtx_inverse[1, 1], 0.0)))
    return {
        "effect": float(beta[1]),
        "se": standard_error,
        "df": int(n_obs - n_par),
        "n_observations": int(n_obs),
        "n_parameters": int(n_par),
    }


def _fit_study(sub: pd.DataFrame, design: ComparativeDesign, exclude_covariates: frozenset[str]) -> dict[str, Any]:
    """Closed-form OLS contrast fit within one study (no p-values)."""
    n = len(sub)
    if n < design.min_observations_per_study:
        raise InferentialComparativeError(
            f"study has {n} observations, below the declared minimum of "
            f"{design.min_observations_per_study}"
        )
    matrix, response, names = _study_arrays(sub, design, exclude_covariates)
    return _ols_effect(matrix, response, names)


# =============================================================================
# Meta-analytic combination (study-level random-effect approximation)
# =============================================================================


def random_effects_summary(effects: Sequence[float], standard_errors: Sequence[float]) -> dict[str, Any]:
    """Combine study-specific effects with the DerSimonian-Laird estimator.

    This is the study-level random-effect approximation: a method-of-moments
    between-study variance (tau-squared) estimated from Cochran's Q, with
    random-effects weights ``1 / (se_i^2 + tau^2)``. Closed-form and
    deterministic; no likelihood optimization. Emits NO p-value.

    Args:
        effects: Per-study effect estimates (at least 2).
        standard_errors: Aligned per-study standard errors, all finite > 0.

    Returns:
        Dict with ``effect``, ``se`` (random-effects combined), the
        fixed-effect estimate, ``q`` (Cochran's Q), ``df``, ``i_squared_percent``
        (0-100), ``tau_squared``, and ``n_studies``.

    Raises:
        ValueError: On empty/mismatched input, fewer than 2 studies, or
            non-positive/non-finite standard errors.
    """
    effect_values = np.asarray(effects, dtype=float)
    se_values = np.asarray(standard_errors, dtype=float)
    if effect_values.ndim != 1 or se_values.ndim != 1 or effect_values.size != se_values.size:
        raise ValueError("effects and standard_errors must be aligned one-dimensional sequences")
    if effect_values.size < 2:
        raise ValueError(f"random-effects combination requires at least 2 studies, got {effect_values.size}")
    if not np.isfinite(effect_values).all() or not np.isfinite(se_values).all():
        raise ValueError("effects and standard_errors must be finite")
    if np.any(se_values <= 0):
        raise ValueError("standard_errors must be strictly positive")

    weights = 1.0 / se_values**2
    fixed_effect = float(np.sum(weights * effect_values) / np.sum(weights))
    q_statistic = float(np.sum(weights * (effect_values - fixed_effect) ** 2))
    df = int(effect_values.size - 1)
    c_statistic = float(np.sum(weights) - np.sum(weights**2) / np.sum(weights))
    tau_squared = max(0.0, (q_statistic - df) / c_statistic)
    i_squared_percent = max(0.0, (q_statistic - df) / q_statistic) * 100.0 if q_statistic > 0 else 0.0

    random_weights = 1.0 / (se_values**2 + tau_squared)
    combined = float(np.sum(random_weights * effect_values) / np.sum(random_weights))
    combined_se = float(np.sqrt(1.0 / np.sum(random_weights)))
    return {
        "effect": combined,
        "se": combined_se,
        "fixed_effect": fixed_effect,
        "q": q_statistic,
        "df": df,
        "i_squared_percent": float(i_squared_percent),
        "tau_squared": float(tau_squared),
        "n_studies": int(effect_values.size),
    }


# =============================================================================
# Per-feature comparative fits
# =============================================================================


def fit_study_effects(
    observations: pd.DataFrame,
    design: ComparativeDesign,
    exclude_studies: Sequence[Any] = (),
    exclude_covariates: Sequence[str] = (),
    feature: Any = None,
) -> list[dict[str, Any]]:
    """Estimate the contrast effect within every study (no p-values).

    Args:
        observations: Long-format DataFrame validated against the design.
        design: The predeclared :class:`ComparativeDesign`.
        exclude_studies: Study labels to leave out (sensitivity refits).
        exclude_covariates: Covariate columns to drop (sensitivity refits).
        feature: Feature label when ``design.feature_col`` is set.

    Returns:
        One dict per study (deterministic sorted order) with ``study``,
        ``effect``, ``se``, ``df``, and ``n_observations``.
    """
    _validate_design(design)
    _validate_observations(observations, design)
    frame = _feature_frame(observations, design, feature)
    excluded_studies = set(exclude_studies)
    excluded_covariates = frozenset(exclude_covariates)
    unknown_covariates = excluded_covariates - set(design.covariate_cols)
    if unknown_covariates:
        raise InferentialComparativeError(
            f"exclude_covariates names columns outside the design: {sorted(unknown_covariates)}"
        )
    fits: list[dict[str, Any]] = []
    for study in sorted(frame[design.study_col].unique(), key=str):
        if study in excluded_studies:
            continue
        sub = frame[frame[design.study_col] == study]
        fit = _fit_study(sub, design, excluded_covariates)
        fits.append({"study": study, **fit})
    if len(fits) < design.min_studies:
        raise InferentialComparativeError(
            f"only {len(fits)} usable studies remain, below the declared minimum of "
            f"{design.min_studies}; the analysis stops (plan section 8)"
        )
    return fits


def fit_comparative_effects(
    observations: pd.DataFrame,
    design: ComparativeDesign,
    exclude_studies: Sequence[Any] = (),
    exclude_covariates: Sequence[str] = (),
    feature: Any = None,
) -> dict[str, Any]:
    """Within-study fits plus the random-effects combination (no p-values).

    Returns:
        Dict with ``studies`` (per-study fits), ``combined`` (the
        :func:`random_effects_summary` dict), and ``n_observations``.
    """
    studies = fit_study_effects(observations, design, exclude_studies, exclude_covariates, feature)
    combined = random_effects_summary(
        [fit["effect"] for fit in studies],
        [fit["se"] for fit in studies],
    )
    return {
        "studies": studies,
        "combined": combined,
        "n_observations": int(sum(fit["n_observations"] for fit in studies)),
    }


# =============================================================================
# Bootstrap uncertainty (seeded)
# =============================================================================


def bootstrap_effect_ci(
    observations: pd.DataFrame,
    design: ComparativeDesign,
    random_seed: int,
    resampling_count: int,
    feature: Any = None,
) -> dict[str, Any]:
    """Percentile bootstrap CI of the combined effect (seeded RNG).

    Biological observations are resampled with replacement within each
    study; each replicate refits the within-study effects and the
    random-effects combination. Replicates whose resampled design is not
    estimable are skipped and counted. Emits NO p-value.

    Raises:
        InferentialComparativeError: If no replicate is estimable.
    """
    if not isinstance(random_seed, int) or isinstance(random_seed, bool) or random_seed < 0:
        raise ValueError(f"random_seed must be a non-negative integer, got {random_seed!r}")
    if not isinstance(resampling_count, int) or isinstance(resampling_count, bool) or resampling_count < 1:
        raise ValueError(f"resampling_count must be a positive integer, got {resampling_count!r}")
    _validate_design(design)
    _validate_observations(observations, design)

    frame = _feature_frame(observations, design, feature)
    studies = sorted(frame[design.study_col].unique(), key=str)
    # Precompute per-study arrays once; the resampling loop is pure numpy.
    arrays_by_study = {
        study: _study_arrays(frame[frame[design.study_col] == study], design, frozenset())
        for study in studies
    }
    rng = np.random.default_rng(random_seed)
    values: list[float] = []
    for _ in range(resampling_count):
        replicate_fits: list[dict[str, Any]] = []
        for study in studies:
            matrix, response, names = arrays_by_study[study]
            pick = rng.integers(0, matrix.shape[0], matrix.shape[0])
            try:
                replicate_fits.append(_ols_effect(matrix[pick], response[pick], names))
            except InferentialComparativeError:
                continue  # resampling artifact; the study is skipped, not silenced
        if len(replicate_fits) < 2:
            continue
        summary = random_effects_summary(
            [fit["effect"] for fit in replicate_fits],
            [fit["se"] for fit in replicate_fits],
        )
        values.append(summary["effect"])
    if not values:
        raise InferentialComparativeError(
            "all bootstrap replicates failed to fit; the design cannot support "
            "bootstrap uncertainty (plan section 8 stopping rule)"
        )
    distribution = np.asarray(values, dtype=float)
    low, high = np.percentile(distribution, [2.5, 97.5])
    return {
        "ci_low": float(low),
        "ci_high": float(high),
        "n_success": int(distribution.size),
        "resampling_count": int(resampling_count),
        "random_seed": int(random_seed),
    }


# =============================================================================
# Gated inferential analysis
# =============================================================================


def run_inferential_comparative_analysis(
    observations: pd.DataFrame,
    design: ComparativeDesign,
    contract: AnalysisProvenance,
    evidence_manifest_frozen: bool = False,
) -> dict[str, Any]:
    """GATED: full inferential comparative analysis under a frozen contract.

    For each feature: within-study OLS contrast fits, the DerSimonian-Laird
    random-effects combination (effect, SE, Cochran's Q, I-squared,
    tau-squared), a seeded percentile bootstrap CI, and the two-sided
    z-test of the combined effect against zero. Raw p-values are adjusted
    only through :func:`declared_inferential_bh_fdr`, which re-validates
    the contract and the tested-feature family.

    Args:
        observations: Long-format DataFrame (one row per biological
            observation, optionally per feature).
        design: The predeclared :class:`ComparativeDesign`.
        contract: A validated :class:`AnalysisProvenance` declared
            ``analysis_role="inferential"`` with a BH-FDR procedure whose
            ``tested_feature_count`` equals the number of features present.
        evidence_manifest_frozen: Must be True; the caller affirms the
            evidence manifest has frozen and inferential output is allowed.

    Returns:
        Dict with ``role="inferential"``, the gate label, contract echoes,
        and ``features``: a DataFrame indexed by feature with columns
        ``effect``, ``se``, ``ci_low``, ``ci_high``, ``bootstrap_n_success``,
        ``z``, ``p_value``, ``p_adj_bh``, ``q``, ``df``, ``i_squared_percent``,
        ``tau_squared``, ``n_studies``, ``n_observations``. The frame carries
        ``attrs["role"] == "inferential"``.

    Raises:
        RuntimeError: If the evidence manifest is not affirmed frozen.
        ProvenanceError: If the contract fails validation.
        StatisticsContractError: On a non-inferential or non-BH contract, a
            tested-feature count mismatch, or a multiplicity failure.
        InferentialComparativeError: On any model that cannot be estimated
            as declared.
    """
    require_inferential_contract(contract, evidence_manifest_frozen)
    _validate_design(design)
    _validate_observations(observations, design)
    features = _feature_labels(observations, design)
    if contract.tested_feature_count != len(features):
        raise StatisticsContractError(
            f"contract declares tested_feature_count={contract.tested_feature_count} but "
            f"the design presents {len(features)} features; refusing to infer a different family"
        )

    rows: dict[Any, dict[str, Any]] = {}
    for feature in features:
        fit = fit_comparative_effects(observations, design, feature=feature)
        bootstrap = bootstrap_effect_ci(
            observations,
            design,
            random_seed=contract.random_seed,
            resampling_count=contract.resampling_count,
            feature=feature,
        )
        combined = fit["combined"]
        z_statistic = combined["effect"] / combined["se"]
        rows[feature if feature is not None else design.response_col] = {
            "effect": combined["effect"],
            "se": combined["se"],
            "ci_low": bootstrap["ci_low"],
            "ci_high": bootstrap["ci_high"],
            "bootstrap_n_success": bootstrap["n_success"],
            "z": float(z_statistic),
            "p_value": float(2.0 * stats.norm.sf(abs(z_statistic))),
            "q": combined["q"],
            "df": combined["df"],
            "i_squared_percent": combined["i_squared_percent"],
            "tau_squared": combined["tau_squared"],
            "n_studies": combined["n_studies"],
            "n_observations": fit["n_observations"],
        }

    frame = pd.DataFrame.from_dict(rows, orient="index")
    frame.index.name = design.feature_col if design.feature_col is not None else "feature"
    frame.attrs["role"] = INFERENTIAL_ROLE

    multiplicity = declared_inferential_bh_fdr(
        frame["p_value"].tolist(),
        contract,
        evidence_manifest_frozen=True,
    )
    frame["p_adj_bh"] = multiplicity["adjusted_p_values"]

    return {
        "role": INFERENTIAL_ROLE,
        "gate": "post-freeze",
        "analysis_id": contract.analysis_id,
        "estimand": contract.estimand,
        "replicate_unit": contract.replicate_unit,
        "random_seed": contract.random_seed,
        "resampling_count": contract.resampling_count,
        "multiple_testing_family": contract.multiple_testing_family,
        "multiple_testing_method": multiplicity["multiple_testing_method"],
        "n_features": len(features),
        "features": frame,
    }


# =============================================================================
# Sensitivity runner (plan section 7)
# =============================================================================


def directional_agreement(expected_direction: str, observed_direction: str) -> bool:
    """Report whether an observed effect movement matches the predeclared one.

    ``"either"`` accepts any movement; any other expected direction agrees
    only with the identical observed direction (``"increase"``,
    ``"decrease"``, or ``"none"``).
    """
    if expected_direction == "either":
        return True
    return expected_direction == observed_direction


def _observed_direction(baseline_effect: float, varied_effect: float) -> str:
    if varied_effect > baseline_effect:
        return "increase"
    if varied_effect < baseline_effect:
        return "decrease"
    return "none"


def _sensitivity_comparisons(
    observations: pd.DataFrame,
    design: ComparativeDesign,
    entry: SensitivityAnalysis,
    features: list[Any],
    baseline_effects: dict[Any, float],
    exclude_studies_for_value: Any,
    exclude_covariates_for_value: Any,
) -> list[dict[str, Any]]:
    comparisons: list[dict[str, Any]] = []
    for value in entry.varied_values:
        for feature in features:
            varied = fit_comparative_effects(
                observations,
                design,
                exclude_studies=exclude_studies_for_value(value),
                exclude_covariates=exclude_covariates_for_value(value),
                feature=feature,
            )
            baseline_effect = baseline_effects[feature]
            varied_effect = varied["combined"]["effect"]
            observed = _observed_direction(baseline_effect, varied_effect)
            comparisons.append(
                {
                    "varied_value": value,
                    "feature": feature if feature is not None else design.response_col,
                    "baseline_effect": baseline_effect,
                    "varied_effect": varied_effect,
                    "effect_change": varied_effect - baseline_effect,
                    "observed_direction": observed,
                    "expected_direction": entry.expected_direction,
                    "directional_agreement": directional_agreement(entry.expected_direction, observed),
                    "n_studies_varied": varied["combined"]["n_studies"],
                }
            )
    return comparisons


def _sensitivity_leave_one_study_out(
    observations: pd.DataFrame,
    design: ComparativeDesign,
    entry: SensitivityAnalysis,
    features: list[Any],
    baseline_effects: dict[Any, float],
) -> list[dict[str, Any]]:
    studies = sorted(observations[design.study_col].unique(), key=str)
    if len(studies) - 1 < 2:
        raise InferentialComparativeError(
            f"leave-one-study-out needs at least 3 studies so a refit keeps 2, got {len(studies)}"
        )
    known = set(studies)

    def exclude_studies_for_value(value: str) -> tuple[Any, ...]:
        if value not in known:
            raise InferentialComparativeError(
                f"sensitivity entry {entry.name!r} varies study {value!r}, which is not "
                f"present in the observations (known: {sorted(known, key=str)})"
            )
        return (value,)

    def exclude_covariates_for_value(value: str) -> frozenset[str]:
        return frozenset()

    return _sensitivity_comparisons(
        observations,
        design,
        entry,
        features,
        baseline_effects,
        exclude_studies_for_value,
        exclude_covariates_for_value,
    )


def _sensitivity_exclude_covariate(
    observations: pd.DataFrame,
    design: ComparativeDesign,
    entry: SensitivityAnalysis,
    features: list[Any],
    baseline_effects: dict[Any, float],
) -> list[dict[str, Any]]:
    known = set(design.covariate_cols)

    def exclude_studies_for_value(value: str) -> tuple[Any, ...]:
        return ()

    def exclude_covariates_for_value(value: str) -> frozenset[str]:
        if value not in known:
            raise InferentialComparativeError(
                f"sensitivity entry {entry.name!r} varies covariate {value!r}, which the "
                f"design does not declare (declared: {sorted(known)})"
            )
        return frozenset({value})

    return _sensitivity_comparisons(
        observations,
        design,
        entry,
        features,
        baseline_effects,
        exclude_studies_for_value,
        exclude_covariates_for_value,
    )


_SENSITIVITY_HANDLERS = {
    "leave_one_study_out": _sensitivity_leave_one_study_out,
    "exclude_covariate": _sensitivity_exclude_covariate,
}


def run_registered_sensitivity_analyses(
    observations: pd.DataFrame,
    design: ComparativeDesign,
    contract: AnalysisProvenance,
    evidence_manifest_frozen: bool = False,
) -> list[dict[str, Any]]:
    """GATED: execute the contract's registered sensitivity analyses.

    Each registered :class:`SensitivityAnalysis` is validated
    (:func:`validate_sensitivity_analysis`), its ``varied_parameter`` must
    name an executable variation (``leave_one_study_out`` or
    ``exclude_covariate``), and each varied refit is compared with the
    primary fit by directional agreement. Sensitivity results carry no
    p-values; they report effects and observed vs expected direction.

    Returns:
        One dict per registered entry with ``role="inferential"``, the
        entry fields, ``n_comparisons``, ``n_agreement``,
        ``fraction_agreement``, and per-comparison records.

    Raises:
        RuntimeError: If the evidence manifest is not affirmed frozen.
        ProvenanceError: If the contract or an entry fails validation.
        StatisticsContractError: On a non-inferential contract or a
            tested-feature count mismatch.
        InferentialComparativeError: On an unsupported ``varied_parameter``
            or a varied value absent from the data/design.
    """
    require_inferential_contract(contract, evidence_manifest_frozen)
    _validate_design(design)
    _validate_observations(observations, design)
    features = _feature_labels(observations, design)
    if contract.tested_feature_count != len(features):
        raise StatisticsContractError(
            f"contract declares tested_feature_count={contract.tested_feature_count} but "
            f"the design presents {len(features)} features; refusing to infer a different family"
        )

    baseline_effects = {
        feature: fit_comparative_effects(observations, design, feature=feature)["combined"]["effect"]
        for feature in features
    }

    results: list[dict[str, Any]] = []
    for entry in contract.sensitivity_analyses:
        validate_sensitivity_analysis(entry)
        handler = _SENSITIVITY_HANDLERS.get(entry.varied_parameter)
        if handler is None:
            raise InferentialComparativeError(
                f"sensitivity entry {entry.name!r} declares varied_parameter="
                f"{entry.varied_parameter!r}, which this runner cannot execute; "
                f"supported parameters: {sorted(_SENSITIVITY_HANDLERS)}"
            )
        comparisons = handler(observations, design, entry, features, baseline_effects)
        n_agreement = sum(1 for comparison in comparisons if comparison["directional_agreement"])
        results.append(
            {
                "role": INFERENTIAL_ROLE,
                "name": entry.name,
                "varied_parameter": entry.varied_parameter,
                "baseline_value": entry.baseline_value,
                "expected_direction": entry.expected_direction,
                "n_comparisons": len(comparisons),
                "n_agreement": n_agreement,
                "fraction_agreement": n_agreement / len(comparisons) if comparisons else 0.0,
                "comparisons": comparisons,
            }
        )
    return results
