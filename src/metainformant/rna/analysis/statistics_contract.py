"""Predeclared statistical contract for cross-species RNA analyses.

Implements the analysis-provenance and descriptive/inferential boundary
required by
``projects/hymenoptera_amalgkit/docs/manuscript/statistical_analysis_plan.md``
(sections 1, 4, 6, 7, 8, and 9):

- :class:`AnalysisProvenance` is the structured record declared before an
  analysis runs: analysis identifier, estimand, biological replicate unit,
  random seed, resampling count, null model, multiple-testing procedure
  (family, method, tested-feature count), analysis role, and software
  versions.
- :class:`SensitivityAnalysis` registers one predeclared robustness check
  (plan section 7): the varied parameter, its baseline and varied values,
  and the direction the primary estimand is expected to move. Registered
  entries are validated fail-closed and rendered as additive
  ``analysis_provenance_sensitivity_*`` lines.
- :class:`ReplicateUnitDeclaration` and :class:`EstimandDeclaration` are the
  structured declaration API for the biological replicate unit and the
  estimand (plan sections 1 and 3): the caller declares the smallest
  defensible unit, its nesting, the independent-replicate counts per
  sampling stratum, and the permitted interpretation boundary; the contract
  validates non-degeneracy fail-closed (no placeholder labels, technical
  replicates never counted as independent, no stratum below the declared
  minimum, inferential roles require at least 2 replicates per unit).
  Declared entries render as additive
  ``analysis_provenance_replicate_unit_*`` and
  ``analysis_provenance_estimand_*`` lines.
- Records may declare the non-analysis roles ``"stopped"``/``"unavailable"``
  (plan section 8) to record a halted or impossible analysis explicitly;
  such a record must not declare multiplicity, cohort-denominator,
  artifact-path, or sensitivity fields that would imply results exist.
- Descriptive outputs (the fingerprint divergence matrix and the
  feature-resampling sensitivity table in
  :mod:`metainformant.rna.analysis.cross_species`) carry
  ``attrs["role"] == "descriptive"``. They are permutation/sensitivity
  scores, never p-values, and never confidence intervals.
- The only inferential path is :func:`declared_inferential_bh_fdr`, which is
  GATED (mirroring ``wilcoxon_duplication_specificity`` in
  ``tissue_specificity.py``): it refuses to run unless the caller passes an
  already-validated contract declared ``inferential`` AND sets
  ``evidence_manifest_frozen=True``. It applies the declared
  Benjamini-Hochberg FDR procedure before returning adjusted values.
- :func:`validate_orthology_profile_invariants` and
  :func:`validate_species_tree_invariants` fail closed with explicit error
  types on orthology-bridge and species-tree violations (plan sections 4
  and 8).
- predeclared MJ-01 comparative designs (:class:`PredeclaredDesign`): the
  study, contrast, and covariate columns (tissue, caste, sex, stage) and
  their exact strata levels are declared up front;
  :func:`validate_predeclared_design` enforces structural non-degeneracy
  and :func:`enforce_predeclared_design` cross-checks the declaration
  against the observations before any fit — unknown covariate names and
  undeclared strata fail closed, and empty/singleton strata raise
  :class:`EmptyStratumError`/:class:`SingletonStratumError` instead of
  being silently dropped.
- role-conditional paired effect sizes (:func:`declared_effect_size` with
  the pure helper :func:`paired_log2fc_effect`): paired log2 fold-change
  with an analytic standard error and, for inferential contracts, a
  seeded percentile bootstrap CI gated post-freeze. Descriptive contracts
  render every inferential field as ``not-applicable``; inferential
  records carry the effect-size method, the paired replicate count, and
  the multiplicity family/method provenance.
- multi-study heterogeneity records (:func:`heterogeneity_record`):
  Cochran's Q, degrees of freedom, I-squared, and tau-squared for
  inferential contracts, fail-closed below 2 studies;
  descriptive-role records are exempt and render ``not-applicable``.
- a leave-one-study-out sensitivity hook
  (:func:`leave_one_study_out_deltas`) recomputing the random-effects
  combined effect per exclusion and returning per-exclusion deltas (pure
  computation, no I/O).

Every validator fails closed: a record with missing or placeholder fields,
or an invariant violation, raises before any result can be produced.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import pandas as pd

__all__ = [
    "AnalysisProvenance",
    "DESCRIPTIVE_ROLE",
    "DesignDeclarationError",
    "EmptyStratumError",
    "EstimandDeclaration",
    "HeterogeneityError",
    "INFERENTIAL_ROLE",
    "PredeclaredDesign",
    "ProvenanceError",
    "ReplicateUnitDeclaration",
    "SensitivityAnalysis",
    "SingletonStratumError",
    "StatisticsContractError",
    "TreeInvariantError",
    "UndeclaredStratumError",
    "UnknownDesignCovariateError",
    "benjamini_hochberg_fdr",
    "declared_effect_size",
    "declared_inferential_bh_fdr",
    "enforce_predeclared_design",
    "heterogeneity_record",
    "leave_one_study_out_deltas",
    "paired_log2fc_effect",
    "render_analysis_provenance_block",
    "render_effect_size_record",
    "result_role",
    "validate_analysis_provenance",
    "validate_estimand_declaration",
    "validate_orthology_profile_invariants",
    "validate_predeclared_design",
    "validate_replicate_unit_declaration",
    "validate_sensitivity_analysis",
    "validate_species_tree_invariants",
]

DESCRIPTIVE_ROLE = "descriptive"
INFERENTIAL_ROLE = "inferential"
# Predeclared non-analysis states (statistical_analysis_plan.md section 8):
# an analysis that was declared but halted before producing results is
# recorded as stopped/unavailable instead of silently disappearing from an
# evidence bundle.
STOPPED_ROLE = "stopped"
UNAVAILABLE_ROLE = "unavailable"

_ALLOWED_ROLES = frozenset({DESCRIPTIVE_ROLE, INFERENTIAL_ROLE, STOPPED_ROLE, UNAVAILABLE_ROLE})
NON_ANALYSIS_ROLES = frozenset({STOPPED_ROLE, UNAVAILABLE_ROLE})
# Predeclared multiplicity procedures (statistical_analysis_plan.md section 6).
_ALLOWED_MT_METHODS = frozenset({"bh-fdr", "benjamini-hochberg", "bonferroni"})
# Strings that make a declared field a placeholder rather than a declaration.
_PLACEHOLDER_STRINGS = frozenset(
    {"", "na", "n/a", "none", "null", "todo", "tbd", "placeholder", "unknown", "pending", "?"}
)

# Predeclared sensitivity-analysis expectation directions (plan section 7):
# what the primary estimand is expected to do under the variation if the
# result is robust. Anything outside this set is not a declaration.
_ALLOWED_SENSITIVITY_DIRECTIONS = frozenset({"increase", "decrease", "either", "none"})


class StatisticsContractError(ValueError):
    """Base class for fail-closed statistical-contract violations."""


class ProvenanceError(StatisticsContractError):
    """An analysis-provenance record is missing or contains placeholder fields."""


class OrthologyInvariantError(StatisticsContractError):
    """An orthology x species presence table violates a declared invariant."""


class TreeInvariantError(StatisticsContractError):
    """A species tree is unrooted, malformed, or has conflicting labels."""


class DesignDeclarationError(StatisticsContractError):
    """A predeclared comparative-design declaration is violated."""


class UnknownDesignCovariateError(DesignDeclarationError):
    """A design column is named outside the predeclared covariate declaration."""


class UndeclaredStratumError(DesignDeclarationError):
    """Observations carry a stratum level the predeclared design did not declare."""


class EmptyStratumError(DesignDeclarationError):
    """A declared stratum level contributes no observations."""


class SingletonStratumError(DesignDeclarationError):
    """A declared stratum level is too thin to support the declared contrast."""


class HeterogeneityError(StatisticsContractError):
    """A multi-study heterogeneity invariant is violated (fail-closed)."""


@dataclass(frozen=True)
class SensitivityAnalysis:
    """Predeclared sensitivity analysis (plan section 7).

    Declares, before the analysis runs, one parameter that is varied around
    a baseline, the exact values it takes, and the direction the primary
    estimand is expected to move under the variation if the result is
    robust. The record is frozen like :class:`AnalysisProvenance`: it is
    part of the predeclared contract, not a post-hoc narrative.
    """

    name: str
    varied_parameter: str
    baseline_value: str
    varied_values: tuple[str, ...]
    expected_direction: str
    notes: str = ""


@dataclass(frozen=True)
class ReplicateUnitDeclaration:
    """Declared biological-replicate unit (plan sections 1 and 3).

    The biological replicate is not automatically a downloaded library:
    public libraries may be pooled individuals, technical replicates,
    repeated measurements, or different biological states. The caller
    declares the smallest defensible biological unit, how it nests (for
    example ``"library nested in study nested in species"``), and the
    number of *independent biological replicates* each sampling stratum
    actually contributes. The declaration is frozen: it is part of the
    predeclared contract, validated for non-degeneracy by
    :func:`validate_replicate_unit_declaration`.
    """

    name: str
    nesting: str
    counts: Mapping[str, int]
    min_independent_replicates: int = 1
    technical_replicates_counted: bool = False


@dataclass(frozen=True)
class EstimandDeclaration:
    """Declared estimand (plan section 1).

    Declares, before the analysis runs, the quantity being estimated, the
    biological contrast it applies to, and the permitted interpretation
    boundary (which claims the analysis may and may not support). It may
    embed the :class:`ReplicateUnitDeclaration` the estimand is defined
    over; an embedded unit that conflicts with the provenance record's own
    unit declaration fails closed.
    """

    name: str
    contrast: str
    permitted_interpretation: str
    replicate_unit: ReplicateUnitDeclaration | None = None


@dataclass(frozen=True)
class PredeclaredDesign:
    """Predeclared comparative design for the MJ-01 inferential lane.

    Maps every design column — the study column, the contrast column, and
    the covariate columns (tissue, caste, sex, stage) — to the exact strata
    levels the analysis is predeclared to encounter. The declaration is
    frozen: it is part of the predeclared contract, not a post-hoc
    description. :func:`validate_predeclared_design` enforces structural
    non-degeneracy, and :func:`enforce_predeclared_design` cross-checks the
    declaration against the observations before any fit runs; unknown
    covariate names and undeclared strata fail closed, and empty/singleton
    strata raise :class:`EmptyStratumError`/:class:`SingletonStratumError`
    instead of being silently dropped.
    """

    covariate_strata: Mapping[str, tuple[str, ...]]


@dataclass(frozen=True)
class AnalysisProvenance:
    """Predeclared analysis-provenance record (plan sections 1, 4, 6, 9).

    The record is frozen: after declaration it is the immutable contract the
    analysis runs under. Descriptive analyses declare
    ``analysis_role="descriptive"``; only contracts declared
    ``analysis_role="inferential"`` with a BH-FDR procedure may be used with
    :func:`declared_inferential_bh_fdr`.

    The optional reporting-contract fields (plan section 9) bind the record
    to the data-root snapshot id, cohort inclusion/exclusion denominators,
    exact artifact paths, the metadata-harmonization review state (plan
    section 3), and the species-tree source and branch-length scale (plan
    section 4); each is validated fail-closed when declared and rendered
    only when declared. ``sensitivity_analyses`` registers frozen
    :class:`SensitivityAnalysis` entries alongside the record. The
    non-analysis roles ``"stopped"``/``"unavailable"`` (plan section 8)
    record a halted analysis and must not declare any field that implies
    results exist.
    """

    analysis_id: str
    estimand: str
    replicate_unit: str
    random_seed: int
    resampling_count: int
    null_model: str
    # Multiplicity provenance defaults to None — the descriptive-role
    # representation (rendered ``not-applicable``). The fields stay optional
    # at construction so legacy and descriptive-role callers keep working
    # unchanged; :func:`validate_analysis_provenance` enforces the
    # ROLE-conditional requirement instead (an inferential record must
    # declare a real family, method, and tested-feature count, and a
    # descriptive record must declare none).
    multiple_testing_family: str | None = None
    multiple_testing_method: str | None = None
    tested_feature_count: int | None = None
    # A record without software versions is never a declared analysis;
    # validate_analysis_provenance refuses an empty mapping for every role.
    software_versions: Mapping[str, str] = field(default_factory=dict)
    analysis_role: str = DESCRIPTIVE_ROLE
    # Reporting-contract fields (plan section 9): bind the analysis to the
    # data-root snapshot it ran against and to the exact cohort denominators
    # and artifact paths. All optional; declared values are validated
    # fail-closed (no placeholders) so a record cannot imply evidence that
    # was not bound.
    data_root_snapshot_id: str | None = None
    cohort_included_count: int | None = None
    cohort_excluded_count: int | None = None
    artifact_paths: Mapping[str, str] | None = None
    # Metadata-harmonization review flag (plan section 3): records the
    # review state of the harmonization table the analysis consumed.
    metadata_harmonization_review: str | None = None
    # Species-tree binding (plan section 4 / methods 'versioned species
    # tree'): source identity plus branch-length scale so a tree-dependent
    # analysis names the exact tree it used.
    species_tree_source: str | None = None
    species_tree_branch_length_scale: str | None = None
    # Sensitivity-analysis registry (plan section 7): predeclared
    # robustness checks that vary one parameter around a baseline. Empty by
    # default; declared entries are validated fail-closed.
    sensitivity_analyses: tuple[SensitivityAnalysis, ...] = ()
    # Structured biological-replicate unit and estimand declarations (plan
    # sections 1 and 3): refinements of the ``replicate_unit`` and
    # ``estimand`` strings above. Optional; declared values are validated
    # fail-closed for non-degeneracy (see
    # :func:`validate_replicate_unit_declaration` and
    # :func:`validate_estimand_declaration`) and rendered additively.
    replicate_unit_declaration: ReplicateUnitDeclaration | None = None
    estimand_declaration: EstimandDeclaration | None = None
    # Predeclared comparative design (MJ-01): the exact design columns
    # (study, contrast, covariates) and their strata levels the inferential
    # comparative lane is predeclared to encounter. Optional; declared
    # values are validated fail-closed by
    # :func:`validate_predeclared_design` and rendered additively; a
    # non-analysis role must not declare one.
    design_declaration: PredeclaredDesign | None = None


def _require_declared(value: object, field: str) -> str:
    """Reject missing or placeholder strings for a declared field."""
    if not isinstance(value, str) or value.strip().lower() in _PLACEHOLDER_STRINGS:
        raise ProvenanceError(
            f"analysis provenance field '{field}' is missing or a placeholder "
            f"(got {value!r}); a declared value is required"
        )
    return value


def validate_sensitivity_analysis(analysis: SensitivityAnalysis) -> None:
    """Fail closed on missing, placeholder, or inconsistent sensitivity fields.

    Raises:
        TypeError: If ``analysis`` is not a :class:`SensitivityAnalysis`.
        ProvenanceError: If the name, varied parameter, or baseline value is
            missing or a placeholder; if ``varied_values`` is not a
            non-empty tuple of declared strings; or if
            ``expected_direction`` is not one of the allowed directions.
    """
    if not isinstance(analysis, SensitivityAnalysis):
        raise TypeError(
            "validate_sensitivity_analysis requires a SensitivityAnalysis record; "
            "ad-hoc dictionaries cannot predeclare a sensitivity check"
        )
    _require_declared(analysis.name, "sensitivity name")
    _require_declared(analysis.varied_parameter, "sensitivity varied_parameter")
    _require_declared(analysis.baseline_value, "sensitivity baseline_value")
    if not isinstance(analysis.varied_values, tuple) or not analysis.varied_values:
        raise ProvenanceError(
            "sensitivity varied_values must be a non-empty tuple of varied values, " f"got {analysis.varied_values!r}"
        )
    for value in analysis.varied_values:
        _require_declared(value, "sensitivity varied_values entry")
    if analysis.expected_direction not in _ALLOWED_SENSITIVITY_DIRECTIONS:
        raise ProvenanceError(
            f"sensitivity expected_direction must be one of "
            f"{sorted(_ALLOWED_SENSITIVITY_DIRECTIONS)}, got {analysis.expected_direction!r}"
        )


def validate_replicate_unit_declaration(
    unit: ReplicateUnitDeclaration,
    *,
    require_inferential_grade: bool = False,
) -> None:
    """Fail closed on missing, placeholder, or degenerate replicate-unit declarations.

    Non-degeneracy is structural (plan sections 1 and 3):

    - ``name`` and ``nesting`` must be declared strings (never placeholders);
    - ``technical_replicates_counted`` must be ``False``: pooled
      individuals, technical replicates, and repeated measurements are not
      independent biological replicates;
    - ``counts`` must be a non-empty mapping of declared stratum labels to
      positive integers, and every stratum must supply at least
      ``min_independent_replicates`` independent biological replicates — a
      stratum with fewer has no defensible replicate structure and cannot
      support the same contrast as replicated strata;
    - ``min_independent_replicates`` must be a positive integer; with
      ``require_inferential_grade=True`` (inferential roles) it must be at
      least 2, because a single replicate cannot define an effect-size
      estimand with uncertainty.

    Raises:
        TypeError: If ``unit`` is not a :class:`ReplicateUnitDeclaration`.
        ProvenanceError: On any violated invariant.
    """
    if not isinstance(unit, ReplicateUnitDeclaration):
        raise TypeError(
            "validate_replicate_unit_declaration requires a ReplicateUnitDeclaration record; "
            "ad-hoc dictionaries cannot predeclare the biological replicate unit"
        )
    _require_declared(unit.name, "replicate unit name")
    _require_declared(unit.nesting, "replicate unit nesting")
    if unit.technical_replicates_counted:
        raise ProvenanceError(
            "replicate unit declares technical_replicates_counted=True; pooled individuals, "
            "technical replicates, and repeated measurements are not independent biological "
            "replicates (plan section 1)"
        )
    if (
        not isinstance(unit.min_independent_replicates, int)
        or isinstance(unit.min_independent_replicates, bool)
        or unit.min_independent_replicates < 1
    ):
        raise ProvenanceError(
            f"min_independent_replicates must be a positive integer, got " f"{unit.min_independent_replicates!r}"
        )
    if require_inferential_grade and unit.min_independent_replicates < 2:
        raise ProvenanceError(
            "an inferential effect-size estimand requires at least 2 independent "
            f"biological replicates per declared unit; got min_independent_replicates="
            f"{unit.min_independent_replicates!r}"
        )
    if not isinstance(unit.counts, Mapping) or not unit.counts:
        raise ProvenanceError(
            "replicate unit counts must be a non-empty mapping of stratum label to "
            f"independent biological replicate count, got {unit.counts!r}"
        )
    degenerate: list[str] = []
    for stratum, count in unit.counts.items():
        _require_declared(stratum, "replicate unit stratum label")
        if not isinstance(count, int) or isinstance(count, bool) or count < 1:
            raise ProvenanceError(f"replicate unit counts[{stratum!r}] must be a positive integer, got {count!r}")
        if count < unit.min_independent_replicates:
            degenerate.append(f"{stratum} ({count} < {unit.min_independent_replicates})")
    if degenerate:
        raise ProvenanceError(
            "replicate unit is degenerate for the declared estimand; strata below the "
            f"declared minimum of {unit.min_independent_replicates} independent biological "
            f"replicates: {', '.join(degenerate)}"
        )


def validate_estimand_declaration(estimand: EstimandDeclaration) -> None:
    """Fail closed on missing, placeholder, or inconsistent estimand declarations.

    The estimand is declared before the analysis runs (plan section 1): the
    quantity, the biological contrast it applies to, and the permitted
    interpretation boundary. If the declaration embeds a
    :class:`ReplicateUnitDeclaration`, it is validated in full.

    Raises:
        TypeError: If ``estimand`` is not an :class:`EstimandDeclaration`.
        ProvenanceError: On any missing or placeholder field.
    """
    if not isinstance(estimand, EstimandDeclaration):
        raise TypeError(
            "validate_estimand_declaration requires an EstimandDeclaration record; "
            "ad-hoc dictionaries cannot predeclare the estimand"
        )
    _require_declared(estimand.name, "estimand name")
    _require_declared(estimand.contrast, "estimand contrast")
    _require_declared(estimand.permitted_interpretation, "estimand permitted_interpretation")
    if estimand.replicate_unit is not None:
        validate_replicate_unit_declaration(estimand.replicate_unit)


def validate_predeclared_design(design: PredeclaredDesign) -> None:
    """Fail closed on missing, placeholder, or degenerate design declarations.

    Structural non-degeneracy for the predeclared MJ-01 comparative design:
    every covariate name must be a declared string (never a placeholder),
    and every strata value must be a non-empty tuple of unique declared
    level strings.

    Raises:
        TypeError: If ``design`` is not a :class:`PredeclaredDesign`.
        ProvenanceError: On any placeholder covariate name, a non-tuple or
            empty strata value, or a placeholder/repeated stratum level.
    """
    if not isinstance(design, PredeclaredDesign):
        raise TypeError(
            "validate_predeclared_design requires a PredeclaredDesign record; "
            "ad-hoc dictionaries cannot predeclare the comparative design"
        )
    if not isinstance(design.covariate_strata, Mapping) or not design.covariate_strata:
        raise ProvenanceError(
            "covariate_strata must be a non-empty mapping of covariate name to "
            f"declared strata levels, got {design.covariate_strata!r}"
        )
    for covariate, strata in design.covariate_strata.items():
        _require_declared(covariate, "predeclared design covariate name")
        if not isinstance(strata, tuple) or not strata:
            raise ProvenanceError(
                f"predeclared design strata for covariate {covariate!r} must be a "
                f"non-empty tuple of declared levels, got {strata!r}"
            )
        seen: set[str] = set()
        for level in strata:
            _require_declared(level, f"predeclared design stratum level for {covariate!r}")
            if level in seen:
                raise ProvenanceError(
                    f"predeclared design strata for covariate {covariate!r} repeat " f"level {level!r}"
                )
            seen.add(level)


def enforce_predeclared_design(
    design: PredeclaredDesign,
    observations: pd.DataFrame,
    *,
    study_col: str,
    contrast_col: str | None = None,
    covariate_cols: Sequence[str] = (),
    min_stratum_observations: int = 2,
) -> dict[str, dict[str, int]]:
    """Cross-check a predeclared design against observations, fail closed.

    Every design column — the study column, the contrast column when given,
    and every covariate column — must exactly match the predeclared
    declaration, must be a column of ``observations``, and every observed
    stratum level must be declared with at least
    ``min_stratum_observations`` observations. Nothing is silently dropped:
    a level observed but not declared raises
    :class:`UndeclaredStratumError`; a declared level with no observations
    raises :class:`EmptyStratumError`; a declared level observed fewer than
    ``min_stratum_observations`` times raises
    :class:`SingletonStratumError`.

    Args:
        design: The predeclared :class:`PredeclaredDesign` declaration.
        observations: Long-format observations table.
        study_col: The observations column grouping observations into studies.
        contrast_col: The observations column holding the contrast levels.
        covariate_cols: The observations columns retained as covariates.
        min_stratum_observations: Minimum observation count per stratum
            level; the default of 2 refuses singleton strata.

    Returns:
        Deterministic ``{column: {level: observation_count}}`` mapping,
        columns and levels sorted.

    Raises:
        TypeError: If ``design`` is not a :class:`PredeclaredDesign` or
            ``observations`` is not a DataFrame.
        ProvenanceError: If ``min_stratum_observations`` is not a positive
            integer, or the declaration itself is degenerate.
        UnknownDesignCovariateError: If the requested design columns do not
            match the declaration, or a declared column is absent from the
            observations.
        UndeclaredStratumError: On an observed level outside the declaration.
        EmptyStratumError: On a declared level with zero observations.
        SingletonStratumError: On a declared level below the minimum count.
    """
    validate_predeclared_design(design)
    if not isinstance(observations, pd.DataFrame):
        raise TypeError(
            "enforce_predeclared_design requires an observations DataFrame; " f"got {type(observations).__name__}"
        )
    if (
        not isinstance(min_stratum_observations, int)
        or isinstance(min_stratum_observations, bool)
        or min_stratum_observations < 1
    ):
        raise ProvenanceError(
            f"min_stratum_observations must be a positive integer, got " f"{min_stratum_observations!r}"
        )
    requested = sorted({str(name) for name in covariate_cols} | {str(study_col)})
    if contrast_col is not None:
        requested.append(str(contrast_col))
        requested = sorted(requested)
    declared = sorted(str(name) for name in design.covariate_strata)
    if requested != declared:
        undeclared = [name for name in requested if name not in declared]
        missing = [name for name in declared if name not in requested]
        raise UnknownDesignCovariateError(
            "the requested design columns do not match the predeclared design: "
            f"columns outside the declaration: {undeclared}; declared columns "
            f"absent from the request: {missing}"
        )
    counts: dict[str, dict[str, int]] = {}
    for name in declared:
        if name not in observations.columns:
            raise UnknownDesignCovariateError(f"predeclared design column {name!r} is not a column of the observations")
        observed = observations[name]
        strata = sorted(str(level) for level in design.covariate_strata[name])
        undeclared_levels = sorted((str(level) for level in observed.unique() if str(level) not in strata))
        if undeclared_levels:
            raise UndeclaredStratumError(
                f"observations carry strata the predeclared design did not declare "
                f"for column {name!r}: {undeclared_levels}; refusing to drop them "
                "silently"
            )
        stratum_counts: dict[str, int] = {}
        for level in strata:
            count = int((observed == level).sum())
            if count == 0:
                raise EmptyStratumError(
                    f"predeclared design stratum {level!r} of column {name!r} has no "
                    "observations; the declaration misstates the cohort"
                )
            if count < min_stratum_observations:
                raise SingletonStratumError(
                    f"predeclared design stratum {level!r} of column {name!r} has "
                    f"{count} observation(s), below the required minimum of "
                    f"{min_stratum_observations}; a stratum this thin cannot support "
                    "the declared contrast and is never silently dropped"
                )
            stratum_counts[level] = count
        counts[name] = stratum_counts
    return counts


def validate_analysis_provenance(record: AnalysisProvenance) -> None:
    """Fail closed on missing, placeholder, or inconsistent provenance fields.

    Raises:
        ProvenanceError: On any missing/placeholder field, a negative seed,
            a role mismatch, or an empty software-version mapping. The
            multiplicity fields are role-conditional and fail closed in both
            directions: ``multiple_testing_family``,
            ``multiple_testing_method``, and ``tested_feature_count`` must
            each be None (or the literal ``'not-applicable'`` for the two
            string fields) exactly when ``analysis_role`` is
            ``'descriptive'`` (no inferential test was performed); a real
            declared family, BH-FDR | Benjamini-Hochberg | Bonferroni
            method, and positive tested-feature count are required exactly
            when the role is ``'inferential'``. Each entry of
            ``sensitivity_analyses`` is validated in full
            (:func:`validate_sensitivity_analysis`), and a record with a
            non-analysis role must not declare any. When declared, the
            structured ``replicate_unit_declaration`` and
            ``estimand_declaration`` are validated in full for
            non-degeneracy (:func:`validate_replicate_unit_declaration`,
            :func:`validate_estimand_declaration`); an inferential role
            additionally requires at least 2 independent biological
            replicates per declared unit, and conflicting embedded units
            are refused. When declared, the ``design_declaration`` is
            validated in full (:func:`validate_predeclared_design`); a
            non-analysis role must not declare one.
        TypeError: If ``record`` is not an :class:`AnalysisProvenance`.
    """
    if not isinstance(record, AnalysisProvenance):
        raise TypeError(
            "validate_analysis_provenance requires an AnalysisProvenance record; "
            "ad-hoc dictionaries cannot bind an analysis to its predeclared contract"
        )
    for provenance_field in (
        "analysis_id",
        "estimand",
        "replicate_unit",
        "null_model",
    ):
        _require_declared(getattr(record, provenance_field), provenance_field)

    if record.analysis_role not in _ALLOWED_ROLES:
        raise ProvenanceError(f"analysis_role must be one of {sorted(_ALLOWED_ROLES)}, got {record.analysis_role!r}")

    # Non-analysis states (plan section 8): a stopped/unavailable record
    # documents that a declared analysis did NOT run to results. It must not
    # carry multiplicity declarations, tested features, or analysis binding
    # fields that would imply results exist.
    if record.analysis_role in NON_ANALYSIS_ROLES:
        for field, value in (
            ("multiple_testing_family", record.multiple_testing_family),
            ("multiple_testing_method", record.multiple_testing_method),
            ("tested_feature_count", record.tested_feature_count),
            ("cohort_included_count", record.cohort_included_count),
            ("cohort_excluded_count", record.cohort_excluded_count),
            ("artifact_paths", record.artifact_paths),
            ("replicate_unit_declaration", record.replicate_unit_declaration),
            ("estimand_declaration", record.estimand_declaration),
            ("design_declaration", record.design_declaration),
        ):
            if value is not None:
                raise ProvenanceError(
                    f"analysis_role={record.analysis_role!r} records a halted or "
                    f"unavailable analysis, so {field}={value!r} must not be declared"
                )
        if record.sensitivity_analyses:
            raise ProvenanceError(
                f"analysis_role={record.analysis_role!r} records a halted or "
                f"unavailable analysis, so sensitivity_analyses="
                f"{record.sensitivity_analyses!r} must not be declared"
            )
        return

    if not isinstance(record.random_seed, int) or isinstance(record.random_seed, bool) or record.random_seed < 0:
        raise ProvenanceError(f"random_seed must be a non-negative integer, got {record.random_seed!r}")
    if (
        not isinstance(record.resampling_count, int)
        or isinstance(record.resampling_count, bool)
        or record.resampling_count < 1
    ):
        raise ProvenanceError(f"resampling_count must be a positive integer, got {record.resampling_count!r}")
    if record.tested_feature_count is None:
        if record.analysis_role != DESCRIPTIVE_ROLE:
            raise ProvenanceError(
                "analysis_role='inferential' must declare the number of tested "
                "features; got tested_feature_count=None"
            )
    elif record.analysis_role == DESCRIPTIVE_ROLE:
        raise ProvenanceError(
            f"analysis_role='descriptive' tests no features, so declaring "
            f"tested_feature_count={record.tested_feature_count!r} misstates the "
            f"analysis; use None (rendered 'not-applicable')"
        )
    elif (
        not isinstance(record.tested_feature_count, int)
        or isinstance(record.tested_feature_count, bool)
        or record.tested_feature_count < 1
    ):
        raise ProvenanceError(
            f"tested_feature_count must be a positive integer for inferential "
            f"analyses, got {record.tested_feature_count!r}"
        )

    # Multiplicity procedures belong to inferential tests only (plan section
    # 6): a descriptive lane applies no procedure, so recording one would
    # misstate the analysis. The symmetry is enforced fail-closed in BOTH
    # directions.
    method = record.multiple_testing_method
    declared_not_applicable = method is None or (isinstance(method, str) and method.strip().lower() == "not-applicable")
    if declared_not_applicable:
        if record.analysis_role != DESCRIPTIVE_ROLE:
            raise ProvenanceError(
                f"analysis_role='inferential' must declare a real multiple-testing "
                f"procedure; got multiple_testing_method={method!r}"
            )
    else:
        if record.analysis_role != INFERENTIAL_ROLE:
            raise ProvenanceError(
                f"analysis_role='descriptive' applies no multiple-testing procedure, so "
                f"declaring multiple_testing_method={method!r} misstates the analysis; use "
                f"None or 'not-applicable' for descriptive lanes"
            )
        declared = _require_declared(method, "multiple_testing_method").strip().lower()
        if declared not in _ALLOWED_MT_METHODS:
            raise ProvenanceError(
                f"multiple_testing_method must be one of {sorted(_ALLOWED_MT_METHODS)}, got "
                f"{record.multiple_testing_method!r}"
            )

    # The multiplicity FAMILY names the tested units of an inferential test
    # (plan section 6); a descriptive lane has none, so declaring one would
    # misstate the analysis. Same fail-closed symmetry as the method.
    family = record.multiple_testing_family
    family_not_applicable = family is None or (isinstance(family, str) and family.strip().lower() == "not-applicable")
    if family_not_applicable:
        if record.analysis_role != DESCRIPTIVE_ROLE:
            raise ProvenanceError(
                f"analysis_role='inferential' must declare a real multiple-testing "
                f"family; got multiple_testing_family={family!r}"
            )
    else:
        if record.analysis_role != INFERENTIAL_ROLE:
            raise ProvenanceError(
                f"analysis_role='descriptive' performs no inferential test, so "
                f"declaring multiple_testing_family={family!r} misstates the analysis; "
                f"use None or 'not-applicable' for descriptive lanes"
            )
        _require_declared(family, "multiple_testing_family")

    if not isinstance(record.software_versions, Mapping) or not record.software_versions:
        raise ProvenanceError("software_versions must be a non-empty mapping of package name to version")
    for name, version in record.software_versions.items():
        _require_declared(name, "software_versions key")
        _require_declared(version, f"software_versions[{name!r}]")
    # Reporting-contract bindings (plan section 9) and the harmonization /
    # species-tree flags are optional: when declared they must be real
    # values, never placeholders. Denominator consistency is checked
    # structurally (non-negative, integers).
    for field in (
        "data_root_snapshot_id",
        "metadata_harmonization_review",
        "species_tree_source",
        "species_tree_branch_length_scale",
    ):
        value = getattr(record, field)
        if value is not None:
            _require_declared(value, field)
    for field in ("cohort_included_count", "cohort_excluded_count"):
        value = getattr(record, field)
        if value is not None and (not isinstance(value, int) or isinstance(value, bool) or value < 0):
            raise ProvenanceError(f"{field} must be a non-negative integer when declared, got {value!r}")
    if record.artifact_paths is not None:
        if not isinstance(record.artifact_paths, Mapping) or not record.artifact_paths:
            raise ProvenanceError("artifact_paths must be a non-empty mapping when declared")
        for name, path in record.artifact_paths.items():
            _require_declared(name, "artifact_paths key")
            _require_declared(path, f"artifact_paths[{name!r}]")
    if not isinstance(record.sensitivity_analyses, tuple):
        raise ProvenanceError(
            "sensitivity_analyses must be a tuple of SensitivityAnalysis records, "
            f"got {type(record.sensitivity_analyses).__name__}"
        )
    for sensitivity in record.sensitivity_analyses:
        validate_sensitivity_analysis(sensitivity)
    # Biological-replicate unit and estimand declarations (plan sections 1
    # and 3): optional structured refinements validated fail-closed for
    # non-degeneracy. An inferential role additionally demands an
    # inferential-grade replicate unit (at least 2 independent biological
    # replicates per declared unit).
    if record.replicate_unit_declaration is not None:
        validate_replicate_unit_declaration(
            record.replicate_unit_declaration,
            require_inferential_grade=record.analysis_role == INFERENTIAL_ROLE,
        )
    if record.estimand_declaration is not None:
        validate_estimand_declaration(record.estimand_declaration)
        embedded = record.estimand_declaration.replicate_unit
        if (
            record.replicate_unit_declaration is not None
            and embedded is not None
            and embedded != record.replicate_unit_declaration
        ):
            raise ProvenanceError(
                "conflicting replicate-unit declarations: the provenance record and the "
                f"estimand declaration name different units ({record.replicate_unit_declaration!r} "
                f"vs {embedded!r})"
            )
    if record.design_declaration is not None:
        validate_predeclared_design(record.design_declaration)


def render_analysis_provenance_block(record: AnalysisProvenance) -> list[str]:
    """Render the provenance record as additive ``key: value`` summary lines.

    The lines use the same ``key: value`` shape as the project's
    ``analysis_summary.txt`` (written by ``run_cross_species_analysis.py``),
    so they can be appended without changing existing keys. Validation runs
    first: a record that would render placeholder provenance fails closed.
    Descriptive lanes render ``not-applicable`` for the multiplicity
    fields; optional reporting bindings render only when declared (artifact
    paths as ``analysis_provenance_artifact_<name>``). Structured
    replicate-unit and estimand declarations render only when declared, as
    ``analysis_provenance_replicate_unit_*`` (name, nesting, declared
    minimum, and per-stratum counts, strata sorted) and
    ``analysis_provenance_estimand_*`` lines. Each registered
    :class:`SensitivityAnalysis` renders as
    ``analysis_provenance_sensitivity_<index>_<field>`` lines, with
    ``notes`` only when non-empty. A declared ``design_declaration``
    renders additively as ``analysis_provenance_design_covariate_<name>``
    lines (covariates sorted, levels in declared order).

    Returns:
        Deterministic list of ``analysis_provenance_*`` lines.
    """
    validate_analysis_provenance(record)
    versions = "; ".join(f"{name}={record.software_versions[name]}" for name in sorted(record.software_versions))
    lines = [
        f"analysis_provenance_role: {record.analysis_role}",
        f"analysis_provenance_analysis_id: {record.analysis_id}",
        f"analysis_provenance_estimand: {record.estimand}",
        f"analysis_provenance_replicate_unit: {record.replicate_unit}",
        f"analysis_provenance_random_seed: {record.random_seed}",
        f"analysis_provenance_resampling_count: {record.resampling_count}",
        f"analysis_provenance_null_model: {record.null_model}",
        "analysis_provenance_multiple_testing_family: "
        + ("not-applicable" if record.multiple_testing_family is None else record.multiple_testing_family),
        "analysis_provenance_multiple_testing_method: "
        + (
            "not-applicable"
            if record.multiple_testing_method is None
            else record.multiple_testing_method.strip().lower()
        ),
        "analysis_provenance_tested_feature_count: "
        + ("not-applicable" if record.tested_feature_count is None else str(record.tested_feature_count)),
    ]
    lines.extend(
        f"analysis_provenance_{name}: {value}"
        for name, value in (
            ("data_root_snapshot_id", record.data_root_snapshot_id),
            ("cohort_included_count", record.cohort_included_count),
            ("cohort_excluded_count", record.cohort_excluded_count),
            ("metadata_harmonization_review", record.metadata_harmonization_review),
            ("species_tree_source", record.species_tree_source),
            ("species_tree_branch_length_scale", record.species_tree_branch_length_scale),
        )
        if value is not None
    )
    if record.artifact_paths is not None:
        lines.extend(
            f"analysis_provenance_artifact_{name}: {record.artifact_paths[name]}"
            for name in sorted(record.artifact_paths)
        )
    if record.replicate_unit_declaration is not None:
        unit = record.replicate_unit_declaration
        lines.extend(
            f"analysis_provenance_replicate_unit_{key}: {value}"
            for key, value in (
                ("name", unit.name),
                ("nesting", unit.nesting),
                ("min_independent_replicates", unit.min_independent_replicates),
            )
        )
        lines.extend(
            f"analysis_provenance_replicate_unit_count_{stratum}: {unit.counts[stratum]}"
            for stratum in sorted(unit.counts)
        )
    if record.estimand_declaration is not None:
        lines.extend(
            f"analysis_provenance_estimand_{key}: {value}"
            for key, value in (
                ("name", record.estimand_declaration.name),
                ("contrast", record.estimand_declaration.contrast),
                ("permitted_interpretation", record.estimand_declaration.permitted_interpretation),
            )
        )
    if record.design_declaration is not None:
        lines.extend(
            f"analysis_provenance_design_covariate_{covariate}: "
            + ",".join(record.design_declaration.covariate_strata[covariate])
            for covariate in sorted(record.design_declaration.covariate_strata)
        )
    for index, sensitivity in enumerate(record.sensitivity_analyses, start=1):
        lines.extend(
            f"analysis_provenance_sensitivity_{index}_{field}: {value}"
            for field, value in (
                ("name", sensitivity.name),
                ("varied_parameter", sensitivity.varied_parameter),
                ("baseline_value", sensitivity.baseline_value),
                ("varied_values", ",".join(sensitivity.varied_values)),
                ("expected_direction", sensitivity.expected_direction),
            )
        )
        if sensitivity.notes:
            lines.append(f"analysis_provenance_sensitivity_{index}_notes: {sensitivity.notes}")
    lines.append(f"analysis_provenance_software_versions: {versions}")
    return lines


def result_role(result: Any) -> str:
    """Return the declared role of a result object, failing closed.

    Descriptive results produced by ``cross_species.compute_fingerprint_*``
    carry ``attrs["role"] == "descriptive"``. The recognized role set also
    covers ``"inferential"`` and the non-analysis states ``"stopped"`` and
    ``"unavailable"``. Any result without an explicit role marker is refused
    rather than being silently treated as either descriptive or inferential.
    """
    attrs = getattr(result, "attrs", None)
    if not isinstance(attrs, Mapping) or "role" not in attrs:
        raise StatisticsContractError("result carries no declared role; refusing to classify unlabeled output")
    role = attrs["role"]
    if role not in _ALLOWED_ROLES:
        raise StatisticsContractError(f"result declares unknown role {role!r}")
    return str(role)


def benjamini_hochberg_fdr(p_values: Sequence[float]) -> list[float]:
    """Apply the Benjamini-Hochberg FDR procedure to raw p-values.

    Args:
        p_values: Non-empty sequence of raw p-values in [0, 1].

    Returns:
        Adjusted q-values in the same order as the input.

    Raises:
        StatisticsContractError: If the sequence is empty or any value is
            not a finite p-value in [0, 1].
    """
    values = list(p_values)
    if not values:
        raise StatisticsContractError("benjamini_hochberg_fdr requires at least one p-value")
    for index, value in enumerate(values):
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise StatisticsContractError(f"p-values must be numeric; index {index} got {value!r}")
        if value != value or value in (float("inf"), float("-inf")):
            raise StatisticsContractError(f"p-value at index {index} is not finite")
        if value < 0.0 or value > 1.0:
            raise StatisticsContractError(f"p-value at index {index} is outside [0, 1]: {value!r}")

    n = len(values)
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    adjusted = [0.0] * n
    running_min = 1.0
    for rank_index in range(n - 1, -1, -1):
        original_index, p_value = indexed[rank_index]
        rank = rank_index + 1
        running_min = min(running_min, p_value * n / rank)
        adjusted[original_index] = min(running_min, 1.0)
    return adjusted


def declared_inferential_bh_fdr(
    p_values: Sequence[float],
    contract: AnalysisProvenance,
    evidence_manifest_frozen: bool = False,
) -> dict[str, Any]:
    """GATED: apply the contract's declared BH-FDR procedure to raw p-values.

    POST-FREEZE USE ONLY. Descriptive permutation and sensitivity scores must
    never reach this function; inferential adjustment requires (a) an explicit
    ``evidence_manifest_frozen=True`` affirmation at every call site, and
    (b) a validated :class:`AnalysisProvenance` declared
    ``analysis_role="inferential"`` with a BH-FDR procedure whose
    ``tested_feature_count`` matches the supplied family exactly. Default
    ``evidence_manifest_frozen=False`` always refuses.

    Returns:
        Dict with ``role: "inferential"``, the declared method, raw and
        BH-adjusted p-values, and the contract fields echoed for provenance.

    Raises:
        RuntimeError: If ``evidence_manifest_frozen`` is not True.
        ProvenanceError: If the contract itself fails validation.
        StatisticsContractError: If the contract is not declared inferential
            with BH-FDR, or the tested-feature count does not match.
    """
    if not evidence_manifest_frozen:
        raise RuntimeError(
            "declared_inferential_bh_fdr is gated for post-freeze use: "
            "refusing to run while the evidence manifest is unfrozen"
        )
    validate_analysis_provenance(contract)
    if contract.analysis_role != INFERENTIAL_ROLE:
        raise StatisticsContractError(
            f"declared_inferential_bh_fdr requires a contract declared "
            f"analysis_role={INFERENTIAL_ROLE!r}, got {contract.analysis_role!r}"
        )
    declared_method = contract.multiple_testing_method
    if declared_method is None or declared_method.strip().lower() not in {"bh-fdr", "benjamini-hochberg"}:
        raise StatisticsContractError(
            f"declared_inferential_bh_fdr requires a BH-FDR procedure, got " f"{contract.multiple_testing_method!r}"
        )

    values = [float(value) for value in p_values]
    if contract.tested_feature_count != len(values):
        raise StatisticsContractError(
            f"contract declares tested_feature_count={contract.tested_feature_count} but "
            f"{len(values)} p-values were supplied; refusing to adjust a different family"
        )
    adjusted = benjamini_hochberg_fdr(values)
    return {
        "role": INFERENTIAL_ROLE,
        "analysis_id": contract.analysis_id,
        "estimand": contract.estimand,
        "replicate_unit": contract.replicate_unit,
        "random_seed": contract.random_seed,
        "null_model": contract.null_model,
        "multiple_testing_family": contract.multiple_testing_family,
        "multiple_testing_method": "bh-fdr",
        "tested_feature_count": len(values),
        "raw_p_values": values,
        "adjusted_p_values": adjusted,
        "gate": "post-freeze",
    }


# =============================================================================
# Paired effect sizes, multi-study heterogeneity, and leave-one-study-out
# sensitivity (MJ-01 inferential lane)
# =============================================================================


_MASK64 = (1 << 64) - 1
_SPLITMIX_GAMMA = 0x9E3779B97F4A7C15
_SPLITMIX_M1 = 0xBF58476DCE4E6B25
_SPLITMIX_M2 = 0x94D049BB133111EB
_Z_975 = 1.959963984540054  # two-sided normal quantile for a 95% analytic CI


def _splitmix64_picks(seed: int, count: int, size: int) -> list[list[int]]:
    """Deterministic splitmix64 resample indices (no external RNG dependency)."""
    state = seed & _MASK64
    picks: list[list[int]] = []
    for _ in range(count):
        replicate: list[int] = []
        for _ in range(size):
            state = (state + _SPLITMIX_GAMMA) & _MASK64
            z = state
            z = ((z ^ (z >> 30)) * _SPLITMIX_M1) & _MASK64
            z = ((z ^ (z >> 27)) * _SPLITMIX_M2) & _MASK64
            replicate.append((z ^ (z >> 31)) % size)
        picks.append(replicate)
    return picks


def _percentile(sorted_values: Sequence[float], percent: float) -> float:
    """Linearly interpolated percentile of a pre-sorted sequence (numpy-compatible)."""
    position = percent / 100.0 * (len(sorted_values) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(sorted_values[lower])
    fraction = position - lower
    return float(sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction)


def _paired_log_ratios(reference_values: Sequence[float], treatment_values: Sequence[float]) -> list[float]:
    """Validate paired observations and return their log2(treatment/reference) ratios.

    Raises:
        StatisticsContractError: On empty or mismatched pairs, non-numeric
            or non-finite values, or non-positive abundances (a log2
            fold-change is undefined at or below zero).
    """
    reference_list = list(reference_values)
    treatment_list = list(treatment_values)
    if not reference_list or len(reference_list) != len(treatment_list):
        raise StatisticsContractError(
            "paired effect sizes require a non-empty, aligned pair of observation "
            f"sequences; got {len(reference_list)} reference vs "
            f"{len(treatment_list)} treatment"
        )
    ratios: list[float] = []
    for index, (reference, treatment) in enumerate(zip(reference_list, treatment_list)):
        for label, value in (("reference", reference), ("treatment", treatment)):
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise StatisticsContractError(
                    f"paired observations must be numeric; index {index} ({label}) got {value!r}"
                )
            if value != value or value in (float("inf"), float("-inf")):
                raise StatisticsContractError(f"paired observation at index {index} ({label}) is not finite")
            if value <= 0:
                raise StatisticsContractError(
                    f"paired observation at index {index} ({label}) is {value!r}; "
                    "log2 fold-change requires strictly positive values"
                )
        ratios.append(math.log2(treatment_list[index] / reference_list[index]))
    return ratios


def paired_log2fc_effect(
    reference_values: Sequence[float],
    treatment_values: Sequence[float],
    *,
    bootstrap_replicates: int = 0,
    random_seed: int = 0,
) -> dict[str, Any]:
    """Paired log2 fold-change with an analytic SE and optional bootstrap CI.

    Contract: the point effect is the mean paired
    ``log2(treatment/reference)`` over aligned pairs; ``se`` is the
    analytic standard error of that mean (``None`` for a single pair,
    which carries no uncertainty); with ``bootstrap_replicates >= 1`` the
    pairs are resampled with replacement through a deterministic
    splitmix64 stream seeded by ``random_seed`` and the 2.5/97.5
    percentiles of the replicate means are returned. Pure computation, no
    I/O, no gating — gating is the caller's contract decision (see
    :func:`declared_effect_size`).

    Raises:
        StatisticsContractError: On empty/mismatched pairs, non-numeric,
            non-finite, or non-positive values, a negative seed, or a
            negative bootstrap replicate count.
    """
    ratios = _paired_log_ratios(reference_values, treatment_values)
    if not isinstance(bootstrap_replicates, int) or isinstance(bootstrap_replicates, bool) or bootstrap_replicates < 0:
        raise StatisticsContractError(
            f"bootstrap_replicates must be a non-negative integer, got {bootstrap_replicates!r}"
        )
    if not isinstance(random_seed, int) or isinstance(random_seed, bool) or random_seed < 0:
        raise StatisticsContractError(f"random_seed must be a non-negative integer, got {random_seed!r}")
    n_pairs = len(ratios)
    mean = sum(ratios) / n_pairs
    se = None
    if n_pairs > 1:
        variance = sum((value - mean) ** 2 for value in ratios) / (n_pairs - 1)
        se = math.sqrt(variance / n_pairs)
    record: dict[str, Any] = {
        "effect": mean,
        "se": se,
        "n_pairs": n_pairs,
        "random_seed": random_seed,
    }
    if bootstrap_replicates > 0:
        replicate_means: list[float] = []
        for picks in _splitmix64_picks(random_seed, bootstrap_replicates, n_pairs):
            replicate_means.append(sum(ratios[pick] for pick in picks) / n_pairs)
        replicate_means.sort()
        record["ci_low"] = _percentile(replicate_means, 2.5)
        record["ci_high"] = _percentile(replicate_means, 97.5)
        record["bootstrap_n_success"] = bootstrap_replicates
        record["bootstrap_replicates"] = bootstrap_replicates
    return record


def declared_effect_size(
    reference_values: Sequence[float],
    treatment_values: Sequence[float],
    contract: AnalysisProvenance,
    evidence_manifest_frozen: bool = False,
) -> dict[str, Any]:
    """Role-conditional paired effect-size record under a predeclared contract.

    Role-conditional gating, mirroring :func:`declared_inferential_bh_fdr`:

    - a validated ``descriptive`` contract produces a point-estimate record
      with every inferential field (method, CI, multiplicity, gate) rendered
      as the literal ``"not-applicable"``; no frozen-manifest affirmation is
      demanded because no inference is performed;
    - an ``inferential`` contract is GATED: it requires an explicit
      ``evidence_manifest_frozen=True`` affirmation (default always
      refuses) and produces a full record carrying the effect-size method,
      the paired replicate count (plus per-stratum counts when the contract
      declares a replicate unit), and the declared multiplicity
      family/method provenance. The bootstrap CI uses the contract's own
      ``random_seed`` and ``resampling_count``;
    - a ``stopped``/``unavailable`` contract produces no record at all.

    Raises:
        ProvenanceError: If the contract fails validation.
        StatisticsContractError: On invalid pairs, a non-analysis role, or
            an inferential contract with a single paired replicate.
        RuntimeError: If an inferential contract is used while the evidence
            manifest is not affirmed frozen.
    """
    validate_analysis_provenance(contract)
    ratios = _paired_log_ratios(reference_values, treatment_values)
    if contract.analysis_role in NON_ANALYSIS_ROLES:
        raise StatisticsContractError(
            f"analysis_role={contract.analysis_role!r} records a halted or unavailable "
            "analysis; no effect-size record may be produced"
        )
    mean = sum(ratios) / len(ratios)
    if contract.analysis_role == DESCRIPTIVE_ROLE:
        analytic_se: Any = "not-applicable"
        if len(ratios) > 1:
            variance = sum((value - mean) ** 2 for value in ratios) / (len(ratios) - 1)
            analytic_se = math.sqrt(variance / len(ratios))
        return {
            "role": DESCRIPTIVE_ROLE,
            "analysis_id": contract.analysis_id,
            "estimand": contract.estimand,
            "effect": mean,
            "se": analytic_se,
            "effect_size_method": "not-applicable",
            "ci_low": "not-applicable",
            "ci_high": "not-applicable",
            "bootstrap_resampling_count": "not-applicable",
            "replicate_count": len(ratios),
            "multiple_testing_family": "not-applicable",
            "multiple_testing_method": "not-applicable",
            "gate": "not-applicable",
        }
    if not evidence_manifest_frozen:
        raise RuntimeError(
            "declared_effect_size is gated for inferential contracts: refusing to "
            "run while the evidence manifest is unfrozen"
        )
    if len(ratios) < 2:
        raise StatisticsContractError(
            "an inferential effect-size record requires at least 2 paired " f"replicates; got {len(ratios)}"
        )
    bootstrap = paired_log2fc_effect(
        reference_values,
        treatment_values,
        bootstrap_replicates=contract.resampling_count,
        random_seed=contract.random_seed,
    )
    replicate_counts: Any = "not-applicable"
    if contract.replicate_unit_declaration is not None:
        replicate_counts = {
            stratum: contract.replicate_unit_declaration.counts[stratum]
            for stratum in sorted(contract.replicate_unit_declaration.counts)
        }
    return {
        "role": INFERENTIAL_ROLE,
        "gate": "post-freeze",
        "analysis_id": contract.analysis_id,
        "estimand": contract.estimand,
        "replicate_unit": contract.replicate_unit,
        "effect_size_method": "paired-log2fc-mean-bootstrap-percentile",
        "effect": bootstrap["effect"],
        "se": bootstrap["se"],
        "ci_low": bootstrap["ci_low"],
        "ci_high": bootstrap["ci_high"],
        "bootstrap_resampling_count": bootstrap["bootstrap_replicates"],
        "bootstrap_n_success": bootstrap["bootstrap_n_success"],
        "random_seed": contract.random_seed,
        "replicate_count": len(ratios),
        "replicate_unit_counts": replicate_counts,
        "multiple_testing_family": contract.multiple_testing_family,
        "multiple_testing_method": contract.multiple_testing_method,
    }


def render_effect_size_record(record: Mapping[str, Any]) -> list[str]:
    """Render an effect-size record as additive ``effect_size_*`` lines.

    Role-conditional rendering: the record's own fields are echoed
    verbatim, so a descriptive record renders its inferential fields
    exactly as ``not-applicable`` while an inferential record carries the
    method, replicate count, and multiplicity family/method provenance.
    Rendering is deterministic (keys sorted).

    Raises:
        StatisticsContractError: If the record is not a mapping, carries no
            ``role`` key, or declares an unknown role.
    """
    if not isinstance(record, Mapping) or "role" not in record:
        raise StatisticsContractError(
            "effect-size record carries no declared role; refusing to render " "unlabeled output"
        )
    role = record["role"]
    if role not in _ALLOWED_ROLES:
        raise StatisticsContractError(f"effect-size record declares unknown role {role!r}")
    return [f"effect_size_{key}: {record[key]}" for key in sorted(record)]


def _heterogeneity_components(effects: Sequence[float], standard_errors: Sequence[float]) -> dict[str, Any]:
    """Cochran's Q, tau-squared, and I-squared for aligned study effects.

    Raises:
        HeterogeneityError: On empty/mismatched input, non-finite effects,
            or non-positive/non-finite standard errors.
    """
    effect_values = [float(value) for value in effects]
    se_values = [float(value) for value in standard_errors]
    if not effect_values or len(effect_values) != len(se_values):
        raise HeterogeneityError(
            f"effects and standard_errors must be non-empty aligned sequences, got "
            f"{len(effect_values)} effects vs {len(se_values)} standard errors"
        )
    for index, (effect, se) in enumerate(zip(effect_values, se_values)):
        if effect != effect or effect in (float("inf"), float("-inf")):
            raise HeterogeneityError(f"effect at index {index} is not finite")
        if se != se or se in (float("inf"), float("-inf")):
            raise HeterogeneityError(f"standard error at index {index} is not finite")
        if se <= 0:
            raise HeterogeneityError(f"standard error at index {index} must be strictly positive, got {se!r}")
    weights = [1.0 / se**2 for se in se_values]
    weight_sum = sum(weights)
    fixed_effect = sum(weight * effect for weight, effect in zip(weights, effect_values)) / weight_sum
    q_statistic = sum(weight * (effect - fixed_effect) ** 2 for weight, effect in zip(weights, effect_values))
    df = len(effect_values) - 1
    c_statistic = weight_sum - sum(weight**2 for weight in weights) / weight_sum
    tau_squared = max(0.0, (q_statistic - df) / c_statistic)
    i_squared = max(0.0, (q_statistic - df) / q_statistic) * 100.0 if q_statistic > 0 else 0.0
    return {
        "fixed_effect": fixed_effect,
        "q": q_statistic,
        "df": df,
        "tau_squared": tau_squared,
        "i_squared_percent": i_squared,
    }


def _dl_combine(effects: Sequence[float], standard_errors: Sequence[float]) -> tuple[float, float]:
    """DerSimonian-Laird random-effects combination (pure math mirror).

    Returns ``(combined_effect, combined_se)``; requires at least 2 studies.

    Raises:
        HeterogeneityError: On invalid or too-thin input.
    """
    effect_values = [float(value) for value in effects]
    if len(effect_values) < 2:
        raise HeterogeneityError(f"random-effects combination requires at least 2 studies, got {len(effect_values)}")
    components = _heterogeneity_components(effect_values, [float(value) for value in standard_errors])
    random_weights = [1.0 / (float(se) ** 2 + components["tau_squared"]) for se in standard_errors]
    combined = sum(weight * effect for weight, effect in zip(random_weights, effect_values)) / sum(random_weights)
    combined_se = math.sqrt(1.0 / sum(random_weights))
    return combined, combined_se


def heterogeneity_record(
    effects: Sequence[float],
    standard_errors: Sequence[float],
    contract: AnalysisProvenance,
) -> dict[str, Any]:
    """Cochran-Q / I-squared heterogeneity record for multi-study strata.

    Role-conditional: an ``inferential`` contract requires at least 2
    studies and fails closed otherwise (a single-study stratum carries no
    between-study heterogeneity and is refused, never approximated); a
    ``descriptive`` contract is exempt and returns a record whose
    heterogeneity fields render as ``not-applicable``; a
    ``stopped``/``unavailable`` contract produces no record at all.

    Raises:
        ProvenanceError: If the contract fails validation.
        StatisticsContractError: On a non-analysis contract.
        HeterogeneityError: On fewer than 2 studies for an inferential
            contract, or invalid effects/standard errors.
    """
    validate_analysis_provenance(contract)
    if contract.analysis_role in NON_ANALYSIS_ROLES:
        raise StatisticsContractError(
            f"analysis_role={contract.analysis_role!r} records a halted or unavailable "
            "analysis; no heterogeneity record may be produced"
        )
    if contract.analysis_role == DESCRIPTIVE_ROLE:
        return {
            "role": DESCRIPTIVE_ROLE,
            "analysis_id": contract.analysis_id,
            "heterogeneity_status": "exempt",
            "cochran_q": "not-applicable",
            "df": "not-applicable",
            "i_squared_percent": "not-applicable",
            "tau_squared": "not-applicable",
            "n_studies": len(list(effects)),
        }
    if len(list(effects)) < 2:
        raise HeterogeneityError(
            f"heterogeneity requires at least 2 studies, got {len(list(effects))}; a "
            "single-study stratum carries no between-study heterogeneity and is "
            "refused rather than silently treated as homogeneous"
        )
    components = _heterogeneity_components(effects, standard_errors)
    return {
        "role": INFERENTIAL_ROLE,
        "analysis_id": contract.analysis_id,
        "heterogeneity_status": "computed",
        "cochran_q": components["q"],
        "df": components["df"],
        "i_squared_percent": components["i_squared_percent"],
        "tau_squared": components["tau_squared"],
        "n_studies": len(list(effects)),
    }


def leave_one_study_out_deltas(studies: Mapping[Any, tuple[float, float]]) -> dict[str, Any]:
    """Leave-one-study-out recompute: per-exclusion combined-effect deltas.

    Pure computation, no I/O: for each study, the DerSimonian-Laird
    combined effect is recomputed from the remaining studies and the delta
    against the full-data combined effect is returned. Requires at least 3
    studies so every exclusion retains at least 2 — the minimum for a
    random-effects recombination.

    Args:
        studies: Mapping of study label to ``(effect, standard_error)``.

    Returns:
        ``{"full_effect": float, "deltas": {label: float}}`` with labels in
        sorted (string) order.

    Raises:
        HeterogeneityError: On fewer than 3 studies, a non-mapping input,
            or non-finite effects / non-positive standard errors.
    """
    if not isinstance(studies, Mapping):
        raise HeterogeneityError(
            f"leave_one_study_out_deltas requires a mapping of study label to "
            f"(effect, standard_error), got {type(studies).__name__}"
        )
    if len(studies) < 3:
        raise HeterogeneityError(
            f"leave-one-study-out recombination requires at least 3 studies so each "
            f"exclusion retains at least 2; got {len(studies)}"
        )
    ordered = sorted(studies.items(), key=lambda item: str(item[0]))
    full_effect, _ = _dl_combine(
        [pair[0] for _, pair in ordered],
        [pair[1] for _, pair in ordered],
    )
    deltas: dict[Any, float] = {}
    for label, _ in ordered:
        remaining = [(effect, se) for other, (effect, se) in ordered if other != label]
        combined, _ = _dl_combine(
            [effect for effect, _ in remaining],
            [se for _, se in remaining],
        )
        deltas[label] = combined - full_effect
    return {"full_effect": full_effect, "deltas": deltas}


# =============================================================================
# Orthology and species-tree invariants (plan sections 4 and 8)
# =============================================================================


def validate_orthology_profile_invariants(
    presence: pd.DataFrame,
    *,
    min_species_fraction: float = 0.5,
    min_species_per_orthogroup: int = 2,
) -> None:
    """Fail closed on orthology x species presence-table violations.

    The presence table has orthogroups as rows and species as columns; a
    truthy value means the orthogroup is present (mapped) in that species.
    Enforced invariants:

    - unique, non-placeholder orthogroup and species labels;
    - presence values restricted to boolean or numeric 0/1 — missing or
      ambiguous states must be explicit, never silently coerced (section 8);
    - every orthogroup maps to at least
      ``max(min_species_per_orthogroup, ceil(min_species_fraction * n_species))``
      species; violations are reported together with explicit reason codes.

    Raises:
        OrthologyInvariantError: On any violated invariant.
    """
    if not isinstance(presence, pd.DataFrame):
        raise OrthologyInvariantError("presence table must be a pandas DataFrame")
    if presence.empty:
        raise OrthologyInvariantError("presence table is empty")
    if presence.shape[1] < 2:
        raise OrthologyInvariantError(f"presence table must cover at least 2 species, got {presence.shape[1]}")
    if not 0.0 <= min_species_fraction <= 1.0:
        raise OrthologyInvariantError(f"min_species_fraction must be within [0, 1], got {min_species_fraction!r}")
    if (
        not isinstance(min_species_per_orthogroup, int)
        or isinstance(min_species_per_orthogroup, bool)
        or min_species_per_orthogroup < 1
    ):
        raise OrthologyInvariantError(
            f"min_species_per_orthogroup must be a positive integer, got " f"{min_species_per_orthogroup!r}"
        )

    duplicated_orthogroups = sorted({label for label in presence.index[presence.index.duplicated(keep=False)].tolist()})
    if duplicated_orthogroups:
        raise OrthologyInvariantError(f"duplicate orthogroup labels: {duplicated_orthogroups}")
    duplicated_species = sorted({label for label in presence.columns[presence.columns.duplicated(keep=False)].tolist()})
    if duplicated_species:
        raise OrthologyInvariantError(f"duplicate species labels: {duplicated_species}")

    bad_labels = [
        label
        for label in list(presence.index) + list(presence.columns)
        if not isinstance(label, str) or label.strip() == "" or label.strip().lower() in _PLACEHOLDER_STRINGS
    ]
    if bad_labels:
        raise OrthologyInvariantError(f"presence table contains missing or placeholder labels: {bad_labels}")

    non_numeric = [
        str(column)
        for column, dtype in presence.dtypes.items()
        if not pd.api.types.is_bool_dtype(dtype) and not pd.api.types.is_numeric_dtype(dtype)
    ]
    if non_numeric:
        raise OrthologyInvariantError(
            "presence values must be boolean or numeric 0/1; non-numeric columns: " + ", ".join(non_numeric)
        )
    numeric = presence.astype(float)
    if numeric.isna().any().any():
        missing = numeric.isna().any(axis=1)
        offenders = [str(label) for label in numeric.index[missing].tolist()]
        raise OrthologyInvariantError(
            "presence table contains missing values; absence must be encoded explicitly as 0, "
            "offending orthogroups: " + ", ".join(offenders)
        )
    unexpected_mask = (numeric != 0.0) & (numeric != 1.0)
    if unexpected_mask.any().any():
        offenders = [str(label) for label in numeric.index[unexpected_mask.any(axis=1)].tolist()]
        raise OrthologyInvariantError(
            "presence values outside {0, 1} are ambiguous and refused, offending orthogroups: " + ", ".join(offenders)
        )

    n_species = numeric.shape[1]
    required = max(min_species_per_orthogroup, math.ceil(min_species_fraction * n_species))
    counts = numeric.sum(axis=1)
    insufficient = counts[counts < required]
    if not insufficient.empty:
        detail = ", ".join(f"{label} ({int(count)}/{n_species} species)" for label, count in insufficient.items())
        raise OrthologyInvariantError(
            f"orthology bridge below declared coverage threshold "
            f"(required {required} of {n_species} species per orthogroup): {detail}"
        )


def _parse_newick(newick: str) -> dict[str, Any]:
    """Parse a plain Newick string (names, branch lengths) into nested dicts.

    Supports unquoted names, optional branch lengths, and optional internal
    labels. Any parse failure raises :class:`TreeInvariantError`.
    """
    text = newick.strip()
    if not text.endswith(";"):
        raise TreeInvariantError("Newick string must terminate with ';'")
    text = text[:-1]
    position = 0

    def parse_node() -> dict[str, Any]:
        nonlocal position
        children: list[dict[str, Any]] = []
        if position < len(text) and text[position] == "(":
            position += 1
            children.append(parse_node())
            while position < len(text) and text[position] == ",":
                position += 1
                children.append(parse_node())
            if position >= len(text) or text[position] != ")":
                raise TreeInvariantError(f"unbalanced parentheses at character {position}")
            position += 1
        start = position
        while position < len(text) and text[position] not in "(),:":
            position += 1
        name = text[start:position].strip()
        if position < len(text) and text[position] == ":":
            position += 1
            length_start = position
            while position < len(text) and text[position] not in "(),":
                position += 1
            length_text = text[length_start:position].strip()
            try:
                float(length_text)
            except ValueError:
                raise TreeInvariantError(f"invalid branch length {length_text!r} in Newick string") from None
        return {"name": name, "children": children}

    root = parse_node()
    if position != len(text):
        raise TreeInvariantError(f"unexpected characters after Newick root at position {position}")
    return root


def validate_species_tree_invariants(
    tree: str | Mapping[str, Any],
    *,
    rooted: bool | None = None,
    require_rooted: bool = True,
    require_bifurcating_root: bool = False,
) -> None:
    """Fail closed on malformed, unlabeled, or undeclared species trees.

    Accepts either a Newick string or the nested-dict representation used by
    ``cross_species.phylogenetic_expression_profile`` (keys ``name``,
    ``children``, ``distance``).

    Rootedness is CALLER-DECLARED PROVENANCE (statistical_analysis_plan.md
    section 5.3 requires the tree source and its provenance to be recorded).
    Plain Newick topology cannot establish biological rootedness: an
    unrooted tree can be serialized behind a bifurcating root, and a
    polytomous root is a valid hard polytomy in a rooted tree. Therefore the
    number of root children is never used to infer rootedness. When
    ``require_rooted`` is set, the caller must declare rootedness through
    ``rooted``:

    - ``rooted=None`` (default) fails closed with :class:`ProvenanceError`:
      an explicit rootedness declaration is required;
    - ``rooted=False`` raises :class:`TreeInvariantError`;
    - ``rooted=True`` passes.

    Enforced invariants:

    - every node is well-formed (leaves and, for dict input, internal nodes
      and the root carry non-placeholder names; internal names may be empty
      only in Newick input);
    - leaf labels are unique (duplicated taxa fail closed);
    - at least two leaves;
    - declared rootedness per ``require_rooted``;
    - ``require_bifurcating_root=True`` additionally asserts the structural
      property that the root node has exactly 2 children (a naming of what
      is actually checked; it does not prove biological rootedness).

    Raises:
        TreeInvariantError: On any violated structural invariant, malformed
            input, or a declared-unrooted tree when ``require_rooted`` is set.
        ProvenanceError: If ``require_rooted`` is set and ``rooted`` was not
            explicitly declared by the caller.
    """
    root: Mapping[str, Any]
    if isinstance(tree, str):
        root = _parse_newick(tree)
        internal_names_optional = True
    elif isinstance(tree, Mapping):
        root = tree
        internal_names_optional = False
    else:
        raise TreeInvariantError("species tree must be a Newick string or a nested-dict tree mapping")

    leaves: list[str] = []
    seen_leaves: dict[str, int] = {}
    root_children: list[Any] | None = None

    def visit(node: Any, depth: int, *, is_root: bool) -> None:
        nonlocal root_children
        if not isinstance(node, Mapping) or "name" not in node:
            raise TreeInvariantError("every tree node must be a mapping with a 'name' key")
        name = node["name"]
        if not isinstance(name, str):
            raise TreeInvariantError("tree node labels must be strings")
        label = name.strip()
        # Dict trees must name every node (phylogenetic_expression_profile
        # reads node["name"]); Newick legitimately allows unnamed internal
        # nodes and roots, so only leaves are held to the named standard.
        if not internal_names_optional and (label == "" or label.lower() in _PLACEHOLDER_STRINGS):
            raise TreeInvariantError(
                f"tree {'root' if is_root else 'internal'} label is missing or a placeholder: {name!r}"
            )
        children = node.get("children", [])
        if not isinstance(children, list):
            raise TreeInvariantError(f"node {name!r} has a non-list 'children' field")
        if is_root:
            root_children = children
        if not children:
            if not isinstance(name, str) or name.strip() == "" or name.strip().lower() in _PLACEHOLDER_STRINGS:
                raise TreeInvariantError(f"tree leaf label is missing or a placeholder: {name!r}")
            leaves.append(name)
            seen_leaves[name] = seen_leaves.get(name, 0) + 1
            return
        for child in children:
            visit(child, depth + 1, is_root=False)

    visit(root, 0, is_root=True)

    duplicates = sorted(label for label, count in seen_leaves.items() if count > 1)
    if duplicates:
        raise TreeInvariantError(f"species tree contains duplicate leaf labels: {duplicates}")
    if len(leaves) < 2:
        raise TreeInvariantError(f"species tree must have at least 2 leaves, got {len(leaves)}")
    assert root_children is not None
    if len(root_children) < 2:
        raise TreeInvariantError(f"species tree root has {len(root_children)} children; at least 2 are required")
    if require_rooted:
        # Rootedness cannot be inferred from topology; demand explicit
        # provenance from the caller (plan section 5.3).
        if rooted is None:
            raise ProvenanceError(
                "require_rooted is set but rootedness was not declared; pass rooted=True or "
                "rooted=False from recorded tree provenance (topology alone cannot establish "
                "biological rooting)"
            )
        if not rooted:
            raise TreeInvariantError("species tree is declared unrooted; a rooted tree is required")
    if require_bifurcating_root and len(root_children) != 2:
        raise TreeInvariantError(
            f"species tree root has {len(root_children)} children; a bifurcating root "
            "requires exactly 2 (note: this checks structure only, not biological rooting)"
        )
