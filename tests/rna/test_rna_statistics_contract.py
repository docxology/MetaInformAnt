"""Tests for the cross-species statistics contract (statistics_contract.py).

Covers the fail-closed contract required by
projects/hymenoptera_amalgkit/docs/manuscript/statistical_analysis_plan.md:

- analysis-provenance records: validation, rendering, placeholder refusal;
- descriptive-vs-inferential separation: fingerprint results carry
  ``attrs["role"] == "descriptive"``; the inferential BH-FDR path is gated;
- biological-replicate unit and estimand declarations (plan sections 1 and
  3): structured caller-declared units validated for non-degeneracy
  (placeholder labels, technical replicates never counted as independent,
  no stratum below the declared minimum, inferential roles require at
  least 2 replicates per unit), with additive rendering;
- orthology x species presence invariants (duplicate labels, missing
  mappings, low replication, incomplete coverage);
- species-tree invariants (declared rootedness provenance, optional
  bifurcating-root structure check, duplicate leaves, malformed input);
- predeclared MJ-01 comparative designs: covariate/strata declaration
  validation and observation enforcement (unknown covariates, undeclared
  strata, empty and singleton strata all fail closed with typed errors);
- role-conditional paired effect-size records (analytic SE, seeded
  bootstrap CI, gated inferential path, not-applicable rendering for
  descriptive roles);
- Cochran-Q / I-squared heterogeneity records (fail-closed below 2
  studies, descriptive roles exempt) and leave-one-study-out deltas.

All fixtures are small deterministic real pandas/numpy data. No mocks,
no network, no live data root.
"""

import math
from typing import Any

import numpy as np
import pandas as pd
import pytest

from metainformant.rna.analysis.cross_species import (
    compute_fingerprint_divergence_matrix,
    compute_fingerprint_stability,
)
from metainformant.rna.analysis.statistics_contract import (
    DESCRIPTIVE_ROLE,
    INFERENTIAL_ROLE,
    AnalysisProvenance,
    EmptyStratumError,
    EstimandDeclaration,
    HeterogeneityError,
    OrthologyInvariantError,
    PredeclaredDesign,
    ProvenanceError,
    ReplicateUnitDeclaration,
    SensitivityAnalysis,
    SingletonStratumError,
    StatisticsContractError,
    TreeInvariantError,
    UndeclaredStratumError,
    UnknownDesignCovariateError,
    benjamini_hochberg_fdr,
    declared_effect_size,
    declared_inferential_bh_fdr,
    enforce_predeclared_design,
    heterogeneity_record,
    leave_one_study_out_deltas,
    paired_log2fc_effect,
    render_analysis_provenance_block,
    render_effect_size_record,
    result_role,
    validate_analysis_provenance,
    validate_estimand_declaration,
    validate_orthology_profile_invariants,
    validate_predeclared_design,
    validate_replicate_unit_declaration,
    validate_sensitivity_analysis,
    validate_species_tree_invariants,
)

SOFTWARE_VERSIONS = {"metainformant": "1.0.0", "numpy": "2.0.0", "python": "3.12"}


def _provenance(**overrides) -> AnalysisProvenance:
    """A valid descriptive provenance record with optional field overrides."""
    fields = dict(
        analysis_id="hymenoptera_fingerprint_v1",
        estimand="pairwise descriptive dissimilarity between common-binned expression distributions",
        replicate_unit="species finalized matrix (one profile per species)",
        random_seed=20260808,
        resampling_count=200,
        null_model="feature resampling with replacement within each species profile",
        multiple_testing_family=None,
        multiple_testing_method=None,
        tested_feature_count=None,
        software_versions=SOFTWARE_VERSIONS,
        analysis_role=DESCRIPTIVE_ROLE,
    )
    fields.update(overrides)
    return AnalysisProvenance(**fields)


def _sensitivity(**overrides: Any) -> SensitivityAnalysis:
    """A valid predeclared sensitivity analysis with optional field overrides."""
    fields: dict[str, Any] = dict(
        name="resampling_fraction_stability",
        varied_parameter="feature_resampling_fraction",
        baseline_value="0.8",
        varied_values=("0.5", "0.6", "0.9"),
        expected_direction="none",
    )
    fields.update(overrides)
    return SensitivityAnalysis(**fields)


def _replicate_unit(**overrides: Any) -> ReplicateUnitDeclaration:
    """A valid declared replicate unit with optional field overrides."""
    fields: dict[str, Any] = dict(
        name="independent biological replicate nested in study",
        nesting="library nested in study nested in species",
        counts={"apis": 3, "cerana": 2},
        min_independent_replicates=2,
    )
    fields.update(overrides)
    return ReplicateUnitDeclaration(**fields)


def _estimand(**overrides: Any) -> EstimandDeclaration:
    """A valid declared estimand with optional field overrides."""
    fields: dict[str, Any] = dict(
        name="pairwise descriptive dissimilarity between common-binned expression distributions",
        contrast="species pairs under the native fingerprint layer",
        permitted_interpretation="exploratory distribution shape; no gene or evolutionary claim",
    )
    fields.update(overrides)
    return EstimandDeclaration(**fields)


def _species_profiles(n_features: int = 60) -> dict[str, pd.Series]:
    """Deterministic small expression profiles with distinct shapes."""
    rng = np.random.default_rng(42)
    return {
        "sp_a": pd.Series(rng.normal(loc=5.0, scale=1.0, size=n_features), index=[f"f{i}" for i in range(n_features)]),
        "sp_b": pd.Series(rng.exponential(scale=2.0, size=n_features), index=[f"f{i}" for i in range(n_features)]),
        "sp_c": pd.Series(rng.uniform(low=0.5, high=9.0, size=n_features), index=[f"f{i}" for i in range(n_features)]),
    }


# =============================================================================
# Analysis provenance: validation and rendering
# =============================================================================


class TestAnalysisProvenanceValidation:
    def test_valid_record_passes(self) -> None:
        validate_analysis_provenance(_provenance())

    def test_inferential_record_passes(self) -> None:
        validate_analysis_provenance(
            _provenance(
                analysis_role=INFERENTIAL_ROLE,
                multiple_testing_family="unordered species pairs",
                multiple_testing_method="bh-fdr",
                tested_feature_count=3,
            )
        )

    @pytest.mark.parametrize(
        "overrides",
        [
            {"analysis_id": ""},
            {"analysis_id": "  "},
            {"estimand": "TBD"},
            {"estimand": "TODO"},
            {"replicate_unit": "n/a"},
            {"null_model": "placeholder"},
            {"multiple_testing_family": "?"},
            {"multiple_testing_family": "unordered species pairs"},
            {"multiple_testing_method": "eyeball"},
            {"multiple_testing_method": "NONE"},
            {"multiple_testing_method": "bh-fdr"},
            {"analysis_role": INFERENTIAL_ROLE, "multiple_testing_family": None},
            {"analysis_role": INFERENTIAL_ROLE, "multiple_testing_family": "not-applicable"},
            {"analysis_role": INFERENTIAL_ROLE, "multiple_testing_method": None},
            {"analysis_role": INFERENTIAL_ROLE, "multiple_testing_method": "not-applicable"},
            {"analysis_role": INFERENTIAL_ROLE, "multiple_testing_method": "eyeball"},
            {"analysis_role": INFERENTIAL_ROLE, "tested_feature_count": None},
            {"tested_feature_count": 0},
            {"tested_feature_count": 351},
            {"software_versions": {}},
            {"software_versions": {"numpy": "NA"}},
            {"analysis_role": "exploratory"},
        ],
    )
    def test_missing_or_placeholder_fields_fail_closed(self, overrides) -> None:
        with pytest.raises(ProvenanceError):
            validate_analysis_provenance(_provenance(**overrides))

    def test_ad_hoc_dict_is_refused(self) -> None:
        with pytest.raises(TypeError):
            validate_analysis_provenance(_provenance().__dict__)  # type: ignore[arg-type]

    def test_render_block_is_additive_key_value_lines(self) -> None:
        lines = render_analysis_provenance_block(_provenance())
        assert all(": " in line for line in lines)
        assert all(line.startswith("analysis_provenance_") for line in lines)
        joined = "\n".join(lines)
        assert "analysis_provenance_random_seed: 20260808" in joined
        assert "analysis_provenance_multiple_testing_family: not-applicable" in joined
        assert "analysis_provenance_multiple_testing_method: not-applicable" in joined
        assert "analysis_provenance_tested_feature_count: not-applicable" in joined
        assert "analysis_provenance_role: descriptive" in joined

    def test_not_applicable_literal_is_accepted_for_descriptive(self) -> None:
        record = _provenance(multiple_testing_method="not-applicable")
        validate_analysis_provenance(record)
        joined = "\n".join(render_analysis_provenance_block(record))
        assert "analysis_provenance_multiple_testing_method: not-applicable" in joined

    def test_render_block_software_versions_deterministic(self) -> None:
        first = render_analysis_provenance_block(_provenance())
        second = render_analysis_provenance_block(_provenance())
        assert first == second
        versions_line = next(line for line in first if line.startswith("analysis_provenance_software_versions"))
        assert versions_line.endswith("metainformant=1.0.0; numpy=2.0.0; python=3.12")

    def test_render_refuses_placeholder_records(self) -> None:
        with pytest.raises(ProvenanceError):
            render_analysis_provenance_block(_provenance(estimand="unknown"))


# =============================================================================
# Sensitivity-analysis registry (plan section 7)
# =============================================================================


class TestSensitivityAnalysisRegistry:
    def test_constructs_and_validates(self) -> None:
        analysis = _sensitivity()
        assert analysis.name == "resampling_fraction_stability"
        assert analysis.varied_parameter == "feature_resampling_fraction"
        assert analysis.baseline_value == "0.8"
        assert analysis.varied_values == ("0.5", "0.6", "0.9")
        assert analysis.expected_direction == "none"
        assert analysis.notes == ""
        validate_sensitivity_analysis(analysis)

    @pytest.mark.parametrize(
        "overrides",
        [
            {"name": ""},
            {"name": "TBD"},
            {"varied_parameter": "n/a"},
            {"baseline_value": "unknown"},
            {"varied_values": ()},
            {"varied_values": ["0.5", "0.9"]},
            {"varied_values": ("0.5", "")},
            {"varied_values": ("0.5", "tbd")},
            {"expected_direction": "higher"},
            {"expected_direction": ""},
        ],
    )
    def test_placeholder_empty_or_bad_direction_fails_closed(self, overrides: dict[str, Any]) -> None:
        with pytest.raises(ProvenanceError):
            validate_sensitivity_analysis(_sensitivity(**overrides))

    def test_ad_hoc_dict_is_refused(self) -> None:
        with pytest.raises(TypeError):
            validate_sensitivity_analysis(_sensitivity().__dict__)  # type: ignore[arg-type]

    def test_non_analysis_role_must_not_declare_sensitivity_analyses(self) -> None:
        with pytest.raises(ProvenanceError) as excinfo:
            validate_analysis_provenance(
                _provenance(
                    analysis_role="stopped",
                    sensitivity_analyses=(_sensitivity(),),
                )
            )
        assert "sensitivity_analyses" in str(excinfo.value)

    def test_non_tuple_registry_is_refused(self) -> None:
        with pytest.raises(ProvenanceError):
            validate_analysis_provenance(_provenance(sensitivity_analyses=[_sensitivity()]))

    def test_render_covers_index_and_fields(self) -> None:
        record = _provenance(
            sensitivity_analyses=(
                _sensitivity(),
                _sensitivity(name="temperature_shift", varied_values=("15c", "25c")),
            )
        )
        lines = render_analysis_provenance_block(record)
        joined = "\n".join(lines)
        assert "analysis_provenance_sensitivity_1_name: resampling_fraction_stability" in joined
        assert "analysis_provenance_sensitivity_1_varied_parameter: feature_resampling_fraction" in joined
        assert "analysis_provenance_sensitivity_1_baseline_value: 0.8" in joined
        assert "analysis_provenance_sensitivity_1_varied_values: 0.5,0.6,0.9" in joined
        assert "analysis_provenance_sensitivity_1_expected_direction: none" in joined
        # Empty notes are omitted entirely; only non-empty notes render.
        assert "analysis_provenance_sensitivity_1_notes" not in joined
        assert "analysis_provenance_sensitivity_2_name: temperature_shift" in joined
        assert "analysis_provenance_sensitivity_2_varied_values: 15c,25c" in joined
        # Sensitivity lines are additive: they keep the standard key prefix.
        assert all(line.startswith("analysis_provenance_sensitivity_") for line in lines if "_sensitivity_" in line)

    def test_notes_rendered_only_when_non_empty(self) -> None:
        record = _provenance(sensitivity_analyses=(_sensitivity(notes="checked on pilot cohort"),))
        joined = "\n".join(render_analysis_provenance_block(record))
        assert "analysis_provenance_sensitivity_1_notes: checked on pilot cohort" in joined

    def test_full_provenance_round_trip(self) -> None:
        record = _provenance(sensitivity_analyses=(_sensitivity(),))
        validate_analysis_provenance(record)
        lines = render_analysis_provenance_block(record)
        assert lines == render_analysis_provenance_block(record)
        assert all(": " in line for line in lines)
        joined = "\n".join(lines)
        assert "analysis_provenance_sensitivity_1_varied_values: 0.5,0.6,0.9" in joined
        assert "analysis_provenance_sensitivity_1_expected_direction: none" in joined


# =============================================================================
# Biological replicate-unit and estimand declarations (plan sections 1, 3)
# =============================================================================


class TestReplicateUnitDeclaration:
    def test_declared_unit_passes(self) -> None:
        validate_replicate_unit_declaration(_replicate_unit())

    def test_descriptive_single_replicate_per_stratum_passes(self) -> None:
        """Descriptive layers may declare one finalized profile per species."""
        unit = _replicate_unit(
            name="species finalized matrix (one profile per species)",
            counts={"apis": 1, "cerana": 1},
            min_independent_replicates=1,
        )
        validate_replicate_unit_declaration(unit)

    @pytest.mark.parametrize(
        "overrides",
        [
            {"name": ""},
            {"name": "TBD"},
            {"nesting": "n/a"},
            {"technical_replicates_counted": True},
            {"counts": {}},
            {"counts": []},
            {"counts": {"apis": 0, "cerana": 2}},
            {"counts": {"apis": 3, "cerana": True}},
            {"counts": {"apis": 3, "TBD": 2}},
            {"min_independent_replicates": 0},
            {"min_independent_replicates": True},
            {"min_independent_replicates": "2"},
        ],
    )
    def test_placeholder_or_degenerate_declarations_fail_closed(self, overrides: dict[str, Any]) -> None:
        with pytest.raises(ProvenanceError):
            validate_replicate_unit_declaration(_replicate_unit(**overrides))

    def test_stratum_below_declared_minimum_is_degenerate(self) -> None:
        """Low replication: a stratum below the declared minimum is refused, naming the stratum."""
        with pytest.raises(ProvenanceError) as excinfo:
            validate_replicate_unit_declaration(_replicate_unit(counts={"apis": 3, "cerana": 1}))
        message = str(excinfo.value)
        assert "degenerate" in message
        assert "cerana" in message

    def test_inferential_grade_requires_two_replicates(self) -> None:
        """An inferential-role estimand cannot rest on a single replicate per unit."""
        with pytest.raises(ProvenanceError) as excinfo:
            validate_replicate_unit_declaration(
                _replicate_unit(min_independent_replicates=1), require_inferential_grade=True
            )
        assert "at least 2" in str(excinfo.value)

    def test_ad_hoc_dict_is_refused(self) -> None:
        with pytest.raises(TypeError):
            validate_replicate_unit_declaration(_replicate_unit().__dict__)  # type: ignore[arg-type]


class TestEstimandDeclaration:
    def test_declared_estimand_passes(self) -> None:
        validate_estimand_declaration(_estimand())

    def test_embedded_replicate_unit_is_validated(self) -> None:
        embedded = _replicate_unit(counts={"apis": 3, "cerana": 1})
        with pytest.raises(ProvenanceError) as excinfo:
            validate_estimand_declaration(_estimand(replicate_unit=embedded))
        assert "cerana" in str(excinfo.value)

    def test_placeholder_fields_fail_closed(self) -> None:
        for overrides in ({"name": "todo"}, {"contrast": ""}, {"permitted_interpretation": "none"}):
            with pytest.raises(ProvenanceError):
                validate_estimand_declaration(_estimand(**overrides))

    def test_ad_hoc_dict_is_refused(self) -> None:
        with pytest.raises(TypeError):
            validate_estimand_declaration(_estimand().__dict__)  # type: ignore[arg-type]


class TestDeclarationProvenanceIntegration:
    def test_descriptive_record_with_declarations_passes_and_renders(self) -> None:
        unit = _replicate_unit()
        record = _provenance(
            replicate_unit_declaration=unit,
            estimand_declaration=_estimand(replicate_unit=unit),
        )
        validate_analysis_provenance(record)
        joined = "\n".join(render_analysis_provenance_block(record))
        assert "analysis_provenance_replicate_unit_name: independent biological replicate nested in study" in joined
        assert "analysis_provenance_replicate_unit_nesting: library nested in study nested in species" in joined
        assert "analysis_provenance_replicate_unit_min_independent_replicates: 2" in joined
        assert "analysis_provenance_replicate_unit_count_apis: 3" in joined
        assert "analysis_provenance_replicate_unit_count_cerana: 2" in joined
        assert (
            "analysis_provenance_estimand_name: pairwise descriptive dissimilarity between "
            "common-binned expression distributions" in joined
        )
        assert "analysis_provenance_estimand_contrast: species pairs under the native fingerprint layer" in joined
        assert (
            "analysis_provenance_estimand_permitted_interpretation: exploratory distribution shape; "
            "no gene or evolutionary claim" in joined
        )

    def test_render_without_declarations_omits_declaration_lines(self) -> None:
        joined = "\n".join(render_analysis_provenance_block(_provenance()))
        assert "analysis_provenance_replicate_unit_name" not in joined
        assert "analysis_provenance_estimand_name" not in joined

    def test_inferential_role_requires_inferential_grade_unit(self) -> None:
        """Role-conditional non-degeneracy: inferential contracts need >=2 replicates."""
        unit = _replicate_unit(min_independent_replicates=1)
        with pytest.raises(ProvenanceError) as excinfo:
            validate_analysis_provenance(
                _provenance(
                    analysis_role=INFERENTIAL_ROLE,
                    multiple_testing_family="ortholog features within species contrasts",
                    multiple_testing_method="bh-fdr",
                    tested_feature_count=3,
                    replicate_unit_declaration=unit,
                )
            )
        assert "at least 2" in str(excinfo.value)

    def test_inferential_role_accepts_grade_two_unit(self) -> None:
        validate_analysis_provenance(
            _provenance(
                analysis_role=INFERENTIAL_ROLE,
                multiple_testing_family="ortholog features within species contrasts",
                multiple_testing_method="bh-fdr",
                tested_feature_count=3,
                replicate_unit_declaration=_replicate_unit(),
            )
        )

    def test_conflicting_embedded_unit_fails_closed(self) -> None:
        record = _provenance(
            replicate_unit_declaration=_replicate_unit(),
            estimand_declaration=_estimand(replicate_unit=_replicate_unit(counts={"apis": 4, "cerana": 2})),
        )
        with pytest.raises(ProvenanceError, match="conflicting replicate-unit declarations"):
            validate_analysis_provenance(record)

    def test_non_analysis_role_must_not_declare_unit_or_estimand(self) -> None:
        with pytest.raises(ProvenanceError) as excinfo:
            validate_analysis_provenance(
                _provenance(analysis_role="stopped", replicate_unit_declaration=_replicate_unit())
            )
        assert "replicate_unit_declaration" in str(excinfo.value)
        with pytest.raises(ProvenanceError) as excinfo:
            validate_analysis_provenance(_provenance(analysis_role="unavailable", estimand_declaration=_estimand()))
        assert "estimand_declaration" in str(excinfo.value)


# =============================================================================
# Descriptive vs inferential separation
# =============================================================================


class TestDescriptiveRoleMarkers:
    def test_fingerprint_divergence_result_is_marked_descriptive(self) -> None:
        profiles = _species_profiles()
        divergence = compute_fingerprint_divergence_matrix(profiles, n_bins=8, min_valid_features=30)
        assert result_role(divergence) == DESCRIPTIVE_ROLE

    def test_fingerprint_stability_result_is_marked_descriptive(self) -> None:
        profiles = _species_profiles()
        divergence = compute_fingerprint_divergence_matrix(profiles, n_bins=8, min_valid_features=30)
        stability = compute_fingerprint_stability(
            profiles,
            divergence,
            n_bins=8,
            min_valid_features=30,
            n_bootstrap=20,
            random_seed=7,
        )
        assert result_role(stability) == DESCRIPTIVE_ROLE
        assert "p_value" not in stability.columns
        assert "confidence" not in " ".join(stability.columns)

    def test_unlabeled_results_are_refused(self) -> None:
        with pytest.raises(StatisticsContractError):
            result_role(pd.DataFrame({"a": [1]}))


class TestBenjaminiHochberg:
    def test_known_adjustment_values(self) -> None:
        assert benjamini_hochberg_fdr([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.04, 0.04])

    def test_adjusted_values_never_below_raw_and_bounded(self) -> None:
        raw = [0.5, 0.1, 0.9, 0.02, 0.7, 0.3]
        adjusted = benjamini_hochberg_fdr(raw)
        for r, q in zip(raw, adjusted):
            assert q >= r
            assert 0.0 <= q <= 1.0

    @pytest.mark.parametrize("bad", [[], [float("nan")], [1.5], [-0.1], [0.3, "x"], [None]])
    def test_invalid_inputs_fail_closed(self, bad) -> None:
        with pytest.raises(StatisticsContractError):
            benjamini_hochberg_fdr(bad)  # type: ignore[arg-type]


class TestGatedInferentialPath:
    def _inferential_contract(self, tested_feature_count: int = 3) -> AnalysisProvenance:
        return _provenance(
            analysis_role=INFERENTIAL_ROLE,
            analysis_id="hymenoptera_confirmatory_v1",
            estimand="effect of caste state on ortholog expression within species",
            replicate_unit="independent biological replicate nested in study",
            resampling_count=1000,
            null_model="permutation of caste labels within species",
            multiple_testing_family="ortholog features within species contrasts",
            multiple_testing_method="bh-fdr",
            tested_feature_count=tested_feature_count,
        )

    def test_gate_refuses_by_default(self) -> None:
        with pytest.raises(RuntimeError, match="gated"):
            declared_inferential_bh_fdr([0.01, 0.02, 0.5], self._inferential_contract())

    def test_gate_refuses_descriptive_contracts(self) -> None:
        with pytest.raises(StatisticsContractError):
            declared_inferential_bh_fdr(
                [0.01, 0.02, 0.5],
                _provenance(),  # descriptive role
                evidence_manifest_frozen=True,
            )

    def test_family_size_must_match_contract(self) -> None:
        with pytest.raises(StatisticsContractError, match="tested_feature_count"):
            declared_inferential_bh_fdr(
                [0.01, 0.02],
                self._inferential_contract(tested_feature_count=3),
                evidence_manifest_frozen=True,
            )

    def test_frozen_inferential_contract_applies_declared_bh(self) -> None:
        result = declared_inferential_bh_fdr(
            [0.01, 0.04, 0.03],
            self._inferential_contract(),
            evidence_manifest_frozen=True,
        )
        assert result["role"] == INFERENTIAL_ROLE
        assert result["multiple_testing_method"] == "bh-fdr"
        assert result["gate"] == "post-freeze"
        assert result["adjusted_p_values"] == pytest.approx([0.03, 0.04, 0.04])
        assert result["random_seed"] == 20260808

    def test_gated_refusal_for_placeholder_contract(self) -> None:
        broken = AnalysisProvenance(
            analysis_id="x",
            estimand="TBD",
            replicate_unit="replicate",
            random_seed=1,
            resampling_count=10,
            null_model="none",
            multiple_testing_family="family",
            multiple_testing_method="bh-fdr",
            tested_feature_count=2,
            software_versions={"metainformant": "1.0.0"},
            analysis_role=INFERENTIAL_ROLE,
        )
        with pytest.raises(ProvenanceError):
            declared_inferential_bh_fdr([0.1, 0.2], broken, evidence_manifest_frozen=True)


# =============================================================================
# Orthology x species presence invariants
# =============================================================================


class TestOrthologyProfileInvariants:
    def _presence(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "apis": [1, 1, 0],
                "cerana": [1, 0, 1],
                "mellifera": [0, 1, 1],
            },
            index=["OG1", "OG2", "OG3"],
        )

    def test_valid_table_passes(self) -> None:
        validate_orthology_profile_invariants(self._presence())

    def test_strict_fraction_passes_for_full_coverage(self) -> None:
        full = self._presence() * 1  # OG1/OG2/OG3 each cover 2 of 3 species
        with pytest.raises(OrthologyInvariantError):
            validate_orthology_profile_invariants(full, min_species_fraction=1.0)
        fully_covered = pd.DataFrame(
            {"a": [1, 1], "b": [1, 1], "c": [1, 1], "d": [1, 1]},
            index=["OG1", "OG2"],
        )
        validate_orthology_profile_invariants(fully_covered, min_species_fraction=1.0)

    def test_duplicate_orthogroup_labels_fail_closed(self) -> None:
        table = pd.DataFrame(
            {"a": [1, 1, 0], "b": [1, 0, 1], "c": [0, 1, 1]},
            index=["OG1", "OG1", "OG2"],
        )
        with pytest.raises(OrthologyInvariantError, match="duplicate orthogroup labels"):
            validate_orthology_profile_invariants(table)

    def test_duplicate_species_labels_fail_closed(self) -> None:
        table = pd.DataFrame(np.ones((2, 3)), index=["OG1", "OG2"], columns=["apis", "apis", "cerana"])
        with pytest.raises(OrthologyInvariantError, match="duplicate species labels"):
            validate_orthology_profile_invariants(table)

    def test_missing_orthology_mapping_fails_closed(self) -> None:
        table = self._presence().astype(float)
        table.iloc[0, 0] = np.nan
        with pytest.raises(OrthologyInvariantError, match="missing values"):
            validate_orthology_profile_invariants(table)

    def test_ambiguous_presence_values_fail_closed(self) -> None:
        table = self._presence().astype(float)
        table.iloc[1, 2] = 0.5
        with pytest.raises(OrthologyInvariantError, match="outside"):
            validate_orthology_profile_invariants(table)

    def test_string_presence_columns_fail_closed(self) -> None:
        table = pd.DataFrame({"a": ["yes", "yes"], "b": ["yes", "no"], "c": ["no", "yes"]}, index=["OG1", "OG2"])
        with pytest.raises(OrthologyInvariantError, match="non-numeric"):
            validate_orthology_profile_invariants(table)

    def test_placeholder_labels_fail_closed(self) -> None:
        table = self._presence()
        table.index = ["OG1", "NA", "OG3"]
        with pytest.raises(OrthologyInvariantError, match="placeholder labels"):
            validate_orthology_profile_invariants(table)

    def test_low_replication_fails_closed(self) -> None:
        table = pd.DataFrame(
            {"a": [1, 1, 0], "b": [0, 1, 0], "c": [0, 1, 0], "d": [1, 1, 0]},
            index=["OG_ok", "OG_ok2", "OG_lonely"],
        )
        with pytest.raises(OrthologyInvariantError, match="OG_lonely"):
            validate_orthology_profile_invariants(table)

    def test_incomplete_coverage_below_declared_threshold_fails_closed(self) -> None:
        table = self._presence()  # every orthogroup covers exactly 2 of 3 species
        with pytest.raises(OrthologyInvariantError, match="coverage threshold"):
            validate_orthology_profile_invariants(table, min_species_fraction=1.0)

    def test_invalid_thresholds_fail_closed(self) -> None:
        with pytest.raises(OrthologyInvariantError):
            validate_orthology_profile_invariants(self._presence(), min_species_fraction=1.5)
        with pytest.raises(OrthologyInvariantError):
            validate_orthology_profile_invariants(self._presence(), min_species_per_orthogroup=0)

    def test_empty_and_single_species_tables_fail_closed(self) -> None:
        with pytest.raises(OrthologyInvariantError):
            validate_orthology_profile_invariants(pd.DataFrame())
        single = pd.DataFrame({"a": [1, 1]}, index=["OG1", "OG2"])
        with pytest.raises(OrthologyInvariantError, match="at least 2 species"):
            validate_orthology_profile_invariants(single)


# =============================================================================
# Species tree invariants
# =============================================================================


class TestSpeciesTreeInvariants:
    def test_rooted_newick_passes(self) -> None:
        validate_species_tree_invariants("((apis,cerana),(mellifera,bombus));", rooted=True)

    def test_unnamed_newick_root_passes(self) -> None:
        validate_species_tree_invariants("(apis,cerana);", rooted=True)

    def test_newick_with_branch_lengths_passes(self) -> None:
        validate_species_tree_invariants("((apis:0.10,cerana:0.12)Apinae:0.20,mellifera_lineage:0.30);", rooted=True)

    def test_rooted_dict_tree_passes(self) -> None:
        tree = {
            "name": "root",
            "children": [
                {"name": "apis", "distance": 0.1},
                {"name": "cerana", "distance": 0.1},
            ],
        }
        validate_species_tree_invariants(tree, rooted=True)

    def test_rooted_polytomy_passes_when_declared_rooted(self) -> None:
        """A polytomous root is a valid hard polytomy, not proof of unrootedness."""
        validate_species_tree_invariants("(apis,cerana,mellifera);", rooted=True)

    def test_undeclared_rootedness_fails_closed(self) -> None:
        """Rootedness is caller-declared provenance; topology cannot establish it."""
        with pytest.raises(ProvenanceError, match="rootedness was not declared"):
            validate_species_tree_invariants("((apis,cerana),(mellifera,bombus));")
        with pytest.raises(ProvenanceError, match="rootedness was not declared"):
            validate_species_tree_invariants({"name": "root", "children": [{"name": "apis"}, {"name": "cerana"}]})

    def test_declared_unrooted_newick_fails_closed(self) -> None:
        with pytest.raises(TreeInvariantError, match="declared unrooted"):
            validate_species_tree_invariants("(apis,cerana,mellifera);", rooted=False)

    def test_declared_unrooted_dict_tree_fails_closed(self) -> None:
        tree = {
            "name": "root",
            "children": [
                {"name": "apis"},
                {"name": "cerana"},
                {"name": "mellifera"},
            ],
        }
        with pytest.raises(TreeInvariantError, match="declared unrooted"):
            validate_species_tree_invariants(tree, rooted=False)

    def test_rootedness_not_required_skips_declaration(self) -> None:
        """With require_rooted=False, no rootedness declaration is demanded."""
        validate_species_tree_invariants("(apis,cerana,mellifera);", require_rooted=False)

    def test_bifurcating_root_check_only_fires_when_requested(self) -> None:
        polytomy = "(apis,cerana,mellifera);"
        validate_species_tree_invariants(polytomy, rooted=True)
        with pytest.raises(TreeInvariantError, match="bifurcating"):
            validate_species_tree_invariants(polytomy, rooted=True, require_bifurcating_root=True)
        validate_species_tree_invariants(
            "((apis,cerana),(mellifera,bombus));", rooted=True, require_bifurcating_root=True
        )

    def test_duplicate_leaf_labels_fail_closed(self) -> None:
        with pytest.raises(TreeInvariantError, match="duplicate leaf labels"):
            validate_species_tree_invariants("((apis,cerana),(apis,mellifera));")

    def test_duplicate_leaf_labels_in_dict_fail_closed(self) -> None:
        tree = {
            "name": "root",
            "children": [{"name": "apis"}, {"name": "apis"}],
        }
        with pytest.raises(TreeInvariantError, match="duplicate leaf labels"):
            validate_species_tree_invariants(tree)

    def test_malformed_newick_fails_closed(self) -> None:
        with pytest.raises(TreeInvariantError):
            validate_species_tree_invariants("(apis,cerana")
        with pytest.raises(TreeInvariantError):
            validate_species_tree_invariants("(apis,cerana)")
        with pytest.raises(TreeInvariantError):
            validate_species_tree_invariants("(apis:oops,cerana);")

    def test_dict_tree_missing_name_fails_closed(self) -> None:
        with pytest.raises(TreeInvariantError):
            validate_species_tree_invariants({"children": [{"name": "apis"}, {"name": "cerana"}]})

    def test_placeholder_labels_fail_closed(self) -> None:
        with pytest.raises(TreeInvariantError, match="placeholder"):
            validate_species_tree_invariants({"name": "root", "children": [{"name": "NA"}, {"name": "apis"}]})
        with pytest.raises(TreeInvariantError, match="placeholder"):
            validate_species_tree_invariants("((apis,cerana),TBD);")

    def test_single_leaf_tree_fails_closed(self) -> None:
        with pytest.raises(TreeInvariantError, match="at least 2"):
            validate_species_tree_invariants("(apis);")

    def test_invalid_input_type_fails_closed(self) -> None:
        with pytest.raises(TreeInvariantError):
            validate_species_tree_invariants(["apis", "cerana"])  # type: ignore[arg-type]


# =============================================================================
# Predeclared comparative design (MJ-01)
# =============================================================================


def _design(**overrides: Any) -> PredeclaredDesign:
    """A valid predeclared comparative design with optional field overrides."""

    fields: dict[str, Any] = dict(
        covariate_strata={
            "study": ("study_a", "study_b"),
            "sex": ("female", "male"),
        }
    )
    fields.update(overrides)
    return PredeclaredDesign(**fields)


def _lane_observations() -> pd.DataFrame:
    """Deterministic observations matching ``_design()`` (no thin strata)."""

    rows = []
    for study in ("study_a", "study_b"):
        for sex in ("female", "male"):
            for _ in range(3):
                rows.append({"study": study, "sex": sex, "response": float(len(rows))})
    return pd.DataFrame(rows)


class TestPredeclaredDesignValidation:
    def test_declared_design_passes(self) -> None:
        validate_predeclared_design(_design())

    @pytest.mark.parametrize(
        "overrides",
        [
            {"covariate_strata": {}},
            {"covariate_strata": []},
            {"covariate_strata": {"study": ()}},
            {"covariate_strata": {"study": ("a", "a")}},
            {"covariate_strata": {"study": ("a", "")}},
            {"covariate_strata": {"study": ("a", "TBD")}},
            {"covariate_strata": {"TBD": ("a",)}},
            {"covariate_strata": {"study": ["a", "b"]}},
            {"covariate_strata": {"study": "a,b"}},
        ],
    )
    def test_placeholder_or_degenerate_declarations_fail_closed(self, overrides: dict[str, Any]) -> None:
        with pytest.raises(ProvenanceError):
            validate_predeclared_design(_design(**overrides))

    def test_ad_hoc_dict_is_refused(self) -> None:
        with pytest.raises(TypeError):
            validate_predeclared_design(_design().__dict__)  # type: ignore[arg-type]


class TestDesignEnforcement:
    def test_matching_design_returns_deterministic_counts(self) -> None:
        counts = enforce_predeclared_design(
            _design(),
            _lane_observations(),
            study_col="study",
            covariate_cols=("sex",),
        )
        assert counts == {
            "sex": {"female": 6, "male": 6},
            "study": {"study_a": 6, "study_b": 6},
        }

    def test_contrast_column_must_be_declared(self) -> None:
        """A contrast column in the request but not the declaration fails closed."""

        frame = _lane_observations().assign(caste=["worker"] * 12)
        with pytest.raises(UnknownDesignCovariateError, match="caste"):
            enforce_predeclared_design(
                _design(),
                frame,
                study_col="study",
                contrast_col="caste",
                covariate_cols=("sex",),
            )

    def test_declared_column_absent_from_observations_refused(self) -> None:
        design = _design(
            covariate_strata={
                "study": ("study_a", "study_b"),
                "sex": ("female", "male"),
                "tissue": ("brain", "muscle"),
            }
        )
        with pytest.raises(UnknownDesignCovariateError, match="tissue"):
            enforce_predeclared_design(
                design,
                _lane_observations(),
                study_col="study",
                covariate_cols=("sex", "tissue"),
            )

    def test_undeclared_stratum_refused(self) -> None:
        frame = _lane_observations()
        frame.loc[0, "sex"] = "intersex"
        with pytest.raises(UndeclaredStratumError, match="intersex"):
            enforce_predeclared_design(
                _design(),
                frame,
                study_col="study",
                covariate_cols=("sex",),
            )

    def test_empty_stratum_refused(self) -> None:
        design = _design(
            covariate_strata={
                "study": ("study_a", "study_b"),
                "sex": ("female", "male", "intersex"),
            }
        )
        with pytest.raises(EmptyStratumError, match="intersex"):
            enforce_predeclared_design(
                design,
                _lane_observations(),
                study_col="study",
                covariate_cols=("sex",),
            )

    def test_singleton_stratum_refused(self) -> None:
        frame = _lane_observations()
        frame.loc[0, "sex"] = "intersex"
        design = _design(
            covariate_strata={
                "study": ("study_a", "study_b"),
                "sex": ("female", "male", "intersex"),
            }
        )
        with pytest.raises(SingletonStratumError, match="intersex"):
            enforce_predeclared_design(
                design,
                frame,
                study_col="study",
                covariate_cols=("sex",),
            )

    def test_declared_minimum_override_tightens_singleton_refusal(self) -> None:
        """With a declared minimum of 4, a two-observation stratum is refused."""

        frame = _lane_observations().drop(index=[7, 8, 9, 10])  # study_b keeps 2 rows
        with pytest.raises(SingletonStratumError, match="study_b"):
            enforce_predeclared_design(
                _design(
                    covariate_strata={
                        "study": ("study_a", "study_b"),
                        "sex": ("female", "male"),
                    }
                ),
                frame,
                study_col="study",
                covariate_cols=("sex",),
                min_stratum_observations=4,
            )

    def test_invalid_minimum_refused(self) -> None:
        with pytest.raises(ProvenanceError):
            enforce_predeclared_design(
                _design(),
                _lane_observations(),
                study_col="study",
                covariate_cols=("sex",),
                min_stratum_observations=0,
            )

    def test_non_dataframe_observations_refused(self) -> None:
        with pytest.raises(TypeError):
            enforce_predeclared_design(
                _design(),
                [{"study": "study_a", "sex": "female"}],
                study_col="study",
                covariate_cols=("sex",),
            )


# =============================================================================
# Role-conditional paired effect sizes
# =============================================================================


def _paired_values() -> tuple[list[float], list[float]]:
    """Paired observations whose log2 ratios are exactly [1, 2, 2, 1]."""

    return ([1.0, 1.0, 2.0, 4.0], [2.0, 4.0, 8.0, 8.0])


def _inferential_provenance(**overrides: Any) -> AnalysisProvenance:
    """A validated inferential provenance record with optional overrides."""

    fields = dict(
        analysis_id="hymenoptera_effect_size_v1",
        estimand="queen-versus-worker paired log2 fold-change per ortholog",
        replicate_unit="independent biological replicate nested in study",
        random_seed=20260922,
        resampling_count=200,
        null_model="caste labels exchangeable within study",
        multiple_testing_family="ortholog features within species contrasts",
        multiple_testing_method="bh-fdr",
        tested_feature_count=4,
        software_versions={"metainformant": "1.0.0"},
        analysis_role=INFERENTIAL_ROLE,
    )
    fields.update(overrides)
    return AnalysisProvenance(**fields)


class TestPairedLog2FcEffect:
    def test_point_effect_and_analytic_se(self) -> None:
        reference, treatment = _paired_values()
        record = paired_log2fc_effect(reference, treatment)
        assert record["effect"] == pytest.approx(1.5)
        assert record["se"] == pytest.approx(1.0 / math.sqrt(12.0))
        assert record["n_pairs"] == 4
        assert "ci_low" not in record and "ci_high" not in record

    def test_bootstrap_ci_is_deterministic_and_brackets_effect(self) -> None:
        reference, treatment = _paired_values()
        first = paired_log2fc_effect(reference, treatment, bootstrap_replicates=100, random_seed=7)
        second = paired_log2fc_effect(reference, treatment, bootstrap_replicates=100, random_seed=7)
        assert first == second
        assert first["ci_low"] <= first["effect"] <= first["ci_high"]
        assert first["bootstrap_n_success"] == 100

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"bootstrap_replicates": -1},
            {"random_seed": -1},
            {"bootstrap_replicates": True},
        ],
    )
    def test_invalid_resampling_options_fail_closed(self, kwargs: dict[str, Any]) -> None:
        reference, treatment = _paired_values()
        with pytest.raises(StatisticsContractError):
            paired_log2fc_effect(reference, treatment, **kwargs)

    def test_invalid_pairs_fail_closed(self) -> None:
        reference, treatment = _paired_values()
        with pytest.raises(StatisticsContractError):
            paired_log2fc_effect(reference, treatment[:3])
        with pytest.raises(StatisticsContractError):
            paired_log2fc_effect([], [])
        with pytest.raises(StatisticsContractError, match="positive"):
            paired_log2fc_effect([0.0, 1.0], [2.0, 3.0])
        with pytest.raises(StatisticsContractError):
            paired_log2fc_effect([float("nan"), 1.0], [2.0, 3.0])


class TestDeclaredEffectSize:
    def test_descriptive_contract_renders_not_applicable(self) -> None:
        reference, treatment = _paired_values()
        record = declared_effect_size(reference, treatment, _provenance())
        assert record["role"] == DESCRIPTIVE_ROLE
        assert record["effect"] == pytest.approx(1.5)
        assert record["replicate_count"] == 4
        for field in (
            "effect_size_method",
            "ci_low",
            "ci_high",
            "bootstrap_resampling_count",
            "multiple_testing_family",
            "multiple_testing_method",
            "gate",
        ):
            assert record[field] == "not-applicable", field

    def test_descriptive_record_renders_not_applicable_lines(self) -> None:
        reference, treatment = _paired_values()
        lines = render_effect_size_record(declared_effect_size(reference, treatment, _provenance()))
        joined = "\n".join(lines)
        assert all(line.startswith("effect_size_") for line in lines)
        assert "effect_size_multiple_testing_family: not-applicable" in joined
        assert "effect_size_effect: 1.5" in joined
        assert lines == render_effect_size_record(declared_effect_size(reference, treatment, _provenance()))

    def test_inferential_gate_refuses_unfrozen_manifest(self) -> None:
        reference, treatment = _paired_values()
        with pytest.raises(RuntimeError, match="gated"):
            declared_effect_size(reference, treatment, _inferential_provenance())

    def test_inferential_single_pair_refused(self) -> None:
        with pytest.raises(StatisticsContractError, match="at least 2 paired"):
            declared_effect_size(
                [1.0],
                [2.0],
                _inferential_provenance(),
                evidence_manifest_frozen=True,
            )

    def test_inferential_frozen_record_carries_provenance(self) -> None:
        reference, treatment = _paired_values()
        record = declared_effect_size(
            reference,
            treatment,
            _inferential_provenance(),
            evidence_manifest_frozen=True,
        )
        assert record["role"] == INFERENTIAL_ROLE
        assert record["gate"] == "post-freeze"
        assert record["effect_size_method"] == "paired-log2fc-mean-bootstrap-percentile"
        assert record["multiple_testing_family"] == "ortholog features within species contrasts"
        assert record["multiple_testing_method"] == "bh-fdr"
        assert record["replicate_count"] == 4
        assert record["bootstrap_resampling_count"] == 200
        assert record["random_seed"] == 20260922
        assert record["ci_low"] <= record["effect"] <= record["ci_high"]
        # Deterministic under the contract's own seed.
        again = declared_effect_size(
            reference,
            treatment,
            _inferential_provenance(),
            evidence_manifest_frozen=True,
        )
        assert (record["ci_low"], record["ci_high"]) == (again["ci_low"], again["ci_high"])

    def test_inferential_record_carries_replicate_unit_counts(self) -> None:
        reference, treatment = _paired_values()
        record = declared_effect_size(
            reference,
            treatment,
            _inferential_provenance(replicate_unit_declaration=_replicate_unit()),
            evidence_manifest_frozen=True,
        )
        assert record["replicate_unit_counts"] == {"apis": 3, "cerana": 2}

    def test_stopped_contract_produces_no_record(self) -> None:
        reference, treatment = _paired_values()
        with pytest.raises(StatisticsContractError, match="halted or unavailable"):
            declared_effect_size(
                reference,
                treatment,
                _provenance(analysis_role="stopped"),
                evidence_manifest_frozen=True,
            )

    def test_render_refuses_unlabeled_or_unknown_roles(self) -> None:
        with pytest.raises(StatisticsContractError):
            render_effect_size_record({"effect": 1.0})
        with pytest.raises(StatisticsContractError, match="unknown role"):
            render_effect_size_record({"role": "exploratory", "effect": 1.0})


# =============================================================================
# Multi-study heterogeneity and leave-one-study-out sensitivity
# =============================================================================


def _study_effects() -> tuple[list[float], list[float]]:
    """Three aligned study effects with standard errors."""

    return ([1.0, 1.0, 1.0], [0.1, 0.2, 0.3])


class TestHeterogeneityRecord:
    def test_inferential_homogeneous_studies_give_zero_heterogeneity(self) -> None:
        effects, ses = _study_effects()
        record = heterogeneity_record(effects, ses, _inferential_provenance())
        assert record["role"] == INFERENTIAL_ROLE
        assert record["heterogeneity_status"] == "computed"
        assert record["cochran_q"] == pytest.approx(0.0)
        assert record["df"] == 2
        assert record["i_squared_percent"] == pytest.approx(0.0)
        assert record["tau_squared"] == pytest.approx(0.0)
        assert record["n_studies"] == 3

    def test_inferential_heterogeneous_studies_give_positive_i_squared(self) -> None:
        effects, ses = [1.0, 2.0], [0.1, 0.1]
        record = heterogeneity_record(effects, ses, _inferential_provenance())
        assert record["cochran_q"] == pytest.approx(50.0)
        assert record["df"] == 1
        assert record["i_squared_percent"] == pytest.approx(98.0)

    def test_inferential_single_study_refused(self) -> None:
        with pytest.raises(HeterogeneityError, match="at least 2 studies"):
            heterogeneity_record([1.0], [0.1], _inferential_provenance())

    def test_descriptive_role_is_exempt(self) -> None:
        effects, ses = _study_effects()
        record = heterogeneity_record(effects, ses, _provenance())
        assert record["role"] == DESCRIPTIVE_ROLE
        assert record["heterogeneity_status"] == "exempt"
        assert record["cochran_q"] == "not-applicable"
        assert record["n_studies"] == 3

    def test_bad_standard_errors_refused(self) -> None:
        effects, ses = _study_effects()
        for bad in ([0.0, 0.2, 0.3], [-0.1, 0.2, 0.3], [float("nan"), 0.2, 0.3]):
            with pytest.raises(HeterogeneityError):
                heterogeneity_record(effects, bad, _inferential_provenance())

    def test_stopped_contract_refused(self) -> None:
        effects, ses = _study_effects()
        with pytest.raises(StatisticsContractError, match="halted or unavailable"):
            heterogeneity_record(effects, ses, _provenance(analysis_role="stopped"))


class TestLeaveOneStudyOut:
    def _studies(self) -> dict[str, tuple[float, float]]:
        return {"a": (1.0, 0.1), "b": (2.0, 0.1), "c": (3.0, 0.1)}

    def test_full_effect_and_per_exclusion_deltas(self) -> None:
        result = leave_one_study_out_deltas(self._studies())
        assert result["full_effect"] == pytest.approx(2.0)
        assert result["deltas"] == {
            "a": pytest.approx(0.5),
            "b": pytest.approx(0.0),
            "c": pytest.approx(-0.5),
        }

    def test_deterministic(self) -> None:
        assert leave_one_study_out_deltas(self._studies()) == leave_one_study_out_deltas(self._studies())

    def test_fewer_than_three_studies_refused(self) -> None:
        with pytest.raises(HeterogeneityError, match="at least 3 studies"):
            leave_one_study_out_deltas({"a": (1.0, 0.1), "b": (2.0, 0.1)})
        with pytest.raises(HeterogeneityError):
            leave_one_study_out_deltas({})

    def test_invalid_inputs_refused(self) -> None:
        with pytest.raises(HeterogeneityError):
            leave_one_study_out_deltas({"a": (1.0, 0.0), "b": (2.0, 0.1), "c": (3.0, 0.1)})
        with pytest.raises(HeterogeneityError):
            leave_one_study_out_deltas([1.0, 2.0, 3.0])  # type: ignore[arg-type]


# =============================================================================
# Design declaration provenance integration
# =============================================================================


class TestDesignDeclarationProvenanceIntegration:
    def test_declared_design_renders_additively(self) -> None:
        record = _provenance(design_declaration=_design())
        validate_analysis_provenance(record)
        joined = "\n".join(render_analysis_provenance_block(record))
        assert "analysis_provenance_design_covariate_sex: female,male" in joined
        assert "analysis_provenance_design_covariate_study: study_a,study_b" in joined

    def test_render_without_design_omits_design_lines(self) -> None:
        joined = "\n".join(render_analysis_provenance_block(_provenance()))
        assert "analysis_provenance_design_covariate_" not in joined

    def test_non_analysis_role_must_not_declare_design(self) -> None:
        with pytest.raises(ProvenanceError) as excinfo:
            validate_analysis_provenance(_provenance(analysis_role="stopped", design_declaration=_design()))
        assert "design_declaration" in str(excinfo.value)
