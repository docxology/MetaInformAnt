"""Tests for the cross-species statistics contract (statistics_contract.py).

Covers the fail-closed contract required by
projects/hymenoptera_amalgkit/docs/manuscript/statistical_analysis_plan.md:

- analysis-provenance records: validation, rendering, placeholder refusal;
- descriptive-vs-inferential separation: fingerprint results carry
  ``attrs["role"] == "descriptive"``; the inferential BH-FDR path is gated;
- orthology x species presence invariants (duplicate labels, missing
  mappings, low replication, incomplete coverage);
- species-tree invariants (declared rootedness provenance, optional
  bifurcating-root structure check, duplicate leaves, malformed input).

All fixtures are small deterministic real pandas/numpy data. No mocks,
no network, no live data root.
"""

from __future__ import annotations

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
    OrthologyInvariantError,
    ProvenanceError,
    StatisticsContractError,
    TreeInvariantError,
    benjamini_hochberg_fdr,
    declared_inferential_bh_fdr,
    render_analysis_provenance_block,
    result_role,
    validate_analysis_provenance,
    validate_orthology_profile_invariants,
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
