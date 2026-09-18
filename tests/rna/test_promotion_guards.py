"""Promotion-boundary guards for downstream analytical contracts.

Companion tests for docs/rna/downstream_contracts.md. They test the
promotion boundary only — whether a descriptive-stage artifact can acquire
biological-inference labels through the public APIs — and reuse the
statistics-contract invariants (validate_analysis_provenance,
result_role, the gated inferential entry points) rather than re-testing
their internals. Internal validator coverage lives in
test_rna_statistics_contract.py; comparative-layer gating in
test_inferential_comparative.py; the Wilcoxon gate in
test_tissue_specificity.py.

All fixtures are small deterministic real pandas data. No mocks, no
network, no live data root; runs clean under
``-m "not network and not external_tool"``.
"""

from typing import Any

import pandas as pd
import pytest

from metainformant.rna.analysis.cross_species import (
    compute_fingerprint_divergence_matrix,
    compute_fingerprint_stability,
)
from metainformant.rna.analysis.inferential_comparative import (
    ComparativeDesign,
    run_inferential_comparative_analysis,
)
from metainformant.rna.analysis.statistics_contract import (
    DESCRIPTIVE_ROLE,
    INFERENTIAL_ROLE,
    AnalysisProvenance,
    ProvenanceError,
    StatisticsContractError,
    declared_inferential_bh_fdr,
    render_analysis_provenance_block,
    result_role,
)

SOFTWARE_VERSIONS = {"metainformant": "1.0.0", "pandas": "2.2.0", "python": "3.12"}

N_FEATURES = 40
FINGERPRINT_KWARGS: dict[str, Any] = {"n_bins": 8, "min_valid_features": 30}


def _descriptive_contract(**overrides: Any) -> AnalysisProvenance:
    fields: dict[str, Any] = dict(
        analysis_id="downstream_contracts_descriptive_v1",
        estimand="pairwise descriptive dissimilarity between common-binned expression distributions",
        replicate_unit="species finalized matrix (one profile per species)",
        random_seed=20260917,
        resampling_count=50,
        null_model="feature resampling with replacement within each species profile",
        multiple_testing_family=None,
        multiple_testing_method=None,
        tested_feature_count=None,
        software_versions=SOFTWARE_VERSIONS,
        analysis_role=DESCRIPTIVE_ROLE,
    )
    fields.update(overrides)
    return AnalysisProvenance(**fields)


def _inferential_contract(**overrides: Any) -> AnalysisProvenance:
    fields: dict[str, Any] = dict(
        analysis_id="downstream_contracts_inferential_v1",
        estimand="random-effects combined contrast coefficient per feature",
        replicate_unit="biological sample within study",
        random_seed=20260917,
        resampling_count=50,
        null_model="two-sided z-test of the combined contrast coefficient against zero",
        multiple_testing_family="feature",
        multiple_testing_method="bh-fdr",
        tested_feature_count=1,
        software_versions=SOFTWARE_VERSIONS,
        analysis_role=INFERENTIAL_ROLE,
    )
    fields.update(overrides)
    return AnalysisProvenance(**fields)


def _design() -> ComparativeDesign:
    return ComparativeDesign(
        response_col="expression",
        contrast_col="condition",
        reference_level="reference",
        treatment_level="treatment",
        study_col="study",
        feature_col=None,
    )


def _observations() -> pd.DataFrame:
    """Minimal balanced long-format observations for one feature."""
    rows: list[dict[str, Any]] = []
    for study in ("study_a", "study_b", "study_c"):
        for i in range(4):
            treated = i % 2 == 0
            rows.append(
                {
                    "study": study,
                    "condition": "treatment" if treated else "reference",
                    "expression": 1.0 if treated else 0.0,
                }
            )
    return pd.DataFrame(rows)


def _species_profiles() -> dict[str, pd.DataFrame]:
    """Three deterministic species profiles of N_FEATURES features each."""
    profiles: dict[str, pd.DataFrame] = {}
    for species_index, species in enumerate(("Apis_mellifera", "Bombus_terrestris", "Nasonia_vitripennis")):
        values = [float(((feature * 7 + species_index * 3) % 11) + 1) for feature in range(N_FEATURES)]
        profiles[species] = pd.DataFrame({"expression": values}, index=[f"feature_{i}" for i in range(N_FEATURES)])
    return profiles


def _descriptive_artifact() -> tuple[pd.DataFrame, pd.DataFrame]:
    divergence = compute_fingerprint_divergence_matrix(_species_profiles(), **FINGERPRINT_KWARGS)
    stability = compute_fingerprint_stability(
        _species_profiles(),
        divergence,
        n_bootstrap=20,
        random_seed=7,
        **FINGERPRINT_KWARGS,
    )
    return divergence, stability


# =============================================================================
# Descriptive artifacts cannot acquire inference labels
# =============================================================================


def test_refused_promotion_attempts_leave_descriptive_labels_untouched() -> None:
    divergence, stability = _descriptive_artifact()
    descriptive = _descriptive_contract()

    # Every promotion entry point refuses a descriptive-stage artifact.
    with pytest.raises(RuntimeError, match="gated"):
        declared_inferential_bh_fdr(
            [float(value) for value in stability["point_estimate"]], descriptive, evidence_manifest_frozen=False
        )
    with pytest.raises(ValueError, match="missing required columns"):
        run_inferential_comparative_analysis(
            divergence, _design(), _inferential_contract(), evidence_manifest_frozen=True
        )

    # The artifact's role marker is unchanged by the refused attempts: no
    # public API rewrites a descriptive artifact into an inferential one.
    assert divergence.attrs["role"] == DESCRIPTIVE_ROLE
    assert stability.attrs["role"] == DESCRIPTIVE_ROLE
    assert result_role(divergence) == DESCRIPTIVE_ROLE
    assert result_role(stability) == DESCRIPTIVE_ROLE


def test_descriptive_matrix_refuses_as_inferential_observations_even_when_frozen() -> None:
    divergence, stability = _descriptive_artifact()
    for artifact in (divergence, stability):
        with pytest.raises(ValueError, match="missing required columns"):
            run_inferential_comparative_analysis(
                artifact, _design(), _inferential_contract(), evidence_manifest_frozen=True
            )


def test_descriptive_scores_cannot_be_redeclared_as_a_p_value_family() -> None:
    _, stability = _descriptive_artifact()
    scores = [float(value) for value in stability["point_estimate"]]
    with pytest.raises(StatisticsContractError, match="descriptive"):
        declared_inferential_bh_fdr(scores, _descriptive_contract(), evidence_manifest_frozen=True)


# =============================================================================
# Promotion-checklist fields fail closed
# =============================================================================


@pytest.mark.parametrize(
    "overrides",
    [
        {"metadata_harmonization_review": "pending"},
        {"metadata_harmonization_review": ""},
        {"species_tree_source": "todo"},
        {"data_root_snapshot_id": "unknown"},
        {"artifact_paths": {"final_matrix": "tbd"}},
    ],
)
def test_placeholder_promotion_fields_refuse_to_render(overrides: dict[str, Any]) -> None:
    record = _descriptive_contract(**overrides)
    with pytest.raises(ProvenanceError):
        render_analysis_provenance_block(record)


def test_declared_promotion_fields_render_as_checklist_lines() -> None:
    record = _descriptive_contract(
        data_root_snapshot_id="snapshot-2026-09-17",
        cohort_included_count=12,
        cohort_excluded_count=3,
        metadata_harmonization_review="reviewed 2026-09-17",
        species_tree_source="treebase:study-12345",
        species_tree_branch_length_scale="substitutions-per-site",
        artifact_paths={"final_matrix": "output/final/finalized_matrix.tsv"},
    )
    joined = "\n".join(render_analysis_provenance_block(record))
    assert "analysis_provenance_metadata_harmonization_review: reviewed 2026-09-17" in joined
    assert "analysis_provenance_data_root_snapshot_id: snapshot-2026-09-17" in joined
    assert "analysis_provenance_cohort_included_count: 12" in joined
    assert "analysis_provenance_artifact_final_matrix: output/final/finalized_matrix.tsv" in joined


def test_halted_role_cannot_claim_promotion_artifacts() -> None:
    halted = _descriptive_contract(
        analysis_role="stopped",
        artifact_paths={"final_matrix": "output/x.tsv"},
    )
    with pytest.raises(ProvenanceError, match="artifact_paths"):
        render_analysis_provenance_block(halted)


# =============================================================================
# Only the frozen gate produces inferential labels
# =============================================================================


def test_only_the_frozen_gate_produces_inferential_labels() -> None:
    divergence, _ = _descriptive_artifact()
    result = run_inferential_comparative_analysis(
        _observations(), _design(), _inferential_contract(), evidence_manifest_frozen=True
    )
    assert result["role"] == INFERENTIAL_ROLE
    assert result["gate"] == "post-freeze"
    features = result["features"]
    assert result_role(features) == INFERENTIAL_ROLE
    # The promoted output is the only labeled-inferential object in play:
    # the descriptive artifact passed through the boundary stays descriptive,
    # and only the gate output carries p-value columns.
    assert divergence.attrs["role"] == DESCRIPTIVE_ROLE
    assert "p_value" not in divergence.columns
    assert "p_value" in features.columns
