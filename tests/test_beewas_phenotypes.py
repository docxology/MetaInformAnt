"""Tests for BeeWAS phenotype workbook curation."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


PIPELINE_DIR = Path(__file__).resolve().parents[1] / "scripts" / "gwas" / "pipelines"
if str(PIPELINE_DIR) not in sys.path:
    sys.path.insert(0, str(PIPELINE_DIR))

from beewas_phenotypes import (  # noqa: E402
    bh_fdr,
    bootstrap_mean_ci,
    build_phenotype_knn,
    build_model_results,
    build_phenotype_statistics,
    build_trait_readiness,
    curate_beewas_phenotypes,
    q_value_label,
    nested_variance_partition,
    normalize_colony,
    permutation_anova_p,
    sample_parts,
    standardized_mean_difference,
)


pytest.importorskip("openpyxl")


def write_fixture_workbook(path: Path) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        pd.DataFrame(
            {
                "Line": ["C1", "M10"],
                "Grafted": [13, 22],
                "In_vivo": [12, 27],
                "In_vitro_W": [8, 10],
                "In_vitro_Q": [10, 10],
                "Worker": [100, 100],
            }
        ).to_excel(writer, sheet_name="0. Experimental design", index=False)

        modified = pd.DataFrame(
            [
                ["This is per-nestmate phenotypic measurements (Queen score, etc) and their DNA extraction plate ID.", None, None, None, None, None, None, None, None, None],
                [None] * 10,
                ["Strain", "colony-by-strain", "Sample ID (plate)", "SampleID_LessIV", "Phenotype", "DaysToEmerge", "QueenScore", "BEE#", "Sample ID (row)", "Sample ID (column)"],
                ["C", "C1", "IV1", 1, "Q", 14, 3, 1, "H", 12],
                ["C", "C1", "IV2", 2, "W", 16, 0, 1, "H", 5],
                ["M", "m10", "IV3", 3, "Q", 12, 2, 1, "A", 1],
            ]
        )
        modified.to_excel(writer, sheet_name="1. Modified_Selected_samples_fo", index=False, header=False)

        pd.DataFrame(
            {
                "colony origin": ["C1", "C1", "m10"],
                "Sample ID (plate)": ["IV1", "IV2", "IV3"],
                "Sample ID (row)": ["H", "H", "A"],
                "Sample ID (column)": [12, 5, 1],
                "Phenotype": ["Q", "W", "Q"],
                "Days to emerge": [14, 16, 12],
                "Score": [3, 0, 2],
                "BEE#": [1, 1, 1],
                "blank": [None, None, None],
                "Line": ["C1", "C1", "M10"],
                "Phenotype.1": ["Q", "W", "Q"],
                "Count": [1, 1, 1],
            }
        ).to_excel(writer, sheet_name="Nestmate_measurements_322", index=False)

        pd.DataFrame(
            {
                "GenomicsID": ["C1G", "C1ITQ", "C1WORK", "M10ITW", "M10IV"],
                "Population": ["C", "C", "C", "M", "M"],
                "Colony": ["C1G", "C1", "C1", "M10", "M10"],
                "Type": [None, None, None, None, None],
                "Phenotype": [None, None, None, None, None],
                "Days to emerge": ["There are 5 biological groups", "ITQ", None, None, None],
                "Score": [None, "In vitro queen", None, None, None],
            }
        ).to_excel(writer, sheet_name="Biological_group_metadata_81", index=False)

        pd.DataFrame(
            {
                "Sample": ["C1G", "C1ITQ", "C1WORK", "M10ITW", "M10IV"],
                "ExtractionID": [1, 2, 3, 4, 5],
                "ng/ul": [385, 593, 300, 410, 220],
                "Elution_Vol": ["25  ul"] * 5,
                "Total ng": [9625, 14825, 7500, 10250, 5500],
            }
        ).to_excel(writer, sheet_name="2. DNA extraction Sample_Summar", index=False)


def test_sample_parts_and_colony_normalization() -> None:
    assert sample_parts("C1G")["colony"] == "C1"
    assert sample_parts("M10ITW")["biological_group_label"] == "In vitro worker"
    assert normalize_colony("m10") == "M10"


def test_curate_beewas_phenotypes_parses_and_writes_outputs(tmp_path: Path) -> None:
    workbook = tmp_path / "phenotypes.xlsx"
    write_fixture_workbook(workbook)
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "sample_id\tlocal_r1\tlocal_r2\n"
        + "\n".join(f"{sid}\t/tmp/{sid}_R1.fastq.gz\t/tmp/{sid}_R2.fastq.gz" for sid in ["C1G", "C1ITQ", "C1WORK", "M10ITW", "M10IV"])
        + "\n",
        encoding="utf-8",
    )

    result = curate_beewas_phenotypes(
        workbook,
        tmp_path / "curated",
        reviewed=True,
        manifest=manifest,
        repo_raw_dir=tmp_path / "repo_raw",
        blue_raw_dir=tmp_path / "blue_raw",
    )

    assert len(result.sample_phenotypes) == 5
    assert result.qc_summary["metadata_dna_exact_match"] is True
    assert result.qc_summary["manifest_sample_exact_match"] is True
    assert result.sample_phenotypes["association_ready"].all()
    assert "C1G" in result.join_report.loc[result.join_report["check"].eq("metadata_colony_column_mismatches"), "left_only"].item()
    assert result.bee_level["colony"].tolist() == ["C1", "C1", "M10"]
    assert result.paths["sample_phenotypes"].exists()
    assert result.paths["trait_dictionary"].exists()
    assert result.paths["variability_summary"].exists()
    assert result.paths["variance_components"].exists()
    assert result.paths["gwas_candidate_traits"].exists()
    assert result.paths["sample_gwas_trait_matrix"].exists()
    assert result.paths["colony_trait_matrix"].exists()
    assert result.paths["trait_readiness"].exists()
    assert result.paths["model_formulae"].exists()
    assert result.paths["model_results"].exists()
    assert result.paths["permutation_tests"].exists()
    assert result.paths["multivariate_trait_matrix"].exists()
    assert result.paths["multivariate_trait_metadata"].exists()
    assert result.paths["multivariate_tests"].exists()
    assert result.paths["multivariate_ordination"].exists()
    assert result.paths["multivariate_loadings"].exists()
    assert result.paths["multivariate_dispersion"].exists()
    assert result.paths["colony_multivariate_effects"].exists()
    assert result.paths["phenotype_knn_neighbors"].exists()
    assert result.paths["phenotype_knn_signal_summary"].exists()
    assert result.paths["phenotype_signal_dimensions"].exists()
    assert result.paths["phenotype_plaintext_statistics_md"].exists()
    assert result.paths["phenotype_plaintext_statistics_txt"].exists()
    assert result.phenotype_statistics_summary["candidate_trait_rows"] > 0
    assert "is_grafted" in set(result.gwas_candidate_traits["trait"])
    assert "is_queen_like" in set(result.gwas_candidate_traits["trait"])
    assert "gwas_trait_mode" in result.sample_gwas_trait_matrix.columns
    assert result.phenotype_statistics_summary["model_profile"] == "mixed"
    assert (tmp_path / "repo_raw" / "source_manifest.json").exists()


def test_sample_review_overrides_mark_only_requested_rows(tmp_path: Path) -> None:
    workbook = tmp_path / "phenotypes.xlsx"
    write_fixture_workbook(workbook)

    result = curate_beewas_phenotypes(
        workbook,
        tmp_path / "curated",
        reviewed=False,
        reviewed_sample_ids=["C1G", "M10ITW"],
        review_notes={"C1G": "Confirmed normalized C1 colony and grafted sample mapping."},
        manifest=None,
        repo_raw_dir=tmp_path / "repo_raw",
        blue_raw_dir=tmp_path / "blue_raw",
        phenotype_stats_profile="compact",
        phenotype_bootstrap=50,
    )

    sample = result.sample_phenotypes.set_index("sample_id")

    assert int(sample["association_ready"].sum()) == 2
    assert bool(sample.loc["C1G", "association_ready"]) is True
    assert bool(sample.loc["M10ITW", "association_ready"]) is True
    assert bool(sample.loc["C1ITQ", "association_ready"]) is False
    assert sample.loc["C1G", "phenotype_mapping_status"] == "curated_workbook_sample_map_reviewed"
    assert sample.loc["C1G", "phenotype_review_scope"] == "sample_review_override"
    assert sample.loc["C1G", "phenotype_review_note"] == "Confirmed normalized C1 colony and grafted sample mapping."
    assert result.qc_summary["reviewed_sample_override_rows"] == 2
    assert result.qc_summary["reviewed_phenotypes"] is False
    assert "phenotype_review_scope" in result.sample_gwas_trait_matrix.columns


def test_bh_fdr_and_bootstrap_are_deterministic() -> None:
    q_values = bh_fdr([0.01, 0.20, np.nan, 0.03])

    assert q_values[0] <= q_values[3] <= q_values[1]
    assert np.isnan(q_values[2])

    first = bootstrap_mean_ci(pd.Series([1.0, 2.0, 3.0, 4.0]), 100, seed=7)
    second = bootstrap_mean_ci(pd.Series([1.0, 2.0, 3.0, 4.0]), 100, seed=7)

    assert first == second
    assert first[0] <= 2.5 <= first[1]
    assert q_value_label(0.004) == "q<=0.01"
    assert standardized_mean_difference(np.array([1.0, 2.0, 3.0]), np.array([5.0, 6.0, 7.0])) < -2


def test_nested_variance_partition_small_and_estimable_cases() -> None:
    small = pd.DataFrame({"population": ["C", "C"], "colony": ["C1", "C1"], "phenotype": ["Q", "Q"], "days_to_emerge": [14, 15]})
    small_row = nested_variance_partition(small, "days_to_emerge")

    assert small_row["n"] == 2
    assert np.isnan(small_row["population_ss_fraction"])

    data = pd.DataFrame(
        {
            "population": ["C", "C", "M", "M", "M", "M"],
            "colony": ["C1", "C2", "M1", "M1", "M2", "M2"],
            "phenotype": ["Q", "Q", "Q", "Q", "Q", "Q"],
            "days_to_emerge": [14, 15, 18, 19, 17, 18],
        }
    )
    row = nested_variance_partition(data, "days_to_emerge")
    fractions = [
        row["population_ss_fraction"],
        row["colony_within_population_ss_fraction"],
        row["residual_ss_fraction"],
    ]

    assert row["n_populations"] == 2
    assert all(np.isfinite(value) for value in fractions)
    assert pytest.approx(sum(fractions), rel=1e-6, abs=1e-6) == 1.0


def test_phenotype_statistics_outputs_expected_schemas(tmp_path: Path) -> None:
    workbook = tmp_path / "phenotypes.xlsx"
    write_fixture_workbook(workbook)

    result = curate_beewas_phenotypes(
        workbook,
        tmp_path / "curated",
        reviewed=False,
        manifest=None,
        repo_raw_dir=tmp_path / "repo_raw",
        blue_raw_dir=tmp_path / "blue_raw",
        phenotype_stats_profile="compact",
        phenotype_bootstrap=50,
        phenotype_min_group_size=1,
        gwas_trait_mode="all",
    )

    assert {
        "data_level",
        "trait",
        "trait_label",
        "group_by",
        "group",
        "n",
        "missing_fraction",
        "effective_analysis_unit_count",
        "mad",
        "mean",
        "cv",
    }.issubset(result.variability_summary.columns)
    assert {"component_type", "trait", "phenotype_class", "test_status", "dominant_variance_component", "welch_anova_p_q_label"}.issubset(result.variance_components.columns)
    assert {
        "trait",
        "trait_label",
        "group_a",
        "group_b",
        "mean_diff_a_minus_b",
        "standardized_mean_diff",
        "effect_magnitude",
        "ci95_crosses_zero",
        "q_value",
        "q_value_label",
        "test_status",
        "effective_unit_warning",
    }.issubset(result.pairwise_contrasts.columns)
    assert {"trait", "trait_type", "role", "trait_tier", "trait_grain", "eligible_after_review"}.issubset(result.gwas_candidate_traits.columns)
    assert {"sample_id", "gwas_trait_mode", "phenotype_reviewed"}.issubset(result.sample_gwas_trait_matrix.columns)
    assert {
        "trait",
        "trait_tier",
        "trait_grain",
        "phenotype_data_level",
        "effective_colony_count",
        "effective_analysis_unit_count",
        "eligible_endpoint_after_review",
        "review_gate_passed",
        "statistically_estimable",
        "qc_or_covariate_only",
        "blocked_before_review",
        "blocked_reason",
        "review_gate_reason",
        "effective_unit_warning",
        "recommended_gwas_model",
    }.issubset(result.trait_readiness.columns)
    assert {"trait", "model_formula", "recommended_gwas_model", "review_status"}.issubset(result.model_formulae.columns)
    assert {"test_family", "trait", "model_status", "q_value", "q_value_label", "ci95_low", "ci95_crosses_zero"}.issubset(result.model_results.columns)
    assert {"trait", "permutation_p_value", "n_permutations", "permutation_p_q_label"}.issubset(result.permutation_tests.columns)
    assert {"fertility_status", "is_queen_like", "is_worker", "strain_fertility_status"}.issubset(result.sample_phenotypes.columns)
    assert set(result.sample_phenotypes["fertility_status"]) >= {"queen_like", "worker_like", "non_queen_like_treatment"}
    assert set(result.bee_level["fertility_status"]) == {"queen_like", "worker_like"}
    assert {"analysis_grain", "entity_id", "fertility_status"}.issubset(result.multivariate_trait_matrix.columns)
    assert {"analysis_grain", "trait", "included", "exclusion_reason", "effective_unit_warning"}.issubset(result.multivariate_trait_metadata.columns)
    assert {"analysis_grain", "test_family", "term", "test_status", "q_value_label"}.issubset(result.multivariate_tests.columns)
    assert {"analysis_grain", "entity_id", "PCoA1"}.issubset(result.multivariate_ordination.columns)
    assert {"analysis_grain", "trait", "axis", "loading", "abs_loading"}.issubset(result.multivariate_loadings.columns)
    assert {"analysis_grain", "term", "group", "mean_distance_to_centroid"}.issubset(result.multivariate_dispersion.columns)
    assert {"colony", "distance_to_strain_centroid", "distance_to_global_centroid"}.issubset(result.colony_multivariate_effects.columns)
    assert {"sample_id", "neighbor_sample_id", "distance", "k_effective", "traits_used", "same_population"}.issubset(result.phenotype_knn_neighbors.columns)
    assert {"dimension", "observed_same_fraction", "expected_same_fraction", "enrichment", "q_value_label", "signal_status"}.issubset(result.phenotype_knn_signal_summary.columns)
    assert {
        "rank",
        "within_top_n",
        "signal_score",
        "signal_family",
        "signal_dimension",
        "trait",
        "trait_label",
        "analysis_grain",
        "q_value_label",
        "map_next",
        "interpretation_guardrail",
    }.issubset(result.phenotype_signal_dimensions.columns)
    assert {"dna_ng_per_ul", "total_dna_ng", "elution_volume_ul"}.issubset(result.sample_gwas_trait_matrix.columns)
    assert result.trait_readiness.loc[result.trait_readiness["trait"].eq("dna_ng_per_ul"), "trait_tier"].item() == "qc_covariate"
    assert bool(result.trait_readiness.loc[result.trait_readiness["trait"].eq("dna_ng_per_ul"), "qc_or_covariate_only"].item()) is True
    assert bool(result.trait_readiness.loc[result.trait_readiness["trait"].eq("is_worker"), "blocked_before_review"].item()) is True
    assert result.trait_readiness.loc[result.trait_readiness["trait"].eq("days_to_emerge_mean_q"), "phenotype_data_level"].item() == "colony_grain"
    assert result.phenotype_statistics_summary["stats_profile"] == "compact"
    assert result.phenotype_statistics_summary["bootstrap_replicates"] == 50
    assert result.phenotype_statistics_summary["review_gated_endpoint_count"] > 0
    assert result.phenotype_statistics_summary["multivariate_permutations"] == 99
    assert result.phenotype_statistics_summary["multivariate_trait_matrix_rows"] == len(result.multivariate_trait_matrix)
    assert result.phenotype_statistics_summary["phenotype_knn_signal_rows"] == len(result.phenotype_knn_signal_summary)
    assert result.phenotype_statistics_summary["phenotype_signal_dimension_rows"] == len(result.phenotype_signal_dimensions)
    assert "Phenotype Statistical Handoff" in result.paths["phenotype_plaintext_statistics_md"].read_text(encoding="utf-8")

    direct = build_phenotype_statistics(
        result.bee_level,
        result.sample_phenotypes,
        reviewed=False,
        stats_profile="compact",
        n_bootstrap=50,
        min_group_size=1,
        gwas_trait_mode="all",
        phenotype_model_profile="mixed",
    )
    assert len(direct) == 22
    assert direct[-1]["candidate_trait_rows"] == len(result.gwas_candidate_traits)
    assert direct[-1]["multivariate_test_rows"] == len(result.multivariate_tests)
    assert direct[-1]["phenotype_signal_dimension_rows"] == len(result.phenotype_signal_dimensions)


def test_trait_readiness_and_model_helpers_are_deterministic(tmp_path: Path) -> None:
    workbook = tmp_path / "phenotypes.xlsx"
    write_fixture_workbook(workbook)
    result = curate_beewas_phenotypes(
        workbook,
        tmp_path / "curated",
        reviewed=False,
        manifest=None,
        repo_raw_dir=tmp_path / "repo_raw",
        blue_raw_dir=tmp_path / "blue_raw",
        phenotype_stats_profile="compact",
        phenotype_bootstrap=50,
        phenotype_min_group_size=1,
        gwas_trait_mode="all",
    )

    readiness = build_trait_readiness(result.sample_phenotypes, result.gwas_candidate_traits, reviewed=False)
    assert readiness["phenotype_reviewed"].eq(False).all()
    assert bool(readiness.loc[readiness["trait"].eq("dna_ng_per_ul"), "eligible_after_review"].item()) is False
    assert bool(readiness.loc[readiness["trait"].eq("dna_ng_per_ul"), "qc_or_covariate_only"].item()) is True
    assert bool(readiness.loc[readiness["trait"].eq("is_worker"), "statistically_estimable"].item()) is False
    assert "weak_binary_balance" in readiness.loc[readiness["trait"].eq("is_worker"), "blocked_reason"].item()
    assert readiness.loc[readiness["trait"].eq("is_worker"), "review_gate_reason"].item() == "human_review_required"
    assert readiness.loc[readiness["trait"].eq("is_worker"), "recommended_gwas_model"].item().startswith("binary/logistic")
    assert readiness.loc[readiness["trait"].eq("days_to_emerge_mean_q"), "duplicated_by_colony"].item()
    assert readiness.loc[readiness["trait"].eq("days_to_emerge_mean_q"), "effective_colony_count"].item() == 2
    assert readiness.loc[readiness["trait"].eq("is_worker"), "effective_analysis_unit_count"].item() == 5

    p1 = permutation_anova_p(
        pd.Series([1.0, 2.0, 3.0, 7.0, 8.0, 9.0]),
        pd.Series(["A", "A", "A", "B", "B", "B"]),
        n_permutations=80,
        seed=11,
    )
    p2 = permutation_anova_p(
        pd.Series([1.0, 2.0, 3.0, 7.0, 8.0, 9.0]),
        pd.Series(["A", "A", "A", "B", "B", "B"]),
        n_permutations=80,
        seed=11,
    )
    assert p1 == p2
    assert np.isfinite(p1[0])
    assert 0.0 <= p1[1] <= 1.0

    model_results, permutation_tests = build_model_results(
        result.bee_level,
        result.colony_trait_matrix,
        phenotype_model_profile="mixed",
        n_bootstrap=50,
        min_group_size=1,
    )
    assert not model_results.empty
    assert "model_status" in model_results.columns
    assert not permutation_tests.empty


def test_phenotype_knn_outputs_are_deterministic(tmp_path: Path) -> None:
    workbook = tmp_path / "phenotypes.xlsx"
    write_fixture_workbook(workbook)
    result = curate_beewas_phenotypes(
        workbook,
        tmp_path / "curated",
        reviewed=False,
        manifest=None,
        repo_raw_dir=tmp_path / "repo_raw",
        blue_raw_dir=tmp_path / "blue_raw",
        phenotype_stats_profile="compact",
        phenotype_bootstrap=50,
        phenotype_min_group_size=1,
    )

    first_neighbors, first_summary = build_phenotype_knn(result.sample_phenotypes, k=2, min_group_size=1, n_permutations=25)
    second_neighbors, second_summary = build_phenotype_knn(result.sample_phenotypes, k=2, min_group_size=1, n_permutations=25)

    pd.testing.assert_frame_equal(first_neighbors, second_neighbors)
    pd.testing.assert_frame_equal(first_summary, second_summary)
    assert first_neighbors["neighbor_rank"].max() == 2
    assert first_summary["q_value_label"].notna().all()
