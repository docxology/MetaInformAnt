"""Tests for the BeeWAS script-local reporting and sampling helpers."""

from __future__ import annotations

import json
import math
import os
import sys
import time
from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


PIPELINE_DIR = Path(__file__).resolve().parents[1] / "scripts" / "gwas" / "pipelines"
if str(PIPELINE_DIR) not in sys.path:
    sys.path.insert(0, str(PIPELINE_DIR))

from beewas_reporting import (  # noqa: E402
    PlotRecord,
    SamplingConfig,
    ValidationConfig,
    VcfMatrix,
    bounded_local_ld_summary,
    build_plot_records,
    compute_pairwise_distance,
    hwe_exact_p,
    in_bed,
    lambda_gc_from_pvalues,
    normal_ci,
    population_genetics_distance_partition,
    population_genetics_group_summary,
    population_genetics_pairwise_fst,
    population_genetics_signal_dimensions,
    qq_envelope,
    read_bed,
    read_site_list,
    sanitize_json,
    select_samples,
    select_variant_rows,
    sidak_threshold,
    target_size,
    validate_output_artifacts,
    validate_real_readiness,
    validate_sampling_config,
    validate_synthetic_summary,
    validation_config_from_args,
    visual_profile_allows,
    write_html_report,
    write_output_manifest,
    write_phenotype_statistics_mosaic,
    write_plot_manifest,
    write_validation_outputs,
)
from analyze_beewas_2026_real import (  # noqa: E402
    build_manifest_status,
    build_sample_processing_progress,
    build_phenotype_genotype_handoff_tables,
    final_genotype_estimators_ready,
    plot_sample_processing_progress,
)


def test_target_size_uses_smaller_count_or_percent() -> None:
    assert target_size(10, 8, 0.5) == 5
    assert target_size(10, 3, 50) == 3
    assert target_size(10, None, 0.25) == 2
    assert target_size(10, None, None) == 10


def test_sample_include_exclude_filters_apply_before_count(tmp_path: Path) -> None:
    include = tmp_path / "include.txt"
    exclude = tmp_path / "exclude.txt"
    include.write_text("A1\nA2\nB1\n", encoding="utf-8")
    exclude.write_text("A2\n", encoding="utf-8")
    config = SamplingConfig(
        include_samples=include,
        exclude_samples=exclude,
        sample_count=2,
        sampling_mode="head",
    )

    selected, table, summary = select_samples(["A1", "A2", "B1", "B2"], config)

    assert selected == ["A1", "B1"]
    assert summary["sample_filter_candidates"] == 2
    assert table.loc[table["sample_id"].eq("A2"), "excluded_by_file"].item() is True
    assert table.loc[table["sample_id"].eq("B2"), "included_by_file"].item() is False


def test_site_list_and_bed_filters_parse_expected_formats(tmp_path: Path) -> None:
    site_list = tmp_path / "sites.tsv"
    site_list.write_text("chr1\t10\nrs7\n", encoding="utf-8")
    bed = tmp_path / "sites.bed"
    bed.write_text("chr1\t0\t20\nchr2\t9\t11\n", encoding="utf-8")

    assert read_site_list(site_list) == {"chr1:10", "chr1", "rs7"}
    intervals = read_bed(bed)
    assert in_bed("chr1", 1, intervals)
    assert in_bed("chr1", 20, intervals)
    assert not in_bed("chr1", 21, intervals)
    assert in_bed("chr2", 10, intervals)


def test_stratified_site_selection_is_reproducible_and_balanced() -> None:
    variants = pd.DataFrame(
        {
            "chrom": ["chr1"] * 4 + ["chr2"] * 4,
            "maf": [0.02, 0.08, 0.16, 0.30, 0.02, 0.08, 0.16, 0.30],
        }
    )
    config = SamplingConfig(site_count=4, sampling_mode="stratified", sampling_seed=17)

    first = select_variant_rows(variants, config)
    second = select_variant_rows(variants, config)
    selected = variants.iloc[first]

    assert first == second
    assert len(first) == 4
    assert set(selected["chrom"]) == {"chr1", "chr2"}


def test_statistical_helpers_handle_small_and_typical_inputs() -> None:
    assert 0 <= hwe_exact_p(10, 5, 10) <= 1
    assert math.isnan(hwe_exact_p(0, 0, 0))
    assert 0.7 < lambda_gc_from_pvalues(np.linspace(0.01, 0.99, 99)) < 1.3
    assert sidak_threshold(0.05, 100) < 0.05

    beta = pd.Series([1.0, -0.5])
    se = pd.Series([0.1, 0.2])
    lower, upper = normal_ci(beta, se)
    assert list(lower.round(3)) == [0.804, -0.892]
    assert list(upper.round(3)) == [1.196, -0.108]

    envelope = qq_envelope(25)
    assert len(envelope) == 25
    assert {"expected_neg_log10_p", "lower_neg_log10_p", "upper_neg_log10_p"}.issubset(envelope.columns)


def test_population_genetics_helpers_are_deterministic_and_status_aware() -> None:
    samples = ["C1G", "C1WORK", "M10G", "M10WORK"]
    variants = pd.DataFrame(
        {
            "variant_index": [0, 1, 2],
            "variant_id": ["v1", "v2", "v3"],
            "chrom": ["chr1", "chr1", "chr1"],
            "pos": [10, 20, 30],
            "ref": ["A", "C", "G"],
            "alt": ["G", "T", "A"],
            "maf": [0.5, 0.5, 0.33],
            "call_rate": [1.0, 1.0, 0.75],
            "missing_rate": [0.0, 0.0, 0.25],
            "is_transition": [True, False, True],
            "hwe_p": [1.0, 0.5, 0.2],
        }
    )
    dosages = np.array(
        [
            [0.0, 0.0, 2.0, 2.0],
            [0.0, 1.0, 1.0, 2.0],
            [np.nan, 0.0, 1.0, 1.0],
        ]
    )
    matrix = VcfMatrix(
        all_samples=samples,
        samples=samples,
        sample_table=pd.DataFrame({"sample_id": samples}),
        variants=variants,
        dosages=dosages,
        imputed_dosages=np.where(np.isfinite(dosages), dosages, np.nanmean(dosages, axis=1)[:, None]),
        sampling_summary={"selected_samples": 4, "selected_sites": 3},
    )
    annotations = pd.DataFrame(
        {
            "sample_id": samples,
            "population": ["C", "C", "M", "M"],
            "colony": ["C1", "C1", "M10", "M10"],
            "biological_group_code": ["G", "WORK", "G", "WORK"],
            "fertility_status": ["queen_like", "worker_like", "queen_like", "worker_like"],
        }
    )

    group_summary = population_genetics_group_summary(matrix, annotations)
    pairwise_fst = population_genetics_pairwise_fst(matrix, annotations, dimensions=("population",), min_sites=2)
    distance_partition = population_genetics_distance_partition(compute_pairwise_distance(matrix), annotations)
    signal_dimensions = population_genetics_signal_dimensions(
        group_summary,
        pairwise_fst,
        distance_partition,
        pd.DataFrame(
            {
                "signal_family": ["model_result"],
                "signal_dimension": ["population"],
                "trait": ["q_fraction"],
                "trait_label": ["Queen fraction"],
                "signal_score": [0.9],
                "map_next": ["Correlate after review."],
            }
        ),
        top_n=5,
    )

    strain_rows = group_summary[group_summary["analysis_dimension"].eq("population")]
    fst_row = pairwise_fst.iloc[0]
    population_partition = distance_partition[
        distance_partition["analysis_dimension"].eq("population")
        & distance_partition["group_value"].eq("all")
    ]

    assert set(strain_rows["group_value"]) == {"C", "M"}
    assert {"observed_heterozygosity", "nucleotide_diversity_pi", "private_alt_site_count", "warning"}.issubset(group_summary.columns)
    assert fst_row["test_status"] == "tested"
    assert fst_row["fst"] > 0
    assert fst_row["sites_used"] == 3
    assert "association remains blocked" in fst_row["warning"]
    assert set(population_partition["scope"]) == {"within", "among"}
    assert population_partition["n_pairs"].min() >= 1
    assert {
        "rank",
        "within_top_n",
        "signal_family",
        "analysis_dimension",
        "comparison",
        "primary_metric",
        "raw_signal_score",
        "signal_score",
        "signal_level",
        "map_next",
        "interpretation_guardrail",
    }.issubset(signal_dimensions.columns)
    assert signal_dimensions["signal_score"].is_monotonic_decreasing
    assert signal_dimensions["signal_score"].between(0, 1).all()
    assert "phenotype_signal_handoff" in set(signal_dimensions["signal_family"])
    handoff = signal_dimensions[signal_dimensions["signal_family"].eq("phenotype_signal_handoff")].iloc[0]
    assert handoff["signal_score"] == 1.0
    assert handoff["raw_signal_score"] == 0.9
    assert "association remains blocked" in signal_dimensions["interpretation_guardrail"].iloc[0]


def test_sanitize_json_removes_nan_and_infinity() -> None:
    cleaned = sanitize_json({"a": np.nan, "b": np.inf, "c": [1, np.float64("-inf")], "d": Path("/tmp/x")})

    assert cleaned == {"a": None, "b": None, "c": [1, None], "d": "/tmp/x"}
    assert "NaN" not in json.dumps(cleaned)
    assert "Infinity" not in json.dumps(cleaned)


def test_visual_profile_gating() -> None:
    assert visual_profile_allows("compact", "summary")
    assert not visual_profile_allows("compact", "per_trait")
    assert visual_profile_allows("standard", "per_trait")
    assert not visual_profile_allows("standard", "full_only")
    assert visual_profile_allows("full", "full_only")


def test_validation_config_cli_and_file_overrides(tmp_path: Path) -> None:
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "validation": {
                    "min_samples": 40,
                    "profiles": {"strict": {"min_sites": 5000, "required_plot_count": 9}},
                }
            }
        ),
        encoding="utf-8",
    )
    args = SimpleNamespace(
        config=config_path,
        validation_profile="strict",
        min_validation_samples=80,
        min_validation_sites=None,
        lambda_gc_min=None,
        lambda_gc_max=None,
        strong_top_rank=2,
        required_plot_count=None,
    )

    config = validation_config_from_args(args)

    assert config.profile == "strict"
    assert config.min_samples == 80
    assert config.min_sites == 5000
    assert config.required_plot_count == 9
    assert config.strong_top_rank == 2


def test_validate_sampling_config_reports_bad_ranges() -> None:
    config = SamplingConfig(sample_pct=150, site_pct=-2, min_maf=0.8, max_missing=1.2, min_call_rate=-0.1)

    checks = validate_sampling_config(config)
    failures = {check.check_id for check in checks if check.status == "fail"}

    assert {"site_pct_range", "min_maf_range", "max_missing_range", "min_call_rate_range"}.issubset(failures)


def test_validate_synthetic_summary_checks_calibration_and_sensitivity() -> None:
    summary = pd.DataFrame(
        [
            {
                "trait": "null_continuous",
                "lambda_gc": 1.0,
                "bonferroni_significant_count": 0,
                "fdr_significant_count": 0,
                "n_truth_variants": 0,
                "best_truth_rank": np.nan,
                "effect_direction_match_rate": np.nan,
            },
            {
                "trait": "spike_continuous_strong",
                "lambda_gc": 1.1,
                "bonferroni_significant_count": 1,
                "fdr_significant_count": 1,
                "n_truth_variants": 1,
                "best_truth_rank": 1,
                "effect_direction_match_rate": 1.0,
            },
            {
                "trait": "pc_confounded_trait_unadjusted",
                "lambda_gc": 2.0,
                "bonferroni_significant_count": 2,
                "fdr_significant_count": 2,
                "n_truth_variants": 0,
                "best_truth_rank": np.nan,
                "effect_direction_match_rate": np.nan,
            },
            {
                "trait": "pc_confounded_trait_pc_adjusted",
                "lambda_gc": 1.2,
                "bonferroni_significant_count": 1,
                "fdr_significant_count": 1,
                "n_truth_variants": 0,
                "best_truth_rank": np.nan,
                "effect_direction_match_rate": np.nan,
            },
        ]
    )

    checks = validate_synthetic_summary(summary, n_samples=60, n_sites=2000, config=ValidationConfig())
    by_id = {check.check_id: check.status for check in checks}

    assert by_id["null_continuous_lambda_gc_calibrated"] == "pass"
    assert by_id["spike_continuous_strong_truth_rank"] == "pass"
    assert by_id["pc_adjustment_lambda_delta"] == "pass"
    assert by_id["synthetic_sample_count"] == "pass"


def test_validate_real_readiness_preserves_blocked_boundary() -> None:
    summary = {
        "complete_pairs": 80,
        "manifest_samples": 81,
        "aligned_crams": 80,
        "gwas_association_status": "blocked: reviewed phenotype mapping required",
    }
    readiness = pd.DataFrame(
        [
            {"check": "association_ready", "pass": False, "note": "blocked intentionally"},
            {"check": "reference_indexed", "pass": True, "note": "ok"},
        ]
    )

    checks = validate_real_readiness(
        summary,
        readiness,
        sampling_summary={"selected_samples": 80, "selected_sites": 10000},
        config=ValidationConfig(min_samples=80, min_sites=10000),
    )
    by_id = {check.check_id: check.status for check in checks}

    assert by_id["real_association_boundary_declared"] == "pass"
    assert by_id["readiness_association_not_premature"] == "pass"


def test_final_genotype_estimator_readiness_reports_not_run_statuses() -> None:
    readiness = pd.DataFrame(
        {
            "estimator": [
                "vcf_counts",
                "sample_join",
                "bcftools_stats",
                "plink_pgen",
                "allele_frequency",
                "missingness",
                "heterozygosity",
                "hwe",
                "pca",
                "king_relatedness_table",
                "king_relatedness_matrix",
                "relationship_matrix",
                "ld_pruned_variants",
                "fst_by_strain",
                "ld_decay",
                "roh_screening",
            ],
            "status": ["not_run"] * 16,
        }
    )

    ready, note = final_genotype_estimators_ready(readiness)

    assert ready is False
    assert "vcf_counts=not_run" in note
    assert "roh_screening=not_run" in note


def test_plot_manifest_html_and_output_manifest(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    plots.mkdir()
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    png = plots / "sampling_flow.png"
    fig.savefig(png)
    plt.close(fig)

    records = [
        PlotRecord(
            key="sampling_flow",
            section="sampling_coverage",
            title="Sampling <Flow>",
            subtitle="Synthetic & real context",
            path="plots/sampling_flow.png",
            n_samples=3,
            n_sites=5,
            source_table="tables/sampling_flow.tsv",
            warning="Controls only <never biological>",
        )
    ]
    table_path = tmp_path / "tables" / "missing_values.tsv"
    table_path.parent.mkdir()
    table_path.write_text("sample_id\tmetric\tnote\nA\t\t\nB\t1.25\tok\n", encoding="utf-8")

    plot_manifest = write_plot_manifest(tmp_path, records)
    html_path = write_html_report(
        tmp_path,
        title="BeeWAS <Report>",
        warning="Synthetic <warning>",
        summary_cards={"samples": 3},
        sections=[{"title": "Missing Preview", "body": "Blank cells stay blank.", "table": table_path}],
        config={"profile": "smoke"},
    )
    output_manifest = write_output_manifest(tmp_path)
    html_text = html_path.read_text(encoding="utf-8")

    assert plot_manifest["plot_count"] == 1
    assert plot_manifest["sections"][0]["id"] == "sampling_coverage"
    assert "Sampling &lt;Flow&gt;" in html_text
    assert "Controls only &lt;never biological&gt;" in html_text
    assert "<td>NaN</td>" not in html_text
    assert "<td>inf</td>" not in html_text.lower()
    assert output_manifest["plot_manifest"] == "plot_manifest.json"
    assert any(row["path"] == "plot_manifest.json" for row in output_manifest["files"])


def test_plot_registry_assigns_phenotype_variability_section_and_source(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    plots.mkdir()
    for name in [
        "phenotype_variance_partition.png",
        "phenotype_pairwise_contrast_forest.png",
        "phenotype_model_readiness_heatmap.png",
        "phenotype_gwas_candidate_overview.png",
        "phenotype_multivariate_pcoa.png",
        "phenotype_multivariate_test_effects.png",
        "phenotype_multivariate_trait_loadings.png",
        "phenotype_colony_multivariate_profile.png",
        "phenotype_multivariate_dispersion.png",
        "population_genetics_heterozygosity_by_group.png",
        "population_genetics_pairwise_fst_heatmap.png",
        "population_genetics_distance_partition.png",
        "population_genetics_signal_overview.png",
    ]:
        fig, ax = plt.subplots()
        ax.plot([0, 1], [1, 2])
        fig.savefig(plots / name)
        plt.close(fig)

    records = build_plot_records(
        tmp_path,
        {},
        context="real cohort QC report",
        n_samples=81,
        n_sites=0,
        warning="Real association remains blocked until reviewed phenotype handoff.",
        source_table_overrides={
            "phenotype_variance_partition": "tables/phenotype_variance_components.tsv",
            "phenotype_pairwise_contrast_forest": "tables/phenotype_pairwise_contrasts.tsv",
            "phenotype_model_readiness_heatmap": "tables/phenotype_trait_readiness.tsv",
            "phenotype_gwas_candidate_overview": "tables/phenotype_gwas_candidate_traits.tsv",
            "phenotype_multivariate_pcoa": "tables/phenotype_multivariate_ordination.tsv",
            "phenotype_multivariate_test_effects": "tables/phenotype_multivariate_tests.tsv",
            "phenotype_multivariate_trait_loadings": "tables/phenotype_multivariate_loadings.tsv",
            "phenotype_colony_multivariate_profile": "tables/phenotype_colony_multivariate_effects.tsv",
            "phenotype_multivariate_dispersion": "tables/phenotype_multivariate_dispersion.tsv",
            "population_genetics_heterozygosity_by_group": "tables/population_genetics_group_summary.tsv",
            "population_genetics_pairwise_fst_heatmap": "tables/population_genetics_pairwise_fst.tsv",
            "population_genetics_distance_partition": "tables/population_genetics_distance_partition.tsv",
            "population_genetics_signal_overview": "tables/population_genetics_signal_dimensions.tsv",
        },
    )
    by_key = {record.key: record for record in records}

    assert by_key["phenotype_variance_partition"].section == "phenotype_variability"
    assert by_key["phenotype_pairwise_contrast_forest"].section == "phenotype_variability"
    assert by_key["phenotype_model_readiness_heatmap"].source_table == "tables/phenotype_trait_readiness.tsv"
    assert by_key["phenotype_gwas_candidate_overview"].source_table == "tables/phenotype_gwas_candidate_traits.tsv"
    assert by_key["phenotype_multivariate_pcoa"].section == "phenotype_variability"
    assert by_key["phenotype_multivariate_test_effects"].source_table == "tables/phenotype_multivariate_tests.tsv"
    assert by_key["phenotype_multivariate_trait_loadings"].source_table == "tables/phenotype_multivariate_loadings.tsv"
    assert by_key["phenotype_colony_multivariate_profile"].source_table == "tables/phenotype_colony_multivariate_effects.tsv"
    assert by_key["phenotype_multivariate_dispersion"].source_table == "tables/phenotype_multivariate_dispersion.tsv"
    assert by_key["population_genetics_heterozygosity_by_group"].section == "population_structure"
    assert by_key["population_genetics_pairwise_fst_heatmap"].source_table == "tables/population_genetics_pairwise_fst.tsv"
    assert by_key["population_genetics_distance_partition"].source_table == "tables/population_genetics_distance_partition.tsv"
    assert by_key["population_genetics_signal_overview"].section == "population_structure"
    assert by_key["population_genetics_signal_overview"].source_table == "tables/population_genetics_signal_dimensions.tsv"
    assert all(record.warning == "Real association remains blocked until reviewed phenotype handoff." for record in by_key.values())


def test_sample_processing_progress_outputs_status_tables_and_plot(tmp_path: Path) -> None:
    base = tmp_path / "beewas"
    tables = tmp_path / "report" / "tables"
    plots = tmp_path / "report" / "plots"
    cram_dir = base / "aligned" / "cram"
    log_dir = base / "logs"
    tmp_dir = base / "tmp"
    for directory in [tables, plots, cram_dir, log_dir, tmp_dir]:
        directory.mkdir(parents=True)
    (cram_dir / "C1G.sorted.cram").write_bytes(b"cram")
    (cram_dir / "C1G.sorted.cram.crai").write_bytes(b"crai")
    (cram_dir / "C1ITQ.sorted.cram.tmp").write_bytes(b"tmp")
    (log_dir / "C1ITQ.bwa.log").write_text("[M::mem_process_seqs] Processed 100 reads\n", encoding="utf-8")
    (log_dir / "C1ITQ.sort.log").write_text("[bam_sort_core] merging\n", encoding="utf-8")
    r14_log = log_dir / "R14G.bwa.log"
    r14_log.write_text("[W::bseq_read] the 2nd file has fewer sequences.\n[gzclose] buffer error\n", encoding="utf-8")
    old_mtime = time.time() - 7200
    os.utime(r14_log, (old_mtime, old_mtime))
    manifest_status = pd.DataFrame(
        [
            {
                "sample_id": "C1G",
                "strain": "C",
                "colony": "C1",
                "library_label": "G",
                "status": "complete_pair",
                "r1_complete": True,
                "r2_complete": True,
                "r1_size_gib": 1.0,
                "r2_size_gib": 1.0,
                "cram_exists": True,
                "crai_exists": True,
                "cram_size_gib": 0.5,
            },
            {
                "sample_id": "C1ITQ",
                "strain": "C",
                "colony": "C1",
                "library_label": "ITQ",
                "status": "complete_pair",
                "r1_complete": True,
                "r2_complete": True,
                "r1_size_gib": 1.0,
                "r2_size_gib": 1.0,
                "cram_exists": False,
                "crai_exists": False,
                "cram_size_gib": 0.0,
            },
            {
                "sample_id": "M10G",
                "strain": "M",
                "colony": "M10",
                "library_label": "G",
                "status": "complete_pair",
                "r1_complete": True,
                "r2_complete": True,
                "r1_size_gib": 1.0,
                "r2_size_gib": 1.0,
                "cram_exists": False,
                "crai_exists": False,
                "cram_size_gib": 0.0,
            },
            {
                "sample_id": "M6ITQ",
                "strain": "M",
                "colony": "M6",
                "library_label": "ITQ",
                "status": "upstream_fastq_replacement_required",
                "r1_complete": False,
                "r2_complete": False,
                "r1_size_gib": 1.0,
                "r2_size_gib": 1.0,
                "r1_integrity_status": "corrupt_gzip",
                "r2_integrity_status": "ok",
                "r1_repair_action": "redownload_or_replace_fastq",
                "r2_repair_action": "none",
                "r1_integrity_message": "unexpected end of file",
                "r2_integrity_message": "",
                "fastq_integrity_checked_mate_count": 2,
                "fastq_integrity_ok_mate_count": 1,
                "fastq_integrity_failed_mate_count": 1,
                "fastq_integrity_failed_mates": "R1",
                "fastq_integrity_repair_actions": "R1:redownload_or_replace_fastq",
                "r1_download_status": "failed_integrity_corrupt_gzip",
                "r2_download_status": "not_checked",
                "r1_download_repair_action": "replace_upstream_fastq",
                "r2_download_repair_action": "none",
                "r1_download_message": "unexpected end of file after direct download",
                "r2_download_message": "",
                "fastq_download_repair_actions": "R1:replace_upstream_fastq",
                "upstream_replacement_mates": "R1",
                "cram_exists": False,
                "crai_exists": False,
                "cram_size_gib": 0.0,
            },
            {
                "sample_id": "R9WORK",
                "strain": "R",
                "colony": "R9",
                "library_label": "WORK",
                "status": "missing_pair",
                "r1_complete": False,
                "r2_complete": False,
                "r1_size_gib": 0.0,
                "r2_size_gib": 0.0,
                "cram_exists": False,
                "crai_exists": False,
                "cram_size_gib": 0.0,
            },
            {
                "sample_id": "R14G",
                "strain": "R",
                "colony": "R14",
                "library_label": "G",
                "status": "complete_pair",
                "r1_complete": True,
                "r2_complete": True,
                "r1_size_gib": 1.0,
                "r2_size_gib": 1.0,
                "cram_exists": False,
                "crai_exists": False,
                "cram_size_gib": 0.0,
            },
        ]
    )
    paths = SimpleNamespace(base=base, tables=tables, plots=plots, output_dir=tmp_path / "report")

    progress, summary = build_sample_processing_progress(paths, manifest_status)
    plot_sample_processing_progress(progress, paths)
    records = build_plot_records(
        paths.output_dir,
        {},
        context="real cohort QC report",
        n_samples=4,
        n_sites=0,
        warning="Real association remains blocked.",
        section_overrides={"sample_processing_progress": "sampling_coverage"},
        source_table_overrides={"sample_processing_progress": "tables/sample_processing_progress.tsv"},
    )
    by_key = {record.key: record for record in records}

    by_sample = progress.set_index("sample_id")
    assert by_sample.loc["C1G", "alignment_stage"] == "cram_and_index_complete"
    assert by_sample.loc["C1ITQ", "alignment_stage"] == "alignment_in_progress"
    assert by_sample.loc["M10G", "alignment_stage"] == "awaiting_alignment"
    assert by_sample.loc["M6ITQ", "alignment_stage"] == "upstream_fastq_replacement_required"
    assert by_sample.loc["M6ITQ", "blocking_reason"] == "validated_download_integrity_failed_replace_source_fastq"
    assert by_sample.loc["M6ITQ", "fastq_integrity_failed_mates"] == "R1"
    assert by_sample.loc["M6ITQ", "fastq_integrity_repair_actions"] == "R1:redownload_or_replace_fastq"
    assert by_sample.loc["M6ITQ", "fastq_download_repair_actions"] == "R1:replace_upstream_fastq"
    assert by_sample.loc["M6ITQ", "upstream_replacement_mates"] == "R1"
    assert by_sample.loc["R14G", "alignment_stage"] == "fastq_integrity_suspected_from_alignment_log"
    assert by_sample.loc["R14G", "blocking_reason"] == "alignment_log_fastq_integrity_suspected_validate_or_redownload"
    assert "bwa_gzip_close_buffer_error" in by_sample.loc["R14G", "alignment_log_fastq_issue_code"]
    assert bool(by_sample.loc["R9WORK", "needs_operator_action"]) is True
    assert summary["complete_cram_crai_pairs"] == 1
    assert summary["alignment_in_progress"] == 1
    assert summary["awaiting_alignment"] == 1
    assert summary["fastq_integrity_failed"] == 0
    assert summary["upstream_fastq_replacement_required"] == 1
    assert summary["fastq_integrity_checked_mates"] == 2
    assert summary["fastq_integrity_ok_mates"] == 1
    assert summary["fastq_integrity_failed_mates"] == 1
    assert summary["repair_action_counts"] == {"redownload_or_replace_fastq": 1}
    assert summary["download_repair_action_counts"] == {"replace_upstream_fastq": 1}
    assert summary["redownload_required_samples"] == []
    assert summary["fastq_integrity_failed_samples"] == []
    assert summary["upstream_replacement_required_samples"] == ["M6ITQ"]
    assert summary["fastq_integrity_suspected_from_logs"] == 1
    assert summary["fastq_integrity_suspected_from_log_samples"] == ["R14G"]
    assert (tables / "sample_processing_progress.tsv").exists()
    report_text = (tables / "sample_processing_progress_report.md").read_text(encoding="utf-8")
    assert report_text.startswith("# BeeWAS Sample Processing Progress")
    assert "FASTQ Integrity Repair Actions" in report_text
    assert "Download Repair Actions" in report_text
    assert "FASTQ integrity failed mates: 1" in report_text
    assert "R1:redownload_or_replace_fastq" in report_text
    assert "R1:replace_upstream_fastq" in report_text
    assert "Upstream FASTQ replacement required: 1" in report_text
    assert "FASTQ Integrity Failures" in report_text
    assert "Alignment-Log FASTQ Integrity Suspects" in report_text
    assert "unexpected end of file" in report_text
    assert "bwa_gzip_close_buffer_error" in report_text
    assert (paths.output_dir / "sample_processing_progress_summary.json").exists()
    assert (plots / "sample_processing_progress.png").exists()
    assert by_key["sample_processing_progress"].section == "sampling_coverage"
    assert by_key["sample_processing_progress"].source_table == "tables/sample_processing_progress.tsv"


def test_manifest_status_uses_fastq_integrity_status_file(tmp_path: Path) -> None:
    base = tmp_path / "beewas"
    tables = tmp_path / "report" / "tables"
    fastq = base / "raw" / "fastq"
    manifests = base / "manifests"
    fastq.mkdir(parents=True)
    tables.mkdir(parents=True)
    manifests.mkdir(parents=True)
    r1 = fastq / "M6ITQ_R1.fastq.gz"
    r2 = fastq / "M6ITQ_R2.fastq.gz"
    r9_r1 = fastq / "R9WORK_R1.fastq.gz"
    r9_r2 = fastq / "R9WORK_R2.fastq.gz"
    r1.write_bytes(b"\x1f\x8btruncated")
    r2.write_bytes(b"\x1f\x8bok")
    manifest = manifests / "manifest.tsv"
    manifest.write_text(
        "sample_id\tlocal_r1\tlocal_r2\n"
        f"M6ITQ\t{r1}\t{r2}\n"
        f"R9WORK\t{r9_r1}\t{r9_r2}\n",
        encoding="utf-8",
    )
    (manifests / "fastq_integrity_status.tsv").write_text(
        "timestamp_utc\tsample_id\tmate\tpath\tintegrity_status\trepair_action\tmessage\n"
        f"2026-06-19T05:00:00+00:00\tM6ITQ\tR1\t{r1}\tcorrupt_gzip\tredownload_or_replace_fastq\tunexpected end of file\n"
        f"2026-06-19T05:00:00+00:00\tM6ITQ\tR2\t{r2}\tok\tnone\t\n"
        f"2026-06-21T20:50:00+00:00\tR9WORK\tR1\t{r9_r1}\tmissing\tredownload_or_replace_fastq\tdirect repair failed\n"
        f"2026-06-21T20:50:00+00:00\tR9WORK\tR2\t{r9_r2}\tmissing\tredownload_or_replace_fastq\tdirect repair failed\n",
        encoding="utf-8",
    )
    (manifests / "google_drive_download_status.tsv").write_text(
        "timestamp_utc\tsample_id\tmate\tpath\tstatus\tsize_bytes\tmessage\n"
        f"2026-06-19T05:10:00+00:00\tM6ITQ\tR1\t{r1}\tfailed_integrity_corrupt_gzip\t0\tunexpected end of file after direct download\n"
        f"2026-06-21T20:51:00+00:00\tR9WORK\tR1\t{r9_r1}\tfailed_rc_1\t0\tcurl --fail -r 50331648-51380223 https://drive.usercontent.google.com/download?id=abc logs/downloads/direct_tmp/R9WORK_R1.fastq.gz\n"
        f"2026-06-21T20:52:00+00:00\tR9WORK\tR2\t{r9_r2}\tfailed_rc_1\t0\tcurl --fail -r 0-0 https://drive.usercontent.google.com/download?id=def logs/downloads/direct_tmp/R9WORK_R2.fastq.gz\n",
        encoding="utf-8",
    )
    paths = SimpleNamespace(base=base, manifest=manifest, tables=tables)

    status = build_manifest_status(paths)

    row = status.iloc[0]
    assert row["status"] == "upstream_fastq_replacement_required"
    assert bool(row["r1_complete"]) is False
    assert row["r1_integrity_status"] == "corrupt_gzip"
    assert row["r1_repair_action"] == "redownload_or_replace_fastq"
    assert row["r1_download_status"] == "failed_integrity_corrupt_gzip"
    assert row["r1_download_repair_action"] == "replace_upstream_fastq"
    assert row["upstream_replacement_mates"] == "R1"
    assert int(row["fastq_integrity_checked_mate_count"]) == 2
    assert int(row["fastq_integrity_ok_mate_count"]) == 1
    assert int(row["fastq_integrity_failed_mate_count"]) == 1
    assert row["fastq_integrity_failed_mates"] == "R1"
    assert row["fastq_integrity_repair_actions"] == "R1:redownload_or_replace_fastq"
    assert row["fastq_download_repair_actions"] == "R1:replace_upstream_fastq"
    r9 = status.set_index("sample_id").loc["R9WORK"]
    assert r9["status"] == "upstream_fastq_replacement_required"
    assert r9["r1_download_status"] == "failed_rc_1"
    assert r9["r1_download_repair_action"] == "replace_upstream_fastq"
    assert r9["r2_download_repair_action"] == "replace_upstream_fastq"
    assert r9["upstream_replacement_mates"] == "R1,R2"
    assert r9["fastq_download_repair_actions"] == "R1:replace_upstream_fastq,R2:replace_upstream_fastq"
    assert (tables / "manifest_download_alignment_status.tsv").exists()
    assert (tables / "google_drive_download_status.tsv").exists()


def test_phenotype_statistics_mosaic_is_manifest_backed(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    plots.mkdir()
    for name in [
        "phenotype_sample_join_coverage.png",
        "phenotype_variance_partition.png",
        "phenotype_multivariate_pcoa.png",
        "sample_het_vs_nonref_scatter.png",
    ]:
        fig, ax = plt.subplots(figsize=(3, 2))
        ax.plot([0, 1, 2], [1, 3, 2])
        ax.set_title(name)
        fig.savefig(plots / name)
        plt.close(fig)

    warning = "Real association remains blocked until reviewed phenotype handoff."
    section_overrides = {
        "phenotype_sample_join_coverage": "real_cohort_qc",
        "phenotype_variance_partition": "phenotype_variability",
        "phenotype_multivariate_pcoa": "phenotype_variability",
        "phenotype_multidimensional_statistics_mosaic": "phenotype_variability",
    }
    source_table_overrides = {
        "phenotype_sample_join_coverage": "tables/phenotype_join_report.tsv",
        "phenotype_variance_partition": "tables/phenotype_variance_components.tsv",
        "phenotype_multivariate_pcoa": "tables/phenotype_multivariate_ordination.tsv",
        "phenotype_multidimensional_statistics_mosaic": "tables/phenotype_graphical_statistics_overview.tsv",
    }
    records = build_plot_records(
        tmp_path,
        {},
        context="real cohort QC report",
        n_samples=81,
        n_sites=1000,
        warning=warning,
        section_overrides=section_overrides,
        source_table_overrides=source_table_overrides,
    )

    mosaic_path, overview = write_phenotype_statistics_mosaic(tmp_path, records, columns=2, width_px=1200, dpi=100)

    assert mosaic_path == tmp_path / "plots" / "phenotype_multidimensional_statistics_mosaic.png"
    assert mosaic_path.exists()
    assert (tmp_path / "tables" / "phenotype_graphical_statistics_overview.tsv").exists()
    assert {
        "mosaic_order",
        "plot_key",
        "section",
        "source_table",
        "image_path",
        "image_width_px",
        "image_height_px",
    }.issubset(overview.columns)
    assert list(overview["plot_key"]) == [
        "phenotype_sample_join_coverage",
        "phenotype_multivariate_pcoa",
        "phenotype_variance_partition",
    ]
    assert "sample_het_vs_nonref_scatter" not in set(overview["plot_key"])
    mosaic_image = plt.imread(mosaic_path)
    assert mosaic_image.shape[1] >= 1000
    assert float(np.var(mosaic_image)) > 0.0

    records_with_mosaic = build_plot_records(
        tmp_path,
        {},
        context="real cohort QC report",
        n_samples=81,
        n_sites=1000,
        warning=warning,
        section_overrides=section_overrides,
        source_table_overrides=source_table_overrides,
    )
    by_key = {record.key: record for record in records_with_mosaic}

    assert by_key["phenotype_multidimensional_statistics_mosaic"].section == "phenotype_variability"
    assert by_key["phenotype_multidimensional_statistics_mosaic"].source_table == "tables/phenotype_graphical_statistics_overview.tsv"
    assert by_key["phenotype_multidimensional_statistics_mosaic"].warning == warning


def test_phenotype_genotype_handoff_tables_are_status_aware(tmp_path: Path) -> None:
    tables = tmp_path / "tables"
    tables.mkdir()
    paths = SimpleNamespace(
        tables=tables,
        phenotype_distance_permutations=25,
        phenotype_stats_profile="compact",
        phenotype_signal_top_n=5,
    )
    samples = ["C1G", "C1ITQ", "M10ITW", "M10IV"]
    phenotype_result = SimpleNamespace(
        sample_phenotypes=pd.DataFrame(
            {
                "sample_id": samples,
                "population": ["C", "C", "M", "M"],
                "colony": ["C1", "C1", "M10", "M10"],
                "biological_group_code": ["G", "ITQ", "ITW", "IV"],
                "fertility_status": ["queen_like", "queen_like", "non_queen_like_treatment", "non_queen_like_treatment"],
                "phenotype_reviewed": [True, True, True, True],
                "phenotype_review_scope": ["sample_review_override"] * 4,
                "association_ready": [True, True, True, True],
                "q_fraction": [0.9, 0.8, 0.2, 0.1],
                "dna_ng_per_ul": [400.0, 380.0, 220.0, 210.0],
            }
        ),
        multivariate_trait_matrix=pd.DataFrame(
            {
                "analysis_grain": ["sample"] * 4,
                "entity_id": samples,
                "sample_id": samples,
                "q_fraction": [1.2, 0.9, -0.8, -1.1],
                "dna_ng_per_ul": [1.0, 0.8, -0.7, -1.1],
            }
        ),
        multivariate_trait_metadata=pd.DataFrame(
            {
                "analysis_grain": ["sample", "sample"],
                "trait": ["q_fraction", "dna_ng_per_ul"],
                "included": [True, True],
            }
        ),
        phenotype_signal_dimensions=pd.DataFrame(
            {
                "trait": ["q_fraction", "dna_ng_per_ul"],
                "rank": [1, 2],
            }
        ),
    )
    genotype_distance = pd.DataFrame(
        [
            [0.0, 0.05, 0.45, 0.50],
            [0.05, 0.0, 0.40, 0.47],
            [0.45, 0.40, 0.0, 0.08],
            [0.50, 0.47, 0.08, 0.0],
        ],
        index=samples,
        columns=samples,
    )
    pca = pd.DataFrame({"sample_id": samples, "PC1": [1.0, 0.8, -0.7, -1.1], "PC2": [0.2, -0.2, 0.3, -0.3]})
    sample_qc = pd.DataFrame({"sample": samples, "het_rate": [0.1, 0.11, 0.2, 0.22], "nonref_rate": [0.3, 0.31, 0.5, 0.52], "singletons": [1, 2, 4, 5]})

    summary = build_phenotype_genotype_handoff_tables(paths, phenotype_result, genotype_distance, pca, sample_qc)

    assert summary["reviewed_subset_samples"] == 4
    assert summary["phenotype_genotype_pair_rows"] == 6
    for name in [
        "phenotype_genotype_relatedness_handoff.tsv",
        "phenotype_genotype_distance_correlations.tsv",
        "phenotype_pc_correlations.tsv",
        "phenotype_genotype_within_among_distance_summary.tsv",
    ]:
        assert (tables / name).exists()
    distance_rows = pd.read_csv(tables / "phenotype_genotype_distance_correlations.tsv", sep="\t")
    pc_rows = pd.read_csv(tables / "phenotype_pc_correlations.tsv", sep="\t")
    assert {"rho", "p_value", "q_value", "q_value_label", "test_status", "warning"}.issubset(distance_rows.columns)
    assert distance_rows["warning"].str.contains("no biological GWAS association claim").all()
    assert {"trait", "genotype_axis", "q_value_label", "test_status"}.issubset(pc_rows.columns)


def test_validation_outputs_and_artifact_checks(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    plots.mkdir()
    fig, ax = plt.subplots()
    ax.plot([0, 1], [1, 0])
    png = plots / "qc.png"
    fig.savefig(png)
    plt.close(fig)
    (tmp_path / "analysis_summary.json").write_text(json.dumps({"ok": True}), encoding="utf-8")
    records = [PlotRecord(key="qc", section="real_cohort_qc", title="QC", subtitle="ok", path="plots/qc.png")]
    write_plot_manifest(tmp_path, records)
    write_html_report(tmp_path, title="Report", warning=None, summary_cards={}, sections=[], config={})

    checks = validate_output_artifacts(tmp_path, records, ValidationConfig(required_plot_count=1))
    summary = write_validation_outputs(tmp_path, checks, validation_config=ValidationConfig(), run_context="unit")

    assert summary["status_counts"]["pass"] >= 3
    assert (tmp_path / "tables" / "run_validation.tsv").exists()
    assert (tmp_path / "run_validation.json").exists()


def test_bounded_local_ld_summary_handles_missing_and_small_inputs() -> None:
    variants = pd.DataFrame(
        {
            "variant_index": [0, 1, 2, 3],
            "variant_id": ["v0", "v1", "v2", "v3"],
            "chrom": ["chr1", "chr1", "chr1", "chr2"],
            "pos": [10, 20, 30, 10],
            "ref": ["A", "A", "A", "A"],
            "alt": ["G", "G", "G", "G"],
        }
    )
    dosages = np.array(
        [
            [0.0, 1.0, 2.0, np.nan],
            [0.0, 1.0, 2.0, 1.0],
            [2.0, 1.0, 0.0, np.nan],
            [0.0, 0.0, 1.0, 1.0],
        ]
    )

    ld = bounded_local_ld_summary(variants, dosages, 1, max_sites=3)
    empty = bounded_local_ld_summary(variants, dosages, 999)

    assert len(ld) == 3
    assert ld.loc[ld["window_variant_index"].eq(1), "ld_r2_to_center"].iloc[0] == 1.0
    assert set(ld["chrom"]) == {"chr1"}
    assert empty.empty
