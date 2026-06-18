"""Tests for the BeeWAS script-local reporting and sampling helpers."""

from __future__ import annotations

import json
import math
import sys
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
    bounded_local_ld_summary,
    build_plot_records,
    hwe_exact_p,
    in_bed,
    lambda_gc_from_pvalues,
    normal_ci,
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
    write_plot_manifest,
    write_validation_outputs,
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

    plot_manifest = write_plot_manifest(tmp_path, records)
    html_path = write_html_report(
        tmp_path,
        title="BeeWAS <Report>",
        warning="Synthetic <warning>",
        summary_cards={"samples": 3},
        sections=[],
        config={"profile": "smoke"},
    )
    output_manifest = write_output_manifest(tmp_path)
    html_text = html_path.read_text(encoding="utf-8")

    assert plot_manifest["plot_count"] == 1
    assert plot_manifest["sections"][0]["id"] == "sampling_coverage"
    assert "Sampling &lt;Flow&gt;" in html_text
    assert "Controls only &lt;never biological&gt;" in html_text
    assert output_manifest["plot_manifest"] == "plot_manifest.json"
    assert any(row["path"] == "plot_manifest.json" for row in output_manifest["files"])


def test_plot_registry_assigns_phenotype_variability_section_and_source(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    plots.mkdir()
    fig, ax = plt.subplots()
    ax.plot([0, 1], [1, 2])
    fig.savefig(plots / "phenotype_variance_partition.png")
    plt.close(fig)

    records = build_plot_records(
        tmp_path,
        {},
        context="real cohort QC report",
        n_samples=81,
        n_sites=0,
        source_table_overrides={"phenotype_variance_partition": "tables/phenotype_variance_components.tsv"},
    )
    record = records[0]

    assert record.section == "phenotype_variability"
    assert record.source_table == "tables/phenotype_variance_components.tsv"


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
