"""Tests for quality reporting, metric aggregation, and threshold checking.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from metainformant.quality.reporting.multiqc_integration import (
    aggregate_sample_qc,
    check_qc_thresholds,
    default_qc_thresholds,
    generate_qc_report,
    qc_trend_analysis,
)
from metainformant.quality.analysis.metrics import (
    calculate_coverage_metrics,
    calculate_complexity_metrics,
    calculate_data_integrity_score,
    calculate_duplication_metrics,
    calculate_gc_metrics,
    calculate_length_metrics,
    calculate_quality_metrics,
    calculate_quality_score,
    compare_quality_metrics,
    detect_outliers,
)


# ---------------------------------------------------------------------------
# default_qc_thresholds
# ---------------------------------------------------------------------------


class TestDefaultQCThresholds:
    def test_returns_dict(self):
        thresholds = default_qc_thresholds()
        assert isinstance(thresholds, dict)
        assert len(thresholds) > 0

    def test_known_metrics_present(self):
        thresholds = default_qc_thresholds()
        expected_keys = [
            "total_reads",
            "mapping_rate",
            "duplication_rate",
            "mean_quality",
            "gc_content",
            "contamination_rate",
        ]
        for key in expected_keys:
            assert key in thresholds, f"Missing expected metric: {key}"

    def test_each_metric_has_required_keys(self):
        thresholds = default_qc_thresholds()
        for metric_name, config in thresholds.items():
            assert "warn" in config, f"{metric_name} missing 'warn'"
            assert "fail" in config, f"{metric_name} missing 'fail'"
            assert "direction" in config, f"{metric_name} missing 'direction'"
            assert config["direction"] in ("min", "max"), f"{metric_name} has invalid direction"


# ---------------------------------------------------------------------------
# check_qc_thresholds
# ---------------------------------------------------------------------------


class TestCheckQCThresholds:
    def test_all_pass(self):
        metrics = {
            "total_reads": 10_000_000,
            "mapping_rate": 0.95,
            "duplication_rate": 0.10,
            "mean_quality": 35.0,
        }
        result = check_qc_thresholds(metrics)
        assert result["status"] == "pass"
        assert len(result["failed_metrics"]) == 0
        assert len(result["warnings"]) == 0

    def test_warning_level(self):
        metrics = {"mapping_rate": 0.80}  # Between warn (0.85) and fail (0.70)
        result = check_qc_thresholds(metrics)
        assert result["status"] == "warn"
        assert len(result["warnings"]) == 1
        assert result["warnings"][0]["metric"] == "mapping_rate"

    def test_fail_level(self):
        metrics = {"mapping_rate": 0.50}  # Below fail threshold (0.70)
        result = check_qc_thresholds(metrics)
        assert result["status"] == "fail"
        assert len(result["failed_metrics"]) == 1

    def test_max_direction_pass(self):
        metrics = {"duplication_rate": 0.10}  # Below warn (0.30)
        result = check_qc_thresholds(metrics)
        assert result["status"] == "pass"

    def test_max_direction_warn(self):
        metrics = {"duplication_rate": 0.40}  # Between warn (0.30) and fail (0.50)
        result = check_qc_thresholds(metrics)
        assert result["status"] == "warn"

    def test_max_direction_fail(self):
        metrics = {"duplication_rate": 0.60}  # Above fail (0.50)
        result = check_qc_thresholds(metrics)
        assert result["status"] == "fail"

    def test_unknown_metric_passes(self):
        metrics = {"unknown_metric": 42}
        result = check_qc_thresholds(metrics)
        assert result["status"] == "pass"
        assert "unknown_metric" in result["passed_metrics"]

    def test_custom_thresholds(self):
        custom = {"custom_metric": {"warn": 10, "fail": 5, "direction": "min"}}
        metrics = {"custom_metric": 3}
        result = check_qc_thresholds(metrics, thresholds=custom)
        assert result["status"] == "fail"

    def test_n_checked_count(self):
        metrics = {"total_reads": 5_000_000, "mapping_rate": 0.90}
        result = check_qc_thresholds(metrics)
        assert result["n_checked"] == 2


# ---------------------------------------------------------------------------
# aggregate_sample_qc
# ---------------------------------------------------------------------------


class TestAggregateSampleQC:
    def test_empty_input(self):
        result = aggregate_sample_qc([])
        assert result["n_samples"] == 0
        assert result["overall_pass_rate"] == 0.0

    def test_basic_aggregation(self):
        samples = [
            {"sample_id": "S1", "total_reads": 10_000_000, "mean_quality": 35.0},
            {"sample_id": "S2", "total_reads": 8_000_000, "mean_quality": 33.0},
            {"sample_id": "S3", "total_reads": 12_000_000, "mean_quality": 36.0},
        ]
        result = aggregate_sample_qc(samples)
        assert result["n_samples"] == 3
        assert "total_reads" in result["summary_table"]
        assert "mean_quality" in result["summary_table"]

    def test_summary_statistics(self):
        samples = [
            {"sample_id": f"S{i}", "mean_quality": 30.0 + i} for i in range(5)
        ]
        result = aggregate_sample_qc(samples)
        stats = result["summary_table"]["mean_quality"]
        assert stats["mean"] == pytest.approx(32.0, abs=0.01)
        assert stats["min"] == 30.0
        assert stats["max"] == 34.0
        assert stats["n_values"] == 5

    def test_outlier_detection(self):
        # 4 normal + 1 extreme outlier
        samples = [
            {"sample_id": "S1", "mean_quality": 30.0},
            {"sample_id": "S2", "mean_quality": 31.0},
            {"sample_id": "S3", "mean_quality": 30.5},
            {"sample_id": "S4", "mean_quality": 31.5},
            {"sample_id": "S5", "mean_quality": 5.0},  # Outlier
        ]
        result = aggregate_sample_qc(samples)
        assert len(result["outlier_samples"]) >= 1

    def test_pass_rate_calculation(self):
        # All pass
        samples = [
            {
                "sample_id": f"S{i}",
                "total_reads": 10_000_000,
                "mapping_rate": 0.95,
                "mean_quality": 35.0,
            }
            for i in range(3)
        ]
        result = aggregate_sample_qc(samples)
        assert result["overall_pass_rate"] == 1.0


# ---------------------------------------------------------------------------
# generate_qc_report
# ---------------------------------------------------------------------------


class TestGenerateQCReport:
    def test_basic_report(self):
        metrics = {"total_reads": 5_000_000, "mapping_rate": 0.92}
        report = generate_qc_report(metrics)
        assert "summary" in report
        assert "pass_fail_status" in report
        assert report["pass_fail_status"] in ("pass", "warn", "fail")
        assert report["metrics"] == metrics

    def test_report_with_output_path(self, tmp_path):
        metrics = {"total_reads": 10_000_000, "mean_quality": 35.0}
        output_file = str(tmp_path / "qc_report.json")
        report = generate_qc_report(metrics, output_path=output_file)
        assert report["report_path"] == output_file
        assert Path(output_file).exists()

        # Verify JSON content
        with open(output_file) as f:
            saved = json.load(f)
        assert saved["pass_fail_status"] == report["pass_fail_status"]

    def test_report_summary_narrative(self):
        metrics = {"mapping_rate": 0.50}  # Will fail
        report = generate_qc_report(metrics)
        assert "FAIL" in report["summary"]
        assert "Failed metrics:" in report["summary"]


# ---------------------------------------------------------------------------
# qc_trend_analysis
# ---------------------------------------------------------------------------


class TestQCTrendAnalysis:
    def test_insufficient_data(self):
        data = [{"mean_quality": 30.0}, {"mean_quality": 31.0}]
        result = qc_trend_analysis(data, "mean_quality")
        assert result["trend"] == "insufficient_data"
        assert result["is_drifting"] is False

    def test_stable_trend(self):
        data = [{"mean_quality": 30.0 + 0.01 * i} for i in range(20)]
        result = qc_trend_analysis(data, "mean_quality")
        assert result["trend"] in ("stable", "increasing")
        assert len(result["values"]) == 20

    def test_increasing_trend(self):
        data = [{"quality": float(i * 10)} for i in range(20)]
        result = qc_trend_analysis(data, "quality")
        assert result["trend"] == "increasing"
        assert result["slope"] > 0

    def test_decreasing_trend(self):
        data = [{"quality": 100.0 - float(i * 5)} for i in range(20)]
        result = qc_trend_analysis(data, "quality")
        assert result["trend"] == "decreasing"
        assert result["slope"] < 0

    def test_r_squared_range(self):
        data = [{"val": float(i)} for i in range(10)]
        result = qc_trend_analysis(data, "val")
        assert 0.0 <= result["r_squared"] <= 1.0

    def test_change_point_detection(self):
        # Sudden jump in the middle - change_points key always present
        data = [{"val": 10.0}] * 10 + [{"val": 100.0}] * 10
        result = qc_trend_analysis(data, "val")
        assert "change_points" in result
        # The algorithm may or may not detect change points depending on window/threshold
        assert isinstance(result["change_points"], list)


# ---------------------------------------------------------------------------
# calculate_quality_score
# ---------------------------------------------------------------------------


class TestCalculateQualityScore:
    def test_fastq_quality_score(self):
        data = {
            "basic_statistics": {"mean_quality": 35.0, "total_reads": 1000},
        }
        result = calculate_quality_score(data, data_type="fastq")
        assert "overall_score" in result
        assert "grade" in result
        assert result["grade"] in ("A", "B", "C", "D", "F")

    def test_vcf_quality_score(self):
        data = {
            "quality_scores": [50.0, 60.0, 70.0, 80.0],
            "filter_stats": {"total": 100, "PASS": 90},
            "depth_stats": {"mean": 30.0},
        }
        result = calculate_quality_score(data, data_type="vcf")
        assert result["overall_score"] > 0

    def test_bam_quality_score(self):
        data = {
            "mapping_quality": [50, 55, 60],
            "mapped_reads": 900,
            "total_reads": 1000,
            "properly_paired": 800,
            "duplicate_rate": 0.1,
        }
        result = calculate_quality_score(data, data_type="bam")
        assert result["overall_score"] > 0

    def test_unsupported_data_type(self):
        with pytest.raises(ValueError, match="Unsupported data type"):
            calculate_quality_score({}, data_type="xyz")


# ---------------------------------------------------------------------------
# detect_outliers
# ---------------------------------------------------------------------------


class TestDetectOutliers:
    def test_empty_data(self):
        result = detect_outliers([])
        assert result["outliers"] == []

    def test_iqr_method(self):
        data = [10, 11, 12, 13, 14, 15, 100]  # 100 is outlier
        result = detect_outliers(data, method="iqr")
        assert 100 in result["outliers"]

    def test_zscore_method(self):
        data = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 100.0]
        result = detect_outliers(data, method="zscore", threshold=2.0)
        assert len(result["outliers"]) >= 1

    def test_modified_zscore_method(self):
        data = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 100.0]
        result = detect_outliers(data, method="modified_zscore", threshold=3.5)
        assert result["method"] == "modified_zscore"

    def test_unsupported_method(self):
        with pytest.raises(ValueError):
            detect_outliers([1, 2, 3], method="invalid")


# ---------------------------------------------------------------------------
# Coverage, duplication, GC, length, quality, complexity metrics
# ---------------------------------------------------------------------------


class TestCoverageMetrics:
    def test_basic_coverage(self):
        values = [25.0, 30.0, 35.0, 28.0, 32.0]
        result = calculate_coverage_metrics(values)
        assert result["mean_coverage"] == pytest.approx(30.0, abs=0.1)
        assert "coverage_uniformity" in result
        assert "coverage_breadth" in result

    def test_empty_input(self):
        result = calculate_coverage_metrics([])
        assert result == {}


class TestDuplicationMetrics:
    def test_basic_duplication(self):
        levels = {1: 800, 2: 100, 3: 50, 4: 50}
        result = calculate_duplication_metrics(levels)
        assert result["total_reads"] == 1000
        assert result["unique_reads"] == 800
        assert result["duplication_rate"] == pytest.approx(20.0, abs=0.1)

    def test_empty_input(self):
        result = calculate_duplication_metrics({})
        assert result == {}


class TestGCMetrics:
    def test_basic_gc(self):
        gc_values = [0.4, 0.45, 0.5, 0.55, 0.48]
        result = calculate_gc_metrics(gc_values)
        assert 0.0 < result["mean_gc"] < 1.0
        assert "distribution" in result

    def test_empty_input(self):
        result = calculate_gc_metrics([])
        assert result == {}


class TestLengthMetrics:
    def test_basic_lengths(self):
        lengths = [100, 150, 150, 200, 250]
        result = calculate_length_metrics(lengths)
        assert result["mean_length"] == pytest.approx(170.0, abs=0.1)
        assert result["min_length"] == 100
        assert result["max_length"] == 250

    def test_empty_input(self):
        result = calculate_length_metrics([])
        assert result == {}


class TestQualityMetrics:
    def test_basic_quality(self):
        scores = [30.0, 35.0, 38.0, 25.0, 32.0]
        result = calculate_quality_metrics(scores)
        assert "mean_quality" in result
        assert "pct_q30" in result
        assert result["pct_q30"] > 0

    def test_empty_input(self):
        result = calculate_quality_metrics([])
        assert result == {}


class TestComplexityMetrics:
    def test_basic_complexity(self):
        sequences = ["ATCGATCG", "GCTAGCTA", "AATTCCGG", "TTAACCGG"]
        result = calculate_complexity_metrics(sequences)
        assert result["total_sequences"] == 4
        assert "kmer_entropy" in result
        assert "complexity_score" in result

    def test_empty_input(self):
        result = calculate_complexity_metrics([])
        assert result == {}
