"""Tests for quality metrics and statistical analysis."""

import numpy as np
import pytest

from metainformant.quality.metrics import (
    calculate_complexity_metrics,
    calculate_coverage_metrics,
    calculate_duplication_metrics,
    calculate_gc_metrics,
    calculate_length_metrics,
    calculate_quality_metrics,
    generate_quality_report,
)


class TestQualityMetrics:
    """Test quality metrics calculation."""

    def test_calculate_quality_metrics_basic(self):
        """Test basic quality metrics calculation."""
        quality_scores = [[30, 35, 32], [28, 31, 29]]

        metrics = calculate_quality_metrics(quality_scores)

        assert "mean_quality" in metrics
        assert "median_quality" in metrics
        assert "pct_q20" in metrics
        assert "pct_q30" in metrics
        assert metrics["mean_quality"] > 0
        assert metrics["pct_q20"] == 100.0  # All scores >= 20

    def test_calculate_quality_metrics_empty(self):
        """Test handling of empty input."""
        metrics = calculate_quality_metrics([])

        assert metrics == {}

    def test_calculate_quality_metrics_single_score(self):
        """Test with single quality score."""
        quality_scores = [[25]]

        metrics = calculate_quality_metrics(quality_scores)

        assert metrics["mean_quality"] == 25.0
        assert metrics["pct_q20"] == 100.0
        assert metrics["pct_q30"] == 0.0  # 25 < 30


class TestGCMetrics:
    """Test GC content metrics calculation."""

    def test_calculate_gc_metrics_basic(self):
        """Test basic GC metrics calculation."""
        gc_content = [45.0, 50.0, 55.0, 48.0]

        metrics = calculate_gc_metrics(gc_content)

        assert "mean_gc" in metrics
        assert "std_gc" in metrics
        assert "gc_bias" in metrics
        assert metrics["mean_gc"] == 49.5
        assert metrics["gc_bias"] == 0.5  # Deviation from expected 50%

    def test_calculate_gc_metrics_empty(self):
        """Test handling of empty GC content."""
        metrics = calculate_gc_metrics([])

        assert metrics == {}

    def test_calculate_gc_metrics_single_value(self):
        """Test with single GC value."""
        gc_content = [60.0]

        metrics = calculate_gc_metrics(gc_content)

        assert metrics["mean_gc"] == 60.0
        assert metrics["gc_bias"] == 10.0  # Deviation from expected 50%


class TestLengthMetrics:
    """Test sequence length metrics calculation."""

    def test_calculate_length_metrics_basic(self):
        """Test basic length metrics calculation."""
        lengths = [100, 150, 200, 125]

        metrics = calculate_length_metrics(lengths)

        assert "mean_length" in metrics
        assert "std_length" in metrics
        assert "pct_short_reads" in metrics
        assert "pct_long_reads" in metrics
        assert metrics["mean_length"] == 143.75
        assert metrics["pct_short_reads"] == 0.0  # All lengths >= 50
        assert metrics["pct_long_reads"] == 0.0  # None of [100, 150, 200, 125] are > 1000

    def test_calculate_length_metrics_empty(self):
        """Test handling of empty lengths."""
        metrics = calculate_length_metrics([])

        assert metrics == {}

    def test_calculate_length_metrics_with_short_reads(self):
        """Test with short reads."""
        lengths = [30, 40, 50, 100]

        metrics = calculate_length_metrics(lengths)

        assert metrics["pct_short_reads"] == 50.0  # 30, 40 < 50


class TestDuplicationMetrics:
    """Test duplication metrics calculation."""

    def test_calculate_duplication_metrics_basic(self):
        """Test basic duplication metrics calculation."""
        duplication_levels = {1: 1000, 2: 500, 3: 200}

        metrics = calculate_duplication_metrics(duplication_levels)

        assert "duplication_rate" in metrics
        assert "unique_rate" in metrics
        assert "mean_duplication_level" in metrics

        total_reads = 1700
        unique_reads = 1000
        expected_duplication_rate = (total_reads - unique_reads) / total_reads * 100
        assert abs(metrics["duplication_rate"] - expected_duplication_rate) < 1e-10

    def test_calculate_duplication_metrics_empty(self):
        """Test handling of empty duplication data."""
        metrics = calculate_duplication_metrics({})

        assert metrics == {}

    def test_calculate_duplication_metrics_all_unique(self):
        """Test with all unique reads."""
        duplication_levels = {1: 1000}

        metrics = calculate_duplication_metrics(duplication_levels)

        assert metrics["duplication_rate"] == 0.0
        assert metrics["unique_rate"] == 100.0


class TestComplexityMetrics:
    """Test complexity metrics calculation."""

    def test_calculate_complexity_metrics_basic(self):
        """Test basic complexity metrics calculation."""
        sequences = ["ATCG", "AAAA", "ACGTACGT"]

        metrics = calculate_complexity_metrics(sequences)

        assert "mean_complexity" in metrics
        assert "low_complexity_rate" in metrics

        # AAAA should have low complexity (4 unique chars / 4 total = 1.0)
        # ATCG should have high complexity (4 unique chars / 4 total = 1.0)
        # ACGTACGT should have high complexity (6 unique chars / 8 total = 0.75)

    def test_calculate_complexity_metrics_empty(self):
        """Test handling of empty sequences."""
        metrics = calculate_complexity_metrics([])

        assert metrics == {}

    def test_calculate_complexity_metrics_identical(self):
        """Test with identical sequences."""
        sequences = ["ATCG", "ATCG", "ATCG"]

        metrics = calculate_complexity_metrics(sequences)

        assert metrics["mean_complexity"] == 1.0  # All sequences have max complexity


class TestCoverageMetrics:
    """Test coverage metrics calculation."""

    def test_calculate_coverage_metrics_basic(self):
        """Test basic coverage metrics calculation."""
        coverage_values = [25.0, 30.0, 35.0, 28.0]

        metrics = calculate_coverage_metrics(coverage_values, target_coverage=30.0)

        assert "mean_coverage" in metrics
        assert "coverage_bias" in metrics
        assert "pct_low_coverage" in metrics
        assert "pct_high_coverage" in metrics

        assert metrics["mean_coverage"] == 29.5
        assert metrics["coverage_bias"] == 0.5  # Deviation from target 30.0

    def test_calculate_coverage_metrics_empty(self):
        """Test handling of empty coverage data."""
        metrics = calculate_coverage_metrics([])

        assert metrics == {}


class TestQualityReport:
    """Test quality report generation."""

    def test_generate_quality_report_basic(self):
        """Test basic quality report generation."""
        quality_data = {
            "quality": {"mean_quality": 35.0, "pct_q30": 85.0},
            "gc_content": {"mean_gc": 48.0, "gc_bias": 2.0},
            "length": {"mean_length": 150.0},
            "duplication": {"duplication_rate": 15.0},
            "complexity": {"mean_complexity": 0.8},
        }

        report = generate_quality_report(quality_data, "Sample_001")

        assert "METAINFORMANT Quality Control Report" in report
        assert "Sample: Sample_001" in report
        assert "Quality Metrics:" in report
        assert "GC Content Metrics:" in report
        assert "Overall Quality Score:" in report

    def test_generate_quality_report_empty(self):
        """Test report generation with no quality data."""
        report = generate_quality_report({}, "Empty_Sample")

        assert "METAINFORMANT Quality Control Report" in report
        assert "Sample: Empty_Sample" in report

    def test_generate_quality_report_excellent_quality(self):
        """Test report generation with excellent quality scores."""
        quality_data = {
            "quality": {"mean_quality": 38.0, "pct_q30": 95.0},
            "gc_content": {"mean_gc": 50.0, "gc_bias": 0.0},
            "duplication": {"duplication_rate": 5.0},
            "complexity": {"mean_complexity": 0.95},
        }

        report = generate_quality_report(quality_data, "Excellent_Sample")

        assert "Overall Quality Score: 10.0" in report
        assert "Assessment: EXCELLENT" in report
