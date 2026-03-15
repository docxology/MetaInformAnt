"""Tests for longread quality summary functions and QCReport not covered by test_longread.py.

Tests generate_qc_summary, QCReport dataclass, and related reporting functions.
All tests use real implementations -- NO MOCKING.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from metainformant.longread.utils.summary import generate_qc_summary
from metainformant.longread.workflow.reporting import (
    QCReport,
    export_report,
    generate_qc_report,
    generate_run_summary,
)

# ---------------------------------------------------------------------------
# generate_qc_summary
# ---------------------------------------------------------------------------


class TestGenerateQcSummary:
    """Tests for generate_qc_summary function."""

    def test_qc_summary_basic(self) -> None:
        """generate_qc_summary computes standard QC stats from read lengths."""
        read_lengths = [1000, 2000, 3000, 5000, 8000, 12000]
        summary = generate_qc_summary(read_lengths)

        assert summary["total_reads"] == 6
        assert summary["total_bases"] == 31000
        assert summary["n50"] > 0
        assert summary["mean_length"] > 0
        assert summary["median_length"] > 0
        assert summary["min_length"] == 1000
        assert summary["max_length"] == 12000
        assert summary["std_dev"] > 0
        assert summary["status"] == "complete"

    def test_qc_summary_empty_reads(self) -> None:
        """generate_qc_summary returns status 'no_data' for empty input."""
        summary = generate_qc_summary([])
        assert summary["total_reads"] == 0
        assert summary["total_bases"] == 0
        assert summary["n50"] == 0
        assert summary["status"] == "no_data"

    def test_qc_summary_with_quality_scores(self) -> None:
        """generate_qc_summary incorporates quality score statistics."""
        read_lengths = [5000, 10000, 15000]
        quality_scores = [18.5, 22.0, 25.3]

        summary = generate_qc_summary(read_lengths, quality_scores=quality_scores)

        assert "mean_quality" in summary
        assert summary["mean_quality"] == pytest.approx(21.93, abs=0.01)
        assert "q7_fraction" in summary
        assert summary["q7_fraction"] == pytest.approx(1.0)
        assert "q10_fraction" in summary
        assert summary["q10_fraction"] == pytest.approx(1.0)
        assert "q20_fraction" in summary
        # Two of three scores >= 20
        assert summary["q20_fraction"] == pytest.approx(2 / 3, abs=0.01)

    def test_qc_summary_with_filter_count(self) -> None:
        """generate_qc_summary includes filter statistics when provided."""
        read_lengths = [2000, 5000, 8000]  # 3 reads that passed
        summary = generate_qc_summary(read_lengths, filtered_count=7)

        assert summary["reads_before_filter"] == 10  # 3 + 7
        assert summary["reads_removed"] == 7
        assert summary["filter_pass_rate"] == pytest.approx(0.3)

    def test_qc_summary_with_adapter_and_chimera_counts(self) -> None:
        """generate_qc_summary records adapter trimming and chimera splitting counts."""
        read_lengths = [3000, 4000]
        summary = generate_qc_summary(
            read_lengths,
            adapter_trimmed=5,
            chimeras_split=2,
        )
        assert summary["adapter_trimmed"] == 5
        assert summary["chimeras_split"] == 2

    def test_qc_summary_single_read(self) -> None:
        """generate_qc_summary works with a single read."""
        summary = generate_qc_summary([50000])
        assert summary["total_reads"] == 1
        assert summary["total_bases"] == 50000
        assert summary["mean_length"] == 50000.0
        assert summary["n50"] == 50000

    def test_qc_summary_realistic_ont_distribution(self) -> None:
        """generate_qc_summary with realistic ONT read length distribution."""
        import random

        rng = random.Random(42)
        # Lognormal distribution typical of ONT
        read_lengths = sorted([int(max(200, rng.lognormvariate(8.5, 1.2))) for _ in range(500)])

        summary = generate_qc_summary(read_lengths)
        assert summary["total_reads"] == 500
        assert summary["total_bases"] > 0
        assert summary["n50"] > 0
        assert summary["n90"] > 0
        assert summary["l50"] >= 1
        assert summary["min_length"] >= 200
        assert summary["status"] == "complete"


# ---------------------------------------------------------------------------
# QCReport dataclass
# ---------------------------------------------------------------------------


class TestQCReport:
    """Tests for QCReport dataclass."""

    def test_qc_report_default_values(self) -> None:
        """QCReport has correct defaults."""
        report = QCReport()
        assert report.sample_name == ""
        assert report.total_reads == 0
        assert report.total_bases == 0
        assert report.n50 == 0
        assert report.mean_length == 0.0
        assert report.mean_quality == 0.0
        assert report.reads_passed_filter == 0
        assert report.filter_pass_rate == 0.0
        assert report.adapter_trimmed == 0
        assert report.length_distribution == {}
        assert report.quality_distribution == {}
        assert report.timestamp == ""

    def test_qc_report_custom_values(self) -> None:
        """QCReport can be initialized with custom values."""
        report = QCReport(
            sample_name="test_sample",
            total_reads=1000,
            total_bases=5_000_000,
            n50=8000,
            mean_length=5000.0,
            mean_quality=15.5,
            reads_passed_filter=950,
            filter_pass_rate=0.95,
            adapter_trimmed=20,
            length_distribution={"mean": 5000.0, "median": 4500.0},
            quality_distribution={"mean_quality": 15.5},
            timestamp="2024-01-01T00:00:00Z",
        )
        assert report.sample_name == "test_sample"
        assert report.total_reads == 1000
        assert report.total_bases == 5_000_000
        assert report.n50 == 8000
        assert report.filter_pass_rate == 0.95
        assert report.length_distribution["mean"] == 5000.0

    def test_qc_report_export_json(self, tmp_path: Path) -> None:
        """QCReport can be exported to JSON via export_report."""
        report = QCReport(
            sample_name="export_test",
            total_reads=500,
            total_bases=2_500_000,
            n50=6000,
            mean_length=5000.0,
            mean_quality=18.0,
            timestamp="2024-06-15T12:00:00Z",
        )
        output_file = tmp_path / "qc_report.json"
        result_path = export_report(report, output_file, format="json")

        assert result_path == output_file
        assert output_file.exists()
        assert output_file.stat().st_size > 0

        import json

        with open(output_file) as f:
            data = json.load(f)
        assert data["sample_name"] == "export_test"
        assert data["total_reads"] == 500
        assert data["n50"] == 6000

    def test_qc_report_export_text(self, tmp_path: Path) -> None:
        """QCReport can be exported to text format."""
        report = QCReport(
            sample_name="text_test",
            total_reads=100,
            total_bases=500_000,
            n50=7000,
        )
        output_file = tmp_path / "qc_report.txt"
        result_path = export_report(report, output_file, format="text")

        assert result_path == output_file
        assert output_file.exists()
        content = output_file.read_text()
        assert "text_test" in content or "Text Test" in content

    def test_qc_report_export_html(self, tmp_path: Path) -> None:
        """QCReport can be exported to HTML format."""
        report = QCReport(
            sample_name="html_test",
            total_reads=200,
            total_bases=1_000_000,
            n50=5000,
            mean_quality=20.0,
        )
        output_file = tmp_path / "qc_report.html"
        result_path = export_report(report, output_file, format="html")

        assert result_path == output_file
        assert output_file.exists()
        content = output_file.read_text()
        assert "<html" in content
        assert "html_test" in content

    def test_export_report_invalid_format(self, tmp_path: Path) -> None:
        """export_report raises ValueError for unsupported format."""
        report = QCReport(sample_name="test")
        with pytest.raises(ValueError, match="Unsupported report format"):
            export_report(report, tmp_path / "report.xyz", format="xyz")


# ---------------------------------------------------------------------------
# generate_qc_report from PipelineResult-like object
# ---------------------------------------------------------------------------


class TestGenerateQcReport:
    """Tests for generate_qc_report function."""

    def test_generate_qc_report_from_summary(self) -> None:
        """generate_qc_report extracts metrics from a pipeline result summary."""

        class FakeResult:
            """Minimal PipelineResult-like object for testing."""

            def __init__(self) -> None:
                self.summary = {
                    "total_reads_input": 100,
                    "reads_after_processing": 90,
                }
                self.steps: list[Any] = []
                self.output_dir = Path("output/qc")

        result = FakeResult()
        report = generate_qc_report(result)

        assert isinstance(report, QCReport)
        assert report.total_reads == 100
        assert report.reads_passed_filter == 90
        assert report.timestamp != ""

    def test_generate_qc_report_empty_summary(self) -> None:
        """generate_qc_report handles empty pipeline result gracefully."""

        class EmptyResult:
            """Minimal PipelineResult-like object with no data."""

            def __init__(self) -> None:
                self.summary: dict[str, Any] = {}
                self.steps: list[Any] = []
                self.output_dir = Path(".")

        result = EmptyResult()
        report = generate_qc_report(result)

        assert isinstance(report, QCReport)
        assert report.total_reads == 0
        assert report.n50 == 0


# ---------------------------------------------------------------------------
# generate_run_summary
# ---------------------------------------------------------------------------


class TestGenerateRunSummary:
    """Tests for generate_run_summary function."""

    def test_run_summary_empty(self) -> None:
        """generate_run_summary with no results returns minimal summary."""
        summary = generate_run_summary([])
        assert summary["num_runs"] == 0
        assert summary["runs"] == []
        assert "timestamp" in summary

    def test_run_summary_multiple_runs(self) -> None:
        """generate_run_summary aggregates across multiple pipeline results."""
        from metainformant.longread.workflow.orchestrator import PipelineResult, PipelineStep

        results = [
            PipelineResult(
                pipeline_name="qc",
                steps=[
                    PipelineStep(name="step1", function=lambda ctx: None, status="completed"),
                    PipelineStep(name="step2", function=lambda ctx: None, status="completed"),
                ],
                success=True,
                total_duration=5.0,
                output_dir=Path("output/run1"),
            ),
            PipelineResult(
                pipeline_name="assembly",
                steps=[
                    PipelineStep(name="step1", function=lambda ctx: None, status="completed"),
                    PipelineStep(name="step2", function=lambda ctx: None, status="failed", error="oops"),
                ],
                success=False,
                total_duration=10.0,
                output_dir=Path("output/run2"),
            ),
        ]

        summary = generate_run_summary(results)
        assert summary["num_runs"] == 2
        assert summary["success_count"] == 1
        assert summary["failure_count"] == 1
        assert summary["total_duration"] == pytest.approx(15.0)
        assert summary["mean_duration"] == pytest.approx(7.5)
        assert len(summary["runs"]) == 2

        # Check first run details
        run1 = summary["runs"][0]
        assert run1["pipeline_name"] == "qc"
        assert run1["success"] is True
        assert run1["steps_completed"] == 2
        assert run1["steps_failed"] == 0

        # Check second run details
        run2 = summary["runs"][1]
        assert run2["pipeline_name"] == "assembly"
        assert run2["success"] is False
        assert run2["steps_failed"] == 1
