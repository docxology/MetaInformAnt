"""Tests for longread methylation summary functions not covered by test_longread.py.

Tests generate_methylation_summary from the utils/summary module.
All tests use real implementations -- NO MOCKING.
"""

from __future__ import annotations

import pytest

from metainformant.longread.utils.summary import generate_methylation_summary


class TestGenerateMethylationSummary:
    """Tests for generate_methylation_summary function."""

    def test_methylation_summary_basic(self) -> None:
        """generate_methylation_summary computes correct basic statistics."""
        summary = generate_methylation_summary(
            total_sites=10000,
            methylated_sites=7000,
        )
        assert summary["modification_type"] == "5mC"
        assert summary["total_sites"] == 10000
        assert summary["methylated_sites"] == 7000
        assert summary["unmethylated_sites"] == 3000
        assert summary["global_methylation_rate"] == pytest.approx(0.7)
        assert summary["status"] == "complete"

    def test_methylation_summary_6ma(self) -> None:
        """generate_methylation_summary handles 6mA modification type."""
        summary = generate_methylation_summary(
            total_sites=5000,
            methylated_sites=500,
            modification_type="6mA",
        )
        assert summary["modification_type"] == "6mA"
        assert summary["global_methylation_rate"] == pytest.approx(0.1)

    def test_methylation_summary_zero_sites(self) -> None:
        """generate_methylation_summary handles zero total sites gracefully."""
        summary = generate_methylation_summary(
            total_sites=0,
            methylated_sites=0,
        )
        assert summary["total_sites"] == 0
        assert summary["methylated_sites"] == 0
        assert summary["unmethylated_sites"] == 0
        assert summary["global_methylation_rate"] == 0.0

    def test_methylation_summary_all_methylated(self) -> None:
        """generate_methylation_summary with 100% methylation."""
        summary = generate_methylation_summary(
            total_sites=8000,
            methylated_sites=8000,
        )
        assert summary["global_methylation_rate"] == pytest.approx(1.0)
        assert summary["unmethylated_sites"] == 0

    def test_methylation_summary_no_methylation(self) -> None:
        """generate_methylation_summary with 0% methylation."""
        summary = generate_methylation_summary(
            total_sites=5000,
            methylated_sites=0,
        )
        assert summary["global_methylation_rate"] == pytest.approx(0.0)
        assert summary["unmethylated_sites"] == 5000

    def test_methylation_summary_with_all_parameters(self) -> None:
        """generate_methylation_summary includes all optional parameters."""
        summary = generate_methylation_summary(
            total_sites=50000,
            methylated_sites=35000,
            modification_type="5mC",
            mean_coverage=25.5,
            regions_analyzed=120,
            differential_sites=45,
        )
        assert summary["total_sites"] == 50000
        assert summary["methylated_sites"] == 35000
        assert summary["mean_coverage"] == 25.5
        assert summary["regions_analyzed"] == 120
        assert summary["differential_sites"] == 45
        assert summary["global_methylation_rate"] == pytest.approx(0.7)

    def test_methylation_summary_coverage_rounding(self) -> None:
        """generate_methylation_summary rounds mean_coverage to 1 decimal."""
        summary = generate_methylation_summary(
            total_sites=1000,
            methylated_sites=500,
            mean_coverage=18.6789,
        )
        assert summary["mean_coverage"] == 18.7

    def test_methylation_summary_rate_precision(self) -> None:
        """generate_methylation_summary rounds rate to 4 decimal places."""
        summary = generate_methylation_summary(
            total_sites=3,
            methylated_sites=1,
        )
        # 1/3 = 0.3333...
        assert summary["global_methylation_rate"] == pytest.approx(0.3333, abs=0.0001)

    def test_methylation_summary_realistic_cpg_islands(self) -> None:
        """generate_methylation_summary with realistic CpG island data."""
        # CpG islands typically show ~70-80% methylation in somatic tissue
        summary = generate_methylation_summary(
            total_sites=28000,
            methylated_sites=21000,
            modification_type="5mC",
            mean_coverage=30.0,
            regions_analyzed=500,
            differential_sites=12,
        )
        assert summary["global_methylation_rate"] == pytest.approx(0.75)
        assert summary["mean_coverage"] == 30.0
        assert summary["regions_analyzed"] == 500
        assert summary["status"] == "complete"
