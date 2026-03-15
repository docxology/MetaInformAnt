"""Tests for pharmacogenomics visualization: all plot functions.

Tests that each plot function creates valid PNG files using realistic
pharmacogenomic data.

NO MOCKING -- all tests use real implementations.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest

from metainformant.pharmacogenomics.visualization.plots import (
    plot_acmg_criteria,
    plot_activity_score_distribution,
    plot_allele_frequencies,
    plot_drug_response_heatmap,
    plot_metabolizer_status,
    plot_population_comparison,
)


class TestPlotMetabolizerStatus:
    """Tests for metabolizer status distribution plot."""

    def test_creates_png(self, tmp_path: Path) -> None:
        phenotypes = {"PM": 5, "IM": 20, "NM": 65, "RM": 8, "UM": 2}
        output = tmp_path / "metabolizer_status.png"
        result = plot_metabolizer_status(phenotypes, output)
        assert result == output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_with_full_names(self, tmp_path: Path) -> None:
        phenotypes = {
            "Poor Metabolizer": 3,
            "Intermediate Metabolizer": 15,
            "Normal Metabolizer": 70,
            "Ultrarapid Metabolizer": 12,
        }
        output = tmp_path / "metabolizer_full.png"
        result = plot_metabolizer_status(phenotypes, output, title="CYP2D6 Status")
        assert output.exists()

    def test_with_float_values(self, tmp_path: Path) -> None:
        phenotypes = {"PM": 0.05, "IM": 0.20, "NM": 0.65, "RM": 0.08, "UM": 0.02}
        output = tmp_path / "metabolizer_freq.png"
        result = plot_metabolizer_status(phenotypes, output)
        assert output.exists()


class TestPlotAlleleFrequencies:
    """Tests for allele frequency distribution plot."""

    def test_creates_bar_chart(self, tmp_path: Path) -> None:
        freqs = {"*1": 0.45, "*2": 0.25, "*4": 0.15, "*10": 0.10, "*17": 0.05}
        output = tmp_path / "allele_freq_bar.png"
        result = plot_allele_frequencies(freqs, "CYP2D6", output, plot_type="bar")
        assert result == output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_creates_pie_chart(self, tmp_path: Path) -> None:
        freqs = {"*1": 0.60, "*2": 0.15, "*3": 0.10, "*17": 0.15}
        output = tmp_path / "allele_freq_pie.png"
        result = plot_allele_frequencies(freqs, "CYP2C19", output, plot_type="pie")
        assert output.exists()

    def test_custom_title(self, tmp_path: Path) -> None:
        freqs = {"*1": 0.80, "*3": 0.20}
        output = tmp_path / "allele_freq_custom.png"
        result = plot_allele_frequencies(freqs, "CYP2C9", output, title="Custom Allele Plot")
        assert output.exists()


class TestPlotActivityScoreDistribution:
    """Tests for activity score distribution histogram."""

    def test_creates_histogram(self, tmp_path: Path) -> None:
        scores = [0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0, 1.5, 2.0, 2.0, 2.0, 2.0, 2.5]
        output = tmp_path / "activity_dist.png"
        result = plot_activity_score_distribution(scores, "CYP2D6", output)
        assert result == output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_single_value(self, tmp_path: Path) -> None:
        scores = [2.0]
        output = tmp_path / "activity_single.png"
        result = plot_activity_score_distribution(scores, "CYP2D6", output)
        assert output.exists()

    def test_cyp2c19_scores(self, tmp_path: Path) -> None:
        scores = [0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.5, 2.5, 3.0]
        output = tmp_path / "activity_cyp2c19.png"
        result = plot_activity_score_distribution(scores, "CYP2C19", output)
        assert output.exists()


class TestPlotDrugResponseHeatmap:
    """Tests for drug-gene response heatmap."""

    def test_creates_heatmap(self, tmp_path: Path) -> None:
        drugs = ["codeine", "clopidogrel", "warfarin", "simvastatin"]
        genes = ["CYP2D6", "CYP2C19", "CYP2C9", "SLCO1B1"]
        phenotypes = {
            ("codeine", "CYP2D6"): "Major",
            ("clopidogrel", "CYP2C19"): "Major",
            ("warfarin", "CYP2C9"): "Moderate",
            ("simvastatin", "SLCO1B1"): "Major",
        }
        output = tmp_path / "drug_heatmap.png"
        result = plot_drug_response_heatmap(drugs, genes, phenotypes, output)
        assert result == output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_with_missing_combinations(self, tmp_path: Path) -> None:
        drugs = ["codeine", "warfarin"]
        genes = ["CYP2D6", "CYP2C9"]
        phenotypes = {("codeine", "CYP2D6"): "Major"}
        output = tmp_path / "heatmap_sparse.png"
        result = plot_drug_response_heatmap(drugs, genes, phenotypes, output)
        assert output.exists()


class TestPlotPopulationComparison:
    """Tests for cross-population allele frequency comparison."""

    def test_creates_grouped_bar(self, tmp_path: Path) -> None:
        allele_freqs = {
            "European": {"*1": 0.60, "*2": 0.20, "*4": 0.15, "*10": 0.05},
            "East Asian": {"*1": 0.30, "*2": 0.10, "*4": 0.01, "*10": 0.50},
            "African": {"*1": 0.50, "*2": 0.25, "*4": 0.07, "*10": 0.08},
        }
        output = tmp_path / "pop_comparison.png"
        result = plot_population_comparison(allele_freqs, "CYP2D6", output)
        assert result == output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_two_populations(self, tmp_path: Path) -> None:
        allele_freqs = {
            "European": {"*1": 0.85, "*2": 0.15},
            "East Asian": {"*1": 0.70, "*3": 0.07, "*2": 0.23},
        }
        output = tmp_path / "pop_two.png"
        result = plot_population_comparison(allele_freqs, "CYP2C19", output)
        assert output.exists()


class TestPlotACMGCriteria:
    """Tests for ACMG criteria evaluation plot."""

    def test_creates_criteria_plot(self, tmp_path: Path) -> None:
        criteria = {
            "PVS1": True,
            "PS1": False,
            "PS3": True,
            "PM1": False,
            "PM2": True,
            "PP3": True,
            "PP5": False,
            "BA1": False,
            "BS1": False,
            "BP4": False,
            "BP7": False,
        }
        output = tmp_path / "acmg_criteria.png"
        result = plot_acmg_criteria(criteria, output)
        assert result == output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_all_met(self, tmp_path: Path) -> None:
        criteria = {"PVS1": True, "PS1": True, "PM2": True, "PP3": True}
        output = tmp_path / "acmg_all_met.png"
        result = plot_acmg_criteria(criteria, output)
        assert output.exists()

    def test_none_met(self, tmp_path: Path) -> None:
        criteria = {"PVS1": False, "PS1": False, "BA1": False, "BP4": False}
        output = tmp_path / "acmg_none.png"
        result = plot_acmg_criteria(criteria, output)
        assert output.exists()

    def test_empty_criteria(self, tmp_path: Path) -> None:
        criteria: dict[str, bool] = {}
        output = tmp_path / "acmg_empty.png"
        result = plot_acmg_criteria(criteria, output)
        assert output.exists()
