"""Tests for GWAS effect size visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas.visualization.statistical.effects import (
    allelic_series_plot,
    effect_direction_plot,
    effect_size_forest_plot,
    functional_enrichment_plot,
)


def test_effect_size_forest_plot(tmp_path: Path) -> None:
    """Test effect size forest plot."""
    # Create sample GWAS results with effect sizes and significant p-values
    results = [
        {"CHROM": "1", "POS": 1000, "BETA": 0.1, "SE": 0.05, "p_value": 1e-10},
        {"CHROM": "1", "POS": 2000, "BETA": -0.2, "SE": 0.03, "p_value": 1e-12},
        {"CHROM": "2", "POS": 1500, "BETA": 0.05, "SE": 0.04, "p_value": 1e-9},
        {"CHROM": "2", "POS": 2500, "BETA": 0.15, "SE": 0.06, "p_value": 1e-8},
    ]

    output_path = tmp_path / "effect_forest.png"
    result = effect_size_forest_plot(results, output_path, top_n=4)

    assert result["status"] == "success"
    assert output_path.exists()


def test_effect_direction_plot(tmp_path: Path) -> None:
    """Test effect direction plot."""
    # Create sample GWAS results with mixed effect directions
    results = [
        {"CHROM": "1", "POS": 1000, "BETA": 0.1, "p_value": 1e-6},
        {"CHROM": "1", "POS": 2000, "BETA": -0.2, "p_value": 1e-7},
        {"CHROM": "2", "POS": 1500, "BETA": 0.05, "p_value": 1e-5},
        {"CHROM": "2", "POS": 2500, "BETA": -0.15, "p_value": 0.01},
    ]

    output_path = tmp_path / "effect_direction.png"
    result = effect_direction_plot(results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()


def test_functional_enrichment_plot(tmp_path: Path) -> None:
    """Test functional enrichment plot with significant variants."""
    # Create sample GWAS results with significant variants
    results = [
        {"CHROM": "1", "POS": 1000, "p_value": 1e-10},
        {"CHROM": "1", "POS": 2000, "p_value": 1e-12},
        {"CHROM": "2", "POS": 5000, "p_value": 1e-8},
        {"CHROM": "2", "POS": 6000, "p_value": 0.01},  # Not significant
    ]

    output_path = tmp_path / "functional_enrichment.png"
    result = functional_enrichment_plot(results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()
    assert result["total_significant_variants"] == 3  # 3 significant variants
    assert result["categories_found"] > 0
    assert "intergenic" in result["category_counts"]  # Default category


def test_functional_enrichment_plot_no_significant(tmp_path: Path) -> None:
    """Test functional enrichment plot with no significant variants."""
    results = [
        {"CHROM": "1", "POS": 1000, "p_value": 0.1},  # Not significant
        {"CHROM": "2", "POS": 2000, "p_value": 0.05},  # Not significant
    ]

    output_path = tmp_path / "functional_enrichment.png"
    result = functional_enrichment_plot(results, output_path)

    assert result["status"] == "failed"
    assert "No significant variants found" in result["error"]


def test_functional_enrichment_plot_empty_results(tmp_path: Path) -> None:
    """Test functional enrichment plot with empty results."""
    output_path = tmp_path / "functional_enrichment.png"
    result = functional_enrichment_plot([], output_path)

    assert result["status"] == "failed"
    assert "No results provided" in result["error"]


def test_functional_enrichment_plot_invalid_threshold(tmp_path: Path) -> None:
    """Test functional enrichment plot with invalid significance threshold."""
    results = [{"CHROM": "1", "POS": 1000, "p_value": 1e-10}]

    output_path = tmp_path / "functional_enrichment.png"

    # Test negative threshold
    with pytest.raises(ValueError, match="significance_threshold must be"):
        functional_enrichment_plot(results, output_path, significance_threshold=-0.1)

    # Test threshold > 1
    with pytest.raises(ValueError, match="significance_threshold must be"):
        functional_enrichment_plot(results, output_path, significance_threshold=1.5)


def test_allelic_series_plot(tmp_path: Path) -> None:
    """Test allelic series plot (placeholder implementation)."""
    # Create sample GWAS results
    results = [
        {"CHROM": "1", "POS": 1000, "p_value": 1e-8},
        {"CHROM": "1", "POS": 1001, "p_value": 1e-6},  # Multiple alleles at locus
    ]

    output_path = tmp_path / "allelic_series.png"
    result = allelic_series_plot(results, output_path, gene_region=("1", 995, 1005))

    # Currently returns skipped status as it's a placeholder
    assert result["status"] == "skipped"
    assert "message" in result


def test_effect_direction_plot_empty(tmp_path: Path) -> None:
    """Test effect direction plot with empty data."""
    results = []

    output_path = tmp_path / "direction_empty.png"
    result = effect_direction_plot(results, output_path)

    assert result["status"] == "failed"
    assert "error" in result
