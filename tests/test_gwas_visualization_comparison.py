"""Tests for GWAS comparison visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas.visualization.visualization_comparison import (
    miami_plot,
    multi_trait_manhattan,
    cross_cohort_forest,
    concordance_plot,
)


def test_miami_plot(tmp_path: Path) -> None:
    """Test Miami plot (back-to-back Manhattan plot) for two traits."""
    # Create sample data for two traits
    trait1_results = [
        {"CHROM": "1", "POS": 1000, "p_value": 1e-8},
        {"CHROM": "1", "POS": 2000, "p_value": 1e-6},
        {"CHROM": "2", "POS": 1500, "p_value": 1e-5},
    ]

    trait2_results = [
        {"CHROM": "1", "POS": 1200, "p_value": 1e-7},
        {"CHROM": "1", "POS": 2200, "p_value": 1e-4},
        {"CHROM": "2", "POS": 1600, "p_value": 1e-6},
    ]

    output_path = tmp_path / "miami_plot.png"
    result = miami_plot(trait1_results, trait2_results, output_path, trait1_name="Trait A", trait2_name="Trait B")

    assert result["status"] == "success"
    assert output_path.exists()
    assert result["trait1_variants"] == 3
    assert result["trait2_variants"] == 3


def test_multi_trait_manhattan(tmp_path: Path) -> None:
    """Test multi-trait Manhattan plot for multiple phenotypes."""
    # Create sample data for multiple traits
    trait_data = {
        "Trait A": [
            {"CHROM": "1", "POS": 1000, "p_value": 1e-8},
            {"CHROM": "1", "POS": 2000, "p_value": 1e-6},
        ],
        "Trait B": [
            {"CHROM": "1", "POS": 1200, "p_value": 1e-7},
            {"CHROM": "1", "POS": 2200, "p_value": 1e-4},
        ],
        "Trait C": [
            {"CHROM": "1", "POS": 1100, "p_value": 1e-5},
            {"CHROM": "1", "POS": 2100, "p_value": 1e-3},
        ],
    }

    output_path = tmp_path / "multi_trait_manhattan.png"
    result = multi_trait_manhattan(trait_data, output_path, title="Multi-Trait GWAS")

    assert result["status"] == "success"
    assert output_path.exists()
    assert result["num_traits"] == 3


def test_concordance_plot(tmp_path: Path) -> None:
    """Test concordance plot between discovery and replication cohorts."""
    # Create sample data for discovery cohort
    discovery_results = [
        {"CHROM": "1", "POS": 1000, "beta": 0.8, "p_value": 1e-8},
        {"CHROM": "1", "POS": 2000, "beta": -0.6, "p_value": 1e-6},
        {"CHROM": "2", "POS": 1500, "beta": 0.4, "p_value": 1e-5},
    ]

    # Create sample data for replication cohort (matching variants)
    replication_results = [
        {"CHROM": "1", "POS": 1000, "beta": 0.7, "p_value": 1e-7},
        {"CHROM": "1", "POS": 2000, "beta": -0.5, "p_value": 1e-5},
        {"CHROM": "2", "POS": 1500, "beta": 0.3, "p_value": 1e-4},
    ]

    output_path = tmp_path / "concordance_plot.png"
    result = concordance_plot(
        discovery_results, replication_results, output_path, title="Discovery vs Replication Concordance"
    )

    assert result["status"] == "success"
    assert output_path.exists()


def test_cross_cohort_forest(tmp_path: Path) -> None:
    """Test cross-cohort forest plot (meta-analysis placeholder)."""
    # Create sample cohort data
    cohort_results = {
        "Cohort A": [
            {"CHROM": "1", "POS": 1000, "beta": 0.8, "p_value": 1e-8},
        ],
        "Cohort B": [
            {"CHROM": "1", "POS": 1000, "beta": 0.7, "p_value": 1e-7},
        ],
    }

    output_path = tmp_path / "forest_plot.png"
    result = cross_cohort_forest(cohort_results, output_path, variant_id="1:1000")

    # Currently returns skipped status as meta-analysis is not fully implemented
    assert result["status"] == "skipped"
    assert "message" in result


def test_miami_plot_empty(tmp_path: Path) -> None:
    """Test Miami plot with empty data."""
    empty_results = []
    sample_results = [{"CHROM": "1", "POS": 1000, "p_value": 1e-8}]

    output_path = tmp_path / "miami_empty.png"
    result = miami_plot(empty_results, sample_results, output_path)

    # Should handle gracefully
    assert result["status"] == "failed"
    assert "error" in result


def test_multi_trait_manhattan_single_trait(tmp_path: Path) -> None:
    """Test multi-trait Manhattan with single trait."""
    trait_data = {
        "Single Trait": [
            {"CHROM": "1", "POS": 1000, "p_value": 1e-8},
            {"CHROM": "1", "POS": 2000, "p_value": 1e-6},
        ]
    }

    output_path = tmp_path / "single_trait.png"
    result = multi_trait_manhattan(trait_data, output_path)

    # Should fail with < 2 traits
    assert result["status"] == "failed"
    assert "error" in result


def test_concordance_plot_perfect_correlation(tmp_path: Path) -> None:
    """Test concordance plot with perfect correlation."""
    # Create perfectly correlated data
    discovery_results = [
        {"CHROM": "1", "POS": 1000, "beta": 0.5, "p_value": 1e-8},
        {"CHROM": "1", "POS": 2000, "beta": -0.3, "p_value": 1e-6},
    ]

    replication_results = [
        {"CHROM": "1", "POS": 1000, "beta": 0.5, "p_value": 1e-7},  # Perfect match
        {"CHROM": "1", "POS": 2000, "beta": -0.3, "p_value": 1e-5},  # Perfect match
    ]

    output_path = tmp_path / "perfect_concordance.png"
    result = concordance_plot(discovery_results, replication_results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()
    assert result["correlation"] > 0.99  # Should be nearly perfect correlation
