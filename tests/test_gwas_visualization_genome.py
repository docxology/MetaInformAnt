"""Tests for GWAS genome-wide visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas.visualization.visualization_genome import (
    chromosome_ideogram,
    circular_manhattan_plot,
    genome_wide_ld_heatmap,
)
from metainformant.gwas.visualization.visualization_genome import manhattan_plot as manhattan_plot_genome


def test_manhattan_plot_genome_basic(tmp_path: Path) -> None:
    """Test genome-wide Manhattan plot generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs2", "p_value": 0.01},
        {"CHROM": "chr2", "POS": 1000, "ID": "rs3", "p_value": 1e-6},
        {"CHROM": "chr2", "POS": 2000, "ID": "rs4", "p_value": 0.05},
    ]

    output_path = tmp_path / "manhattan_genome.png"
    result = manhattan_plot_genome(results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_circular_manhattan_plot(tmp_path: Path) -> None:
    """Test circular Manhattan plot generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs2", "p_value": 0.01},
        {"CHROM": "chr2", "POS": 1000, "ID": "rs3", "p_value": 1e-6},
    ]

    output_path = tmp_path / "circular_manhattan.png"
    result = circular_manhattan_plot(results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()


def test_chromosome_ideogram(tmp_path: Path) -> None:
    """Test chromosome ideogram generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8},
        {"CHROM": "chr2", "POS": 2000, "ID": "rs2", "p_value": 1e-6},
    ]

    output_path = tmp_path / "ideogram.png"
    result = chromosome_ideogram(results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()


def test_genome_wide_ld_heatmap(tmp_path: Path) -> None:
    """Test genome-wide LD heatmap generation."""
    # Create mock LD data
    ld_data = [
        {"CHROM": "chr1", "POS1": 1000, "POS2": 2000, "r2": 0.5},
        {"CHROM": "chr1", "POS1": 1000, "POS2": 3000, "r2": 0.3},
        {"CHROM": "chr2", "POS1": 1000, "POS2": 2000, "r2": 0.4},
    ]

    output_path = tmp_path / "ld_heatmap.png"
    result = genome_wide_ld_heatmap(ld_data, output_path)

    assert result["status"] == "success"
    assert output_path.exists()


def test_manhattan_plot_genome_empty_results(tmp_path: Path) -> None:
    """Test Manhattan plot with empty results."""
    output_path = tmp_path / "empty.png"
    result = manhattan_plot_genome([], output_path)

    assert result["status"] == "failed"
    assert "error" in result


def test_manhattan_plot_genome_with_thresholds(tmp_path: Path) -> None:
    """Test Manhattan plot with custom significance thresholds."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-9},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs2", "p_value": 1e-6},
        {"CHROM": "chr2", "POS": 1000, "ID": "rs3", "p_value": 0.01},
    ]

    output_path = tmp_path / "manhattan_thresholds.png"
    result = manhattan_plot_genome(
        results,
        output_path,
        significance_threshold=1e-8,
        suggestiveness_threshold=1e-5,
    )

    assert result["status"] == "success"
    assert output_path.exists()
