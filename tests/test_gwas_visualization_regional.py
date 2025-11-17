"""Tests for GWAS regional visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas.visualization_regional import (
    gene_annotation_plot,
    recombination_rate_plot,
    regional_plot as regional_plot_detailed,
    regional_ld_plot,
)


def test_regional_plot_detailed(tmp_path: Path) -> None:
    """Test detailed regional plot generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8, "beta": 0.5},
        {"CHROM": "chr1", "POS": 1500, "ID": "rs2", "p_value": 1e-6, "beta": 0.3},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs3", "p_value": 0.01, "beta": 0.1},
    ]

    output_path = tmp_path / "regional_detailed.png"
    result = regional_plot_detailed(
        results, "chr1", 500, 2500, output_path
    )

    assert result["status"] == "success"
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_regional_ld_plot(tmp_path: Path) -> None:
    """Test regional LD plot generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8},
        {"CHROM": "chr1", "POS": 1500, "ID": "rs2", "p_value": 1e-6},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs3", "p_value": 0.01},
    ]

    # Mock LD matrix
    ld_matrix = np.array([[1.0, 0.8, 0.5], [0.8, 1.0, 0.6], [0.5, 0.6, 1.0]])

    output_path = tmp_path / "regional_ld.png"
    result = regional_ld_plot(
        results, "chr1", 500, 2500, ld_matrix, output_path
    )

    assert result["status"] == "success"
    assert output_path.exists()


def test_gene_annotation_plot(tmp_path: Path) -> None:
    """Test gene annotation plot generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs2", "p_value": 1e-6},
    ]

    # Mock gene annotations
    genes = [
        {"name": "GENE1", "start": 800, "end": 1200, "strand": "+"},
        {"name": "GENE2", "start": 1800, "end": 2200, "strand": "-"},
    ]

    output_path = tmp_path / "gene_annotation.png"
    result = gene_annotation_plot(
        results, "chr1", 500, 2500, genes, output_path
    )

    assert result["status"] == "success"
    assert output_path.exists()


def test_recombination_rate_plot(tmp_path: Path) -> None:
    """Test recombination rate plot generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs2", "p_value": 1e-6},
    ]

    # Mock recombination rates
    recombination_data = [
        {"POS": 500, "rate": 1.0},
        {"POS": 1500, "rate": 1.5},
        {"POS": 2500, "rate": 2.0},
    ]

    output_path = tmp_path / "recombination_rate.png"
    result = recombination_rate_plot(
        results, "chr1", 500, 2500, recombination_data, output_path
    )

    assert result["status"] == "success"
    assert output_path.exists()


def test_regional_plot_detailed_empty_region(tmp_path: Path) -> None:
    """Test regional plot with empty region."""
    results = [
        {"CHROM": "chr2", "POS": 1000, "ID": "rs1", "p_value": 1e-8},
    ]

    output_path = tmp_path / "empty_region.png"
    result = regional_plot_detailed(
        results, "chr1", 500, 2500, output_path
    )

    # Should handle empty region gracefully
    assert "status" in result

