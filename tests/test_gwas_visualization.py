"""Tests for GWAS visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas import manhattan_plot, qq_plot, regional_plot
from metainformant.core.io import write_tsv


def test_manhattan_plot_basic(tmp_path: Path) -> None:
    """Test generating Manhattan plot from results."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "REF": "A", "ALT": "G", "p_value": 1e-8},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs2", "REF": "T", "ALT": "C", "p_value": 0.01},
        {"CHROM": "chr2", "POS": 1000, "ID": "rs3", "REF": "G", "ALT": "A", "p_value": 1e-6},
        {"CHROM": "chr2", "POS": 2000, "ID": "rs4", "REF": "C", "ALT": "T", "p_value": 0.05},
    ]

    output_path = tmp_path / "manhattan.png"

    result = manhattan_plot(results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_manhattan_plot_from_file(tmp_path: Path) -> None:
    """Test generating Manhattan plot from results file."""
    results_file = tmp_path / "results.tsv"
    results_data = [
        ["CHROM", "POS", "ID", "REF", "ALT", "p_value"],
        ["chr1", "1000", "rs1", "A", "G", "1e-8"],
        ["chr1", "2000", "rs2", "T", "C", "0.01"],
        ["chr2", "1000", "rs3", "G", "A", "1e-6"],
    ]
    write_tsv(results_data, results_file)

    output_path = tmp_path / "manhattan.png"
    result = manhattan_plot(results_file, output_path)

    assert result["status"] == "success"
    assert output_path.exists()


def test_qq_plot_basic(tmp_path: Path) -> None:
    """Test generating Q-Q plot from p-values."""
    # Generate p-values (mix of null and significant)
    np.random.seed(42)
    null_pvalues = np.random.uniform(0.01, 1.0, 1000).tolist()
    significant_pvalues = np.random.uniform(1e-8, 1e-5, 10).tolist()
    pvalues = null_pvalues + significant_pvalues

    output_path = tmp_path / "qq.png"

    result = qq_plot(pvalues, output_path)

    assert result["status"] == "success"
    assert output_path.exists()
    assert "lambda_gc" in result
    assert result["lambda_gc"] > 0


def test_qq_plot_from_file(tmp_path: Path) -> None:
    """Test generating Q-Q plot from results file."""
    results_file = tmp_path / "results.tsv"
    results_data = [
        ["CHROM", "POS", "p_value"],
        ["chr1", "1000", "0.001"],
        ["chr1", "2000", "0.01"],
        ["chr1", "3000", "0.05"],
        ["chr1", "4000", "0.1"],
        ["chr1", "5000", "0.5"],
    ]
    write_tsv(results_data, results_file)

    output_path = tmp_path / "qq.png"
    result = qq_plot(results_file, output_path)

    assert result["status"] == "success"
    assert output_path.exists()


def test_regional_plot_basic(tmp_path: Path) -> None:
    """Test generating regional association plot."""
    results = [
        {"CHROM": "chr1", "POS": 100000, "ID": "rs1", "REF": "A", "ALT": "G", "p_value": 1e-6},
        {"CHROM": "chr1", "POS": 100100, "ID": "rs2", "REF": "T", "ALT": "C", "p_value": 1e-5},
        {"CHROM": "chr1", "POS": 100200, "ID": "rs3", "REF": "G", "ALT": "A", "p_value": 0.01},
        {"CHROM": "chr1", "POS": 100300, "ID": "rs4", "REF": "C", "ALT": "T", "p_value": 0.05},
    ]

    output_path = tmp_path / "regional.png"
    region = "chr1:100000-100400"

    result = regional_plot(results, region, output_path)

    assert result["status"] == "success"
    assert output_path.exists()
    assert result["region"] == region


def test_regional_plot_no_variants(tmp_path: Path) -> None:
    """Test regional plot with no variants in region."""
    results = [
        {"CHROM": "chr1", "POS": 100000, "ID": "rs1", "REF": "A", "ALT": "G", "p_value": 0.01},
    ]

    output_path = tmp_path / "regional.png"
    region = "chr2:200000-200100"  # Different chromosome

    result = regional_plot(results, region, output_path)

    assert result["status"] == "failed"
    assert "error" in result


def test_visualization_empty_results(tmp_path: Path) -> None:
    """Test visualization functions with empty results."""
    output_path = tmp_path / "empty.png"

    result = manhattan_plot([], output_path)
    assert result["status"] == "failed"

    result = qq_plot([], output_path)
    assert result["status"] == "failed"

