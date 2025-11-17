"""Tests for GWAS statistical visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas.visualization_statistical import (
    lambda_gc_plot,
    power_plot,
    qq_plot as qq_plot_statistical,
    qq_plot_stratified,
    volcano_plot,
)


def test_qq_plot_statistical_basic(tmp_path: Path) -> None:
    """Test statistical Q-Q plot generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs2", "p_value": 0.01},
        {"CHROM": "chr2", "POS": 1000, "ID": "rs3", "p_value": 1e-6},
        {"CHROM": "chr2", "POS": 2000, "ID": "rs4", "p_value": 0.05},
    ]

    output_path = tmp_path / "qq_statistical.png"
    result = qq_plot_statistical(results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_qq_plot_stratified(tmp_path: Path) -> None:
    """Test stratified Q-Q plot generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8, "MAF": 0.05},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs2", "p_value": 0.01, "MAF": 0.1},
        {"CHROM": "chr2", "POS": 1000, "ID": "rs3", "p_value": 1e-6, "MAF": 0.2},
    ]

    output_path = tmp_path / "qq_stratified.png"
    result = qq_plot_stratified(results, output_path, stratify_by="MAF")

    assert result["status"] == "success"
    assert output_path.exists()


def test_lambda_gc_plot(tmp_path: Path) -> None:
    """Test lambda GC plot generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs2", "p_value": 0.01},
        {"CHROM": "chr2", "POS": 1000, "ID": "rs3", "p_value": 1e-6},
    ]

    output_path = tmp_path / "lambda_gc.png"
    result = lambda_gc_plot(results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()
    assert "lambda_gc" in result


def test_volcano_plot(tmp_path: Path) -> None:
    """Test volcano plot generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8, "beta": 0.5},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs2", "p_value": 0.01, "beta": 0.2},
        {"CHROM": "chr2", "POS": 1000, "ID": "rs3", "p_value": 1e-6, "beta": -0.3},
    ]

    output_path = tmp_path / "volcano.png"
    result = volcano_plot(results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()


def test_power_plot(tmp_path: Path) -> None:
    """Test power plot generation."""
    output_path = tmp_path / "power.png"
    result = power_plot(
        output_path,
        sample_sizes=[100, 500, 1000, 5000],
        effect_sizes=[0.1, 0.2, 0.5],
        alpha=5e-8,
    )

    assert result["status"] == "success"
    assert output_path.exists()


def test_qq_plot_statistical_empty(tmp_path: Path) -> None:
    """Test Q-Q plot with empty results."""
    output_path = tmp_path / "empty_qq.png"
    result = qq_plot_statistical([], output_path)

    assert result["status"] == "failed"
    assert "error" in result

