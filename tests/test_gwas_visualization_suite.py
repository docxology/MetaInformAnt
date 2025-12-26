"""Tests for GWAS visualization suite functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas.visualization_suite import (
    generate_all_plots,
)


def test_generate_all_plots(tmp_path: Path) -> None:
    """Test generating all plots for a GWAS analysis."""
    # Create sample GWAS results file
    association_file = tmp_path / "gwas_results.tsv"
    with open(association_file, 'w') as f:
        f.write("CHROM\tPOS\tp_value\tBETA\tSE\n")
        f.write("1\t1000\t1e-8\t0.2\t0.05\n")
        f.write("1\t2000\t1e-6\t-0.1\t0.04\n")
        f.write("2\t1500\t1e-5\t0.15\t0.06\n")

    output_dir = tmp_path / "gwas_plots"
    output_dir.mkdir()

    result = generate_all_plots(
        association_results=association_file,
        output_dir=output_dir,
        significance_threshold=1e-6,
    )

    assert result["status"] == "completed"
    assert result["num_plots_generated"] >= 5  # Should generate several plots

    # Check that plot files were created
    plot_files = list(output_dir.glob("*.png"))
    assert len(plot_files) >= 5


def test_generate_all_plots_minimal(tmp_path: Path) -> None:
    """Test generating plots with minimal data."""
    # Create minimal GWAS results file
    association_file = tmp_path / "minimal_results.tsv"
    with open(association_file, 'w') as f:
        f.write("CHROM\tPOS\tp_value\n")
        f.write("1\t1000\t1e-7\n")
        f.write("2\t2000\t1e-6\n")

    output_dir = tmp_path / "minimal_plots"
    output_dir.mkdir()

    result = generate_all_plots(
        association_results=association_file,
        output_dir=output_dir,
    )

    assert result["status"] == "completed"
    assert result["num_plots_generated"] >= 3  # Should generate basic plots


def test_generate_all_plots_empty_results(tmp_path: Path) -> None:
    """Test generating plots with empty results."""
    association_results = []

    output_dir = tmp_path / "empty_plots"

    result = generate_all_plots(
        association_results=association_results,
        output_dir=output_dir,
    )

    # Should handle gracefully
    assert isinstance(result, dict)
    assert "status" in result
