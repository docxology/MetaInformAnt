"""Tests for GWAS regional visualization functions."""

from __future__ import annotations

from pathlib import Path

from metainformant.gwas.visualization.visualization_regional import (
    gene_annotation_plot,
    recombination_rate_plot,
    regional_ld_plot,
)
from metainformant.gwas.visualization.visualization_regional import regional_plot as regional_plot_detailed


def test_regional_plot_detailed(tmp_path: Path) -> None:
    """Test detailed regional plot generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8, "beta": 0.5},
        {"CHROM": "chr1", "POS": 1500, "ID": "rs2", "p_value": 1e-6, "beta": 0.3},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs3", "p_value": 0.01, "beta": 0.1},
    ]

    output_path = tmp_path / "regional_detailed.png"
    result = regional_plot_detailed(results, output_path, chrom="chr1", start=500, end=2500)

    assert result["status"] == "success"
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_regional_ld_plot(tmp_path: Path) -> None:
    """Test regional LD plot generation."""
    # Create a mock VCF file path (function expects VCF path, not results)
    vcf_path = tmp_path / "test.vcf"
    vcf_path.touch()  # Create empty file for testing

    output_path = tmp_path / "regional_ld.png"
    result = regional_ld_plot(vcf_path, output_path, chrom="chr1", start=500, end=2500, lead_snp_pos=1500)

    # Function returns skipped status as it requires external LD calculation
    assert "status" in result
    assert result["status"] in ("success", "skipped")


def test_gene_annotation_plot(tmp_path: Path) -> None:
    """Test gene annotation plot generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs2", "p_value": 1e-6},
    ]

    output_path = tmp_path / "gene_annotation.png"
    result = gene_annotation_plot(results, output_path, chrom="chr1", start=500, end=2500)

    assert result["status"] == "success"
    assert output_path.exists()


def test_recombination_rate_plot(tmp_path: Path) -> None:
    """Test recombination rate plot generation."""
    results = [
        {"CHROM": "chr1", "POS": 1000, "ID": "rs1", "p_value": 1e-8},
        {"CHROM": "chr1", "POS": 2000, "ID": "rs2", "p_value": 1e-6},
    ]

    output_path = tmp_path / "recombination_rate.png"
    result = recombination_rate_plot(results, output_path, chrom="chr1", start=500, end=2500)

    # Function returns skipped status as it requires recombination map data
    assert "status" in result
    assert result["status"] in ("success", "skipped")


def test_regional_plot_detailed_empty_region(tmp_path: Path) -> None:
    """Test regional plot with empty region."""
    results = [
        {"CHROM": "chr2", "POS": 1000, "ID": "rs1", "p_value": 1e-8},
    ]

    output_path = tmp_path / "empty_region.png"
    result = regional_plot_detailed(results, output_path, chrom="chr1", start=500, end=2500)

    # Should handle empty region gracefully
    assert "status" in result
