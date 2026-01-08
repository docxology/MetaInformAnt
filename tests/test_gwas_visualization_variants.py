"""Tests for GWAS variant-specific visualization functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas.visualization.visualization_variants import (
    maf_distribution,
    variant_density_plot,
    hwe_deviation_plot,
    missingness_plot,
    transition_transversion_plot,
)


def test_maf_distribution(tmp_path: Path) -> None:
    """Test MAF distribution plot."""
    # Create sample GWAS results with MAF data
    results = [
        {"MAF": maf} for maf in np.random.beta(2, 5, 100)
    ]

    output_path = tmp_path / "maf_dist.png"
    result = maf_distribution(results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()


def test_variant_density_plot(tmp_path: Path) -> None:
    """Test variant density plot."""
    # Create sample variant position data
    variants = [
        {"CHROM": "1", "POS": 1000},
        {"CHROM": "1", "POS": 2000},
        {"CHROM": "1", "POS": 1500},
        {"CHROM": "2", "POS": 3000},
    ]

    output_path = tmp_path / "variant_density.png"
    result = variant_density_plot(variants, output_path, window_size=1000)

    assert result["status"] == "success"
    assert output_path.exists()

    assert result["status"] == "success"
    assert output_path.exists()


def test_hwe_deviation_plot(tmp_path: Path) -> None:
    """Test HWE deviation plot."""
    # Create sample GWAS results with HWE p-values
    results = [
        {"HWE_P": hwe_p} for hwe_p in np.random.uniform(0, 1, 100)
    ]

    output_path = tmp_path / "hwe_deviation.png"
    result = hwe_deviation_plot(results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()


def test_missingness_plot_file_not_found(tmp_path: Path) -> None:
    """Test missingness plot with non-existent VCF file."""
    output_path = tmp_path / "missingness.png"
    result = missingness_plot("nonexistent.vcf", output_path)

    assert result["status"] == "failed"
    assert "VCF parsing failed" in result["error"]


def test_missingness_plot_invalid_params(tmp_path: Path) -> None:
    """Test missingness plot with invalid parameters."""
    output_path = tmp_path / "missingness.png"

    # Test with invalid by_sample parameter type
    with pytest.raises(ValueError, match="by_sample must be of type"):
        missingness_plot("dummy.vcf", output_path, by_sample="invalid")

    # Test with empty VCF path
    with pytest.raises(ValueError, match="vcf_path.*cannot be empty"):
        missingness_plot("", output_path)

    # Test with None output path
    with pytest.raises(ValueError, match="output_path.*cannot be empty"):
        missingness_plot("dummy.vcf", None)
    assert result["status"] == "skipped"
    assert "message" in result


def test_transition_transversion_plot(tmp_path: Path) -> None:
    """Test transition/transversion ratio plot."""
    # Create sample GWAS results with REF/ALT alleles
    bases = ['A', 'C', 'G', 'T']
    results = []
    for _ in range(100):
        ref = np.random.choice(bases)
        alt = np.random.choice([b for b in bases if b != ref])
        results.append({"REF": ref, "ALT": alt})

    output_path = tmp_path / "tstv.png"
    result = transition_transversion_plot(results, output_path)

    assert result["status"] == "success"
    assert output_path.exists()
    assert "ts_tv_ratio" in result


def test_maf_distribution_empty(tmp_path: Path) -> None:
    """Test MAF distribution with empty data."""
    maf_data = np.array([])

    output_path = tmp_path / "maf_empty.png"
    result = maf_distribution(maf_data, output_path)

    # Should handle gracefully
    assert isinstance(result, dict)
    assert "status" in result
