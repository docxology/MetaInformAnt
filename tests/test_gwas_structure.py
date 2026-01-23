"""Tests for GWAS population structure analysis."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas import compute_kinship_matrix, compute_pca, estimate_population_structure, parse_vcf_full


def test_compute_pca_basic() -> None:
    """Test PCA computation on genotype matrix."""
    # Create simple genotype matrix: 5 samples, 10 variants
    # Samples 0-2 are similar, samples 3-4 are similar
    genotypes = [
        [0, 0, 1, 1, 0, 0, 1, 1, 0, 0],  # Sample 0
        [0, 0, 1, 1, 0, 0, 1, 1, 0, 0],  # Sample 1 (similar to 0)
        [0, 0, 1, 1, 0, 0, 1, 1, 0, 0],  # Sample 2 (similar to 0)
        [2, 2, 1, 1, 2, 2, 1, 1, 2, 2],  # Sample 3 (different)
        [2, 2, 1, 1, 2, 2, 1, 1, 2, 2],  # Sample 4 (similar to 3)
    ]

    result = compute_pca(genotypes, n_components=3)

    assert result["status"] == "success"
    assert "pcs" in result
    assert len(result["pcs"]) == 5  # 5 samples
    assert len(result["pcs"][0]) == 3  # 3 components
    assert "explained_variance_ratio" in result
    assert len(result["explained_variance_ratio"]) == 3


def test_compute_pca_empty_matrix() -> None:
    """Test PCA with empty genotype matrix."""
    result = compute_pca([], n_components=3)
    assert result["status"] == "failed"
    assert "error" in result


def test_compute_pca_with_missing_data() -> None:
    """Test PCA with missing genotypes."""
    genotypes = [
        [0, 1, -1, 2, 0],  # Sample 0 with missing
        [1, 1, 1, 1, 1],
        [2, 0, 0, 2, 2],
    ]

    result = compute_pca(genotypes, n_components=2)
    assert result["status"] == "success"
    # Missing data should be imputed


def test_compute_kinship_matrix_vanraden() -> None:
    """Test kinship matrix computation using VanRaden method."""
    # Small genotype matrix
    genotypes = [
        [0, 1, 2, 0, 1],  # Sample 0
        [0, 1, 2, 0, 1],  # Sample 1 (identical to 0)
        [2, 1, 0, 2, 1],  # Sample 2 (different)
    ]

    result = compute_kinship_matrix(genotypes, method="vanraden")

    assert result["status"] == "success"
    assert "kinship_matrix" in result
    kinship = np.array(result["kinship_matrix"])
    assert kinship.shape == (3, 3)

    # Samples 0 and 1 should have high kinship (they're identical)
    assert kinship[0, 1] > 0.5  # High relatedness
    assert abs(kinship[0, 1] - kinship[1, 0]) < 1e-6  # Symmetric


def test_compute_kinship_matrix_astle() -> None:
    """Test kinship computation using Astle-Balding method."""
    genotypes = [
        [0, 1, 2],
        [0, 1, 2],
        [2, 1, 0],
    ]

    result = compute_kinship_matrix(genotypes, method="astle")
    assert result["status"] == "success"
    assert "kinship_matrix" in result


def test_compute_kinship_matrix_yang() -> None:
    """Test kinship computation using Yang method."""
    genotypes = [
        [0, 1, 2],
        [0, 1, 2],
        [2, 1, 0],
    ]

    result = compute_kinship_matrix(genotypes, method="yang")
    assert result["status"] == "success"
    assert "kinship_matrix" in result


def test_compute_kinship_matrix_invalid_method() -> None:
    """Test kinship computation with invalid method."""
    genotypes = [[0, 1, 2], [1, 2, 0]]

    result = compute_kinship_matrix(genotypes, method="invalid")
    assert result["status"] == "failed"
    assert "error" in result


def test_estimate_population_structure_integration(tmp_path: Path) -> None:
    """Test full population structure estimation from VCF."""
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
chr1	100	rs1	A	G	60	PASS	.	GT	0/1	1/1	0/0	0/1	0/0
chr1	200	rs2	T	C	80	PASS	.	GT	0/0	0/1	1/1	0/0	0/1
chr1	300	rs3	G	A	50	PASS	.	GT	1/1	0/0	0/1	1/1	0/0
chr1	400	rs4	C	T	70	PASS	.	GT	0/1	0/1	0/0	0/1	1/1
chr1	500	rs5	A	T	90	PASS	.	GT	0/0	1/1	0/1	0/0	1/1
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    structure_dir = tmp_path / "structure"
    structure_dir.mkdir()

    config = {
        "compute_pca": True,
        "n_components": 3,
        "compute_relatedness": True,
        "kinship_method": "vanraden",
    }

    result = estimate_population_structure(vcf_file, config, output_dir=structure_dir)

    assert result["status"] == "success"
    assert "pca" in result or "kinship" in result

    # Check output files if generated
    if structure_dir.exists():
        assert (structure_dir / "structure_summary.json").exists()
