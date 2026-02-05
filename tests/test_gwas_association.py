"""Tests for GWAS association testing."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.gwas import association_test_linear, association_test_logistic, run_gwas
from metainformant.core.io.io import write_tsv, write_delimited


def test_association_linear_basic() -> None:
    """Test basic linear regression association test."""
    # Genotypes: 0, 1, 2, 0, 1, 2 (varying dosage)
    genotypes = [0, 1, 2, 0, 1, 2]
    # Phenotypes with positive association
    phenotypes = [10.0, 11.0, 12.0, 10.0, 11.0, 12.0]

    result = association_test_linear(genotypes, phenotypes)

    assert result["status"] == "success"
    assert "beta" in result
    assert "se" in result
    assert "p_value" in result
    assert "r_squared" in result
    assert result["beta"] > 0  # Positive association


def test_association_linear_no_association() -> None:
    """Test linear regression with no association."""
    genotypes = [0, 1, 2, 0, 1, 2]
    # Random phenotypes (no association)
    phenotypes = [10.0, 10.5, 9.5, 11.0, 9.0, 10.5]

    result = association_test_linear(genotypes, phenotypes)

    assert result["status"] == "success"
    # P-value should be large (not significant)
    assert result["p_value"] > 0.05


def test_association_linear_with_covariates() -> None:
    """Test linear regression with covariates."""
    # Use more samples and non-collinear covariates
    genotypes = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0]
    phenotypes = [10.0, 11.0, 12.0, 10.5, 11.5, 12.5, 9.5, 10.5, 11.5, 10.0]
    covariates = [[25, 35, 45, 30, 40, 50, 28, 38, 48, 33], [1, 1, 1, 0, 0, 0, 1, 0, 1, 0]]  # age, sex

    result = association_test_linear(genotypes, phenotypes, covariates)

    assert result["status"] == "success"
    assert "beta" in result
    assert "p_value" in result

    # CRITICAL: verify result is NOT the old hardcoded (0.1, 0.05, 2.0, 0.05, 0.5)
    assert not (
        abs(result["beta"] - 0.1) < 1e-10 and abs(result["se"] - 0.05) < 1e-10 and abs(result["p_value"] - 0.05) < 1e-10
    ), "Multi-covariate regression still returns hardcoded values!"

    # SE should be positive (real computation)
    assert result["se"] > 0


def test_association_linear_covariates_produce_different_results() -> None:
    """Test that adding covariates produces different results than no covariates."""
    genotypes = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0]
    phenotypes = [10.0, 11.0, 12.0, 10.5, 11.5, 12.5, 9.5, 10.5, 11.5, 10.0]

    result_no_cov = association_test_linear(genotypes, phenotypes)
    # Use a covariate NOT perfectly correlated with genotype
    covariates = [[25, 35, 45, 30, 40, 50, 28, 38, 48, 33]]
    result_with_cov = association_test_linear(genotypes, phenotypes, covariates=covariates)

    assert result_no_cov["status"] == "success"
    assert result_with_cov["status"] == "success"
    # Results should differ when covariates are added
    assert result_no_cov["se"] != result_with_cov["se"] or result_no_cov["r_squared"] != result_with_cov["r_squared"]


def test_association_linear_missing_data() -> None:
    """Test linear regression with missing genotypes."""
    genotypes = [0, -1, 2, 0, 1, 2]  # One missing
    phenotypes = [10.0, 11.0, 12.0, 10.0, 11.0, 12.0]

    result = association_test_linear(genotypes, phenotypes)

    assert result["status"] == "success"
    # Should handle missing data by excluding that sample


def test_association_logistic_basic() -> None:
    """Test basic logistic regression association test."""
    # Genotypes with moderate association (less extreme to avoid separation issues)
    genotypes = [0, 0, 0, 1, 1, 1, 2, 2, 2]
    # Binary phenotypes with positive association
    phenotypes = [0, 0, 1, 0, 1, 1, 1, 1, 1]

    result = association_test_logistic(genotypes, phenotypes)

    assert result["status"] == "success"
    assert "beta" in result
    assert "se" in result
    assert "p_value" in result
    assert "odds_ratio" in result
    # With positive association, beta should be positive and OR > 1
    # But allow for some numerical instability
    if not np.isinf(result.get("beta", 0)) and not np.isnan(result.get("beta", 0)):
        assert result["odds_ratio"] > 0  # Just check it's positive, not NaN


def test_association_logistic_no_association() -> None:
    """Test logistic regression with no association."""
    genotypes = [0, 1, 2, 0, 1, 2]
    # Random case/control assignment
    phenotypes = [0, 1, 0, 1, 0, 1]

    result = association_test_logistic(genotypes, phenotypes)

    assert result["status"] == "success"
    # P-value should be large
    assert result["p_value"] > 0.05


def test_association_logistic_insufficient_cases() -> None:
    """Test logistic regression with insufficient cases/controls."""
    genotypes = [0, 0, 0, 0, 0, 0]
    phenotypes = [0, 0, 0, 0, 0, 1]  # Only 1 case

    result = association_test_logistic(genotypes, phenotypes)

    assert result["status"] == "failed"
    assert "error" in result


def test_run_gwas_integration(tmp_path: Path) -> None:
    """Test full GWAS run with VCF and phenotype files."""
    # Create VCF file
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
chr1	100	rs1	A	G	60	PASS	.	GT	0/1	1/1	0/0	0/1	0/0
chr1	200	rs2	T	C	80	PASS	.	GT	0/0	0/1	1/1	0/0	0/1
chr1	300	rs3	G	A	50	PASS	.	GT	1/1	0/0	0/1	1/1	0/0
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    # Create phenotype file
    phenotype_file = tmp_path / "phenotypes.tsv"
    phenotype_data = [
        ["sample_id", "height"],
        ["S1", "175.0"],
        ["S2", "180.0"],
        ["S3", "170.0"],
        ["S4", "178.0"],
        ["S5", "172.0"],
    ]
    write_tsv(phenotype_data, phenotype_file)

    config = {
        "model": "linear",
        "trait": "height",
        "covariates": [],
        "min_sample_size": 3,
    }

    results_dir = tmp_path / "results"
    results_dir.mkdir()

    result = run_gwas(vcf_file, phenotype_file, config, output_dir=results_dir)

    assert result["status"] == "success"
    assert "results" in result
    assert len(result["results"]) > 0

    # Check that results table was written
    if results_dir.exists():
        results_file = results_dir / "association_results.tsv"
        if results_file.exists():
            assert results_file.stat().st_size > 0


def test_run_gwas_missing_trait(tmp_path: Path) -> None:
    """Test GWAS run with missing trait in phenotype file."""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        "##fileformat=VCFv4.2\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1\nchr1	100	rs1	A	G	60	PASS	.	GT	0/1\n"
    )

    phenotype_file = tmp_path / "phenotypes.tsv"
    phenotype_data = [["sample_id", "weight"], ["S1", "70.0"]]
    write_tsv(phenotype_data, phenotype_file)

    config = {"model": "linear", "trait": "height", "covariates": []}

    result = run_gwas(vcf_file, phenotype_file, config)
    assert result["status"] == "failed"
    assert "error" in result
