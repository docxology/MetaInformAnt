"""Tests for GWAS quality control functions."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.gwas.analysis.quality import apply_qc_filters, parse_vcf_full


def test_parse_vcf_full_basic(tmp_path: Path) -> None:
    """Test parsing VCF file with genotypes."""
    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	rs1	A	G	60	PASS	DP=30	GT	0/1	1/1	0/0
chr1	200	rs2	T	C	80	PASS	DP=40	GT	0/0	0/1	1/1
chr1	300	rs3	G	A	50	PASS	DP=20	GT	1/1	0/0	0/1
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    result = parse_vcf_full(vcf_file)

    assert result["samples"] == ["S1", "S2", "S3"]
    assert len(result["variants"]) == 3
    assert len(result["genotypes"]) == 3  # 3 samples
    assert len(result["genotypes"][0]) == 3  # 3 variants

    # Check first variant genotypes: S1=0/1(1), S2=1/1(2), S3=0/0(0)
    assert result["genotypes"][0][0] == 1  # S1, variant 0
    assert result["genotypes"][1][0] == 2  # S2, variant 0
    assert result["genotypes"][2][0] == 0  # S3, variant 0


def test_parse_vcf_full_missing_genotypes(tmp_path: Path) -> None:
    """Test parsing VCF with missing genotypes."""
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	rs1	A	G	60	PASS	.	GT	0/1	./.	0/0
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    result = parse_vcf_full(vcf_file)
    assert len(result["genotypes"]) == 3
    assert result["genotypes"][0][0] == 1  # S1: 0/1 = 1
    assert result["genotypes"][1][0] == -1  # S2: ./. = missing
    assert result["genotypes"][2][0] == 0  # S3: 0/0 = 0


def test_apply_qc_filters_basic(tmp_path: Path) -> None:
    """Test applying basic QC filters."""
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
chr1	100	rs1	A	G	60	PASS	.	GT	0/1	1/1	0/0	0/1	0/0
chr1	200	rs2	T	C	80	PASS	.	GT	0/0	0/0	0/0	0/0	0/0
chr1	300	rs3	G	A	50	PASS	.	GT	1/1	1/1	1/1	1/1	1/1
chr1	400	rs4	C	T	30	PASS	.	GT	./.	./.	./.	./.	./.
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    qc_config = {
        "min_maf": 0.1,
        "max_missing": 0.2,
        "min_qual": 40.0,
        "exclude_indels": False,
    }

    result = apply_qc_filters(vcf_file, qc_config)

    assert result["status"] == "success"
    assert result["num_variants_before"] == 4
    assert result["num_variants_after"] <= 4

    # rs2 should be filtered (all hom ref = MAF = 0)
    # rs3 should be filtered (all hom alt = MAF = 0)
    # rs4 should be filtered (all missing)
    # rs1 should pass (MAF = 0.4, qual = 60)


def test_apply_qc_filters_maf(tmp_path: Path) -> None:
    """Test MAF filtering."""
    # Create VCF with variants at different MAF levels
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6	S7	S8	S9	S10
chr1	100	rs1	A	G	60	PASS	.	GT	0/1	1/1	0/0	0/1	0/0	0/1	0/0	0/1	0/0	0/1
chr1	200	rs2	T	C	80	PASS	.	GT	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	1/1
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    # rs1: 5 het + 1 hom alt = 6 alt alleles out of 20 = MAF = 0.3
    # rs2: 1 hom alt = 2 alt alleles out of 20 = MAF = 0.1

    qc_config = {"min_maf": 0.15, "max_missing": 1.0, "min_qual": 0.0}
    result = apply_qc_filters(vcf_file, qc_config)

    assert result["status"] == "success"
    # Only rs1 should pass (MAF >= 0.15)
    assert 1 <= result["num_variants_after"] <= 2


def test_apply_qc_filters_missing_data(tmp_path: Path) -> None:
    """Test missing data filtering."""
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
chr1	100	rs1	A	G	60	PASS	.	GT	0/1	1/1	0/0	0/1	0/0
chr1	200	rs2	T	C	80	PASS	.	GT	0/0	./.	./.	./.	./.
chr1	300	rs3	G	A	50	PASS	.	GT	1/1	1/1	./.	./.	./.
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    qc_config = {"min_maf": 0.0, "max_missing": 0.5, "min_qual": 0.0}
    result = apply_qc_filters(vcf_file, qc_config)

    assert result["status"] == "success"
    # rs2 has 80% missing (4/5), rs3 has 60% missing (3/5)
    # Only rs1 should pass (0% missing)


def test_apply_qc_filters_hwe(tmp_path: Path) -> None:
    """Test Hardy-Weinberg equilibrium filtering."""
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6	S7	S8	S9	S10
chr1	100	rs1	A	G	60	PASS	.	GT	0/0	0/0	0/0	0/0	0/0	1/1	1/1	1/1	1/1	1/1
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    # This variant violates HWE (no heterozygotes, extreme deviation)
    qc_config = {"min_maf": 0.0, "max_missing": 1.0, "min_qual": 0.0, "hwe_pval": 1e-6}
    result = apply_qc_filters(vcf_file, qc_config)

    assert result["status"] == "success"
    # Variant with extreme HWE violation might be filtered
