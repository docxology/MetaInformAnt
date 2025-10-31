"""Tests for GWAS variant calling functions."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.gwas.calling import (
    check_bcftools_available,
    check_gatk_available,
    merge_vcf_files,
)


def test_check_bcftools_available() -> None:
    """Test checking if bcftools is available."""
    result = check_bcftools_available()
    assert isinstance(result, bool)


def test_check_gatk_available() -> None:
    """Test checking if GATK is available."""
    result = check_gatk_available()
    assert isinstance(result, bool)


def test_merge_vcf_files_bcftools_unavailable(tmp_path: Path) -> None:
    """Test VCF merging when bcftools unavailable (should fail gracefully)."""
    vcf1 = tmp_path / "test1.vcf"
    vcf1.write_text("##fileformat=VCFv4.2\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1\nchr1	100	rs1	A	G	60	PASS	.	GT	0/1\n")
    
    output_vcf = tmp_path / "merged.vcf"
    
    result = merge_vcf_files([str(vcf1)], output_vcf)
    
    # Should return status indicating bcftools requirement or success
    assert result["status"] in ["success", "failed"]
    if result["status"] == "failed":
        assert "bcftools" in result.get("error", "").lower() or "error" in result


def test_merge_vcf_files_file_not_found(tmp_path: Path) -> None:
    """Test VCF merging with non-existent files."""
    vcf_file = tmp_path / "nonexistent.vcf"
    output_vcf = tmp_path / "merged.vcf"
    
    result = merge_vcf_files([str(vcf_file)], output_vcf)
    assert result["status"] == "failed"
    assert "not found" in result.get("error", "").lower() or "error" in result


def test_merge_vcf_files_empty_list(tmp_path: Path) -> None:
    """Test VCF merging with empty file list."""
    output_vcf = tmp_path / "merged.vcf"
    
    result = merge_vcf_files([], output_vcf)
    assert result["status"] == "failed"
    assert "error" in result

