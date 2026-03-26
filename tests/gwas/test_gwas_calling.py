"""Tests for GWAS variant calling functions."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from metainformant.gwas.analysis.calling import (
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


def test_merge_vcf_files_single_file_error(tmp_path: Path) -> None:
    """Test VCF merging with a single file — verifies real error handling.

    When bcftools IS available: merge of 1 file raises RuntimeError
    (bcftools requires >=2 files without --force-single).
    When bcftools is NOT available: raises FileNotFoundError.
    Both are correct error handling behaviors.
    """
    vcf1 = tmp_path / "test1.vcf.gz"
    # Write a minimal valid .vcf.gz (plain text works for the existence check)
    vcf1.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\trs1\tA\tG\t60\tPASS\t.\tGT\t0/1\n"
    )

    output_vcf = tmp_path / "merged.vcf"

    # Should raise an error regardless of bcftools availability:
    # - FileNotFoundError: bcftools not installed
    # - RuntimeError: bcftools merge fails
    # - CalledProcessError: tabix indexing fails on non-BGZF file
    with pytest.raises((FileNotFoundError, RuntimeError, subprocess.CalledProcessError)):
        merge_vcf_files([str(vcf1)], output_vcf)


def test_merge_vcf_files_file_not_found(tmp_path: Path) -> None:
    """Test VCF merging with non-existent files."""
    vcf_file = tmp_path / "nonexistent.vcf"
    output_vcf = tmp_path / "merged.vcf"

    with pytest.raises(FileNotFoundError, match="VCF file not found"):
        merge_vcf_files([str(vcf_file)], output_vcf)


def test_merge_vcf_files_empty_list(tmp_path: Path) -> None:
    """Test VCF merging with empty file list."""
    output_vcf = tmp_path / "merged.vcf"

    with pytest.raises(ValueError, match="Must provide at least one VCF file"):
        merge_vcf_files([], output_vcf)
