"""Tests for DNA variant VCF parsing."""

from __future__ import annotations

from pathlib import Path

from metainformant.dna import variants


def test_parse_vcf_minimal(tmp_path: Path) -> None:
    """Test parsing a minimal VCF file."""
    vcf_text = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t1\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1
chr1\t2\trs2\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/1
"""
    p = tmp_path / "small.vcf"
    p.write_text(vcf_text)
    parsed = variants.parse_vcf(p)
    # Function returns 'total_variants' not 'num_variants'
    assert parsed["total_variants"] == 2
    assert parsed["samples"] == ["S1", "S2"]
    assert "variants" in parsed


def test_parse_vcf_returns_dict(tmp_path: Path) -> None:
    """Test that parse_vcf returns a dictionary."""
    vcf_text = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1\t.\tA\tG\t.\tPASS\t.
"""
    p = tmp_path / "minimal.vcf"
    p.write_text(vcf_text)
    parsed = variants.parse_vcf(p)
    assert isinstance(parsed, dict)
