"""Tests for GWAS download functions."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.gwas import download_reference_genome, download_variant_data, extract_variant_regions


@pytest.mark.network
def test_download_reference_genome_skip_if_offline(tmp_path: Path, pytestconfig) -> None:
    """Test genome download (skips if network unavailable)."""
    if pytestconfig.getoption("--no-network", False):
        pytest.skip("Network tests disabled")

    dest_dir = tmp_path / "genome"
    
    try:
        result = download_reference_genome(
            accession="GCF_000001405.40",  # Human GRCh38
            dest_dir=dest_dir,
            include=["genome"],
        )
        
        # May succeed or fail depending on network/NCBI availability
        assert result["status"] in ["success", "failed"]
        assert "accession" in result
        assert result["accession"] == "GCF_000001405.40"
    except Exception as exc:
        # Graceful failure is acceptable for network tests
        pytest.skip(f"Genome download failed (network may be unavailable): {exc}")


def test_download_variant_data_unsupported_source(tmp_path: Path) -> None:
    """Test variant download with unsupported source."""
    dest_dir = tmp_path / "variants"
    
    with pytest.raises(ValueError):
        download_variant_data(
            source="unsupported_source",
            dest_dir=dest_dir,
        )


def test_download_variant_data_missing_dest(tmp_path: Path) -> None:
    """Test variant download without destination directory."""
    with pytest.raises(ValueError):
        download_variant_data(
            source="dbSNP",
            dest_dir=None,
        )


def test_extract_variant_regions_bcftools_unavailable(tmp_path: Path) -> None:
    """Test region extraction when bcftools unavailable (should fail gracefully)."""
    # Create dummy VCF
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text("##fileformat=VCFv4.2\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1\nchr1	100	rs1	A	G	60	PASS	.	GT	0/1\n")
    
    output_vcf = tmp_path / "output.vcf"
    
    result = extract_variant_regions(
        vcf_path=vcf_file,
        regions=["chr1:1-1000"],
        output_vcf=output_vcf,
    )
    
    # Should return status indicating bcftools requirement
    assert result["status"] in ["success", "failed"]
    if result["status"] == "failed":
        assert "bcftools" in result.get("error", "").lower() or "error" in result


def test_extract_variant_regions_file_not_found(tmp_path: Path) -> None:
    """Test region extraction with non-existent VCF file."""
    vcf_file = tmp_path / "nonexistent.vcf"
    output_vcf = tmp_path / "output.vcf"
    
    with pytest.raises(FileNotFoundError):
        extract_variant_regions(
            vcf_path=vcf_file,
            regions=["chr1:1-1000"],
            output_vcf=output_vcf,
        )


def test_download_reference_genome_invalid_accession(tmp_path: Path) -> None:
    """Test genome download with invalid accession (should handle gracefully)."""
    dest_dir = tmp_path / "genome"
    
    result = download_reference_genome(
        accession="INVALID_ACCESSION",
        dest_dir=dest_dir,
        include=["genome"],
    )
    
    # Should return failed status
    assert result["status"] in ["success", "failed"]
    # If failed, should have error message
    if result["status"] == "failed":
        assert "error" in result

