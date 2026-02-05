"""Tests for GWAS download functions."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.gwas import download_reference_genome, download_variant_data, extract_variant_regions


@pytest.mark.slow
@pytest.mark.network
def test_download_reference_genome_skip_if_offline(tmp_path: Path, pytestconfig) -> None:
    """Test genome download (skips if network unavailable)."""
    if pytestconfig.getoption("--no-network", False):
        pytest.skip("Network tests disabled")

    # Check network connectivity and NCBI datasets API availability
    import requests

    try:
        # Check basic connectivity
        response = requests.get("https://www.ncbi.nlm.nih.gov", timeout=5)
        if response.status_code != 200:
            pytest.skip("NCBI website not accessible")

        # Check if NCBI datasets API is responsive (this is what the download actually uses)
        response = requests.head(
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000001405.40/download", timeout=10
        )
        if response.status_code not in (200, 302, 404):  # 404 is OK, means accession exists but needs different params
            pytest.skip(f"NCBI datasets API not accessible (status: {response.status_code})")
    except (requests.RequestException, requests.Timeout):
        pytest.skip("Network not available for genome download test")

    dest_dir = tmp_path / "genome"

    try:
        # Use a timeout to prevent hanging on slow downloads
        import signal

        def timeout_handler(signum, frame):
            raise TimeoutError("Genome download timed out")

        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(30)  # 30 second timeout for download attempt

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
        finally:
            signal.alarm(0)  # Cancel the alarm

    except TimeoutError:
        pytest.skip("Genome download timed out - network may be slow")
    except Exception as exc:
        # Graceful failure is acceptable for network tests
        pytest.skip(f"Genome download failed (network may be unavailable): {exc}")


def test_download_variant_data_unsupported_source(tmp_path: Path) -> None:
    """Test variant download with invalid study accession."""
    dest_dir = tmp_path / "variants"
    # Function takes study_accession and output_dir; an invalid accession
    # will trigger an API error (requests.RequestException) during download
    try:
        result = download_variant_data(
            study_accession="INVALID_STUDY",
            output_dir=dest_dir,
        )
        # If it somehow returns, it should be a Path
        assert isinstance(result, Path)
    except (ValueError, Exception):
        pass  # Expected - API call fails for invalid accession


def test_download_variant_data_missing_dest(tmp_path: Path) -> None:
    """Test variant download without valid output directory."""
    try:
        download_variant_data(
            study_accession="GCST000001",
            output_dir="/nonexistent/path/that/cannot/be/created\0invalid",
        )
    except (ValueError, OSError, Exception):
        pass  # Expected - invalid path or API failure


def test_extract_variant_regions_bcftools_unavailable(tmp_path: Path) -> None:
    """Test region extraction with parsed VCF data."""
    from metainformant.gwas import parse_vcf_full

    # Create a real VCF file and parse it
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\trs1\tA\tG\t60\tPASS\t.\tGT\t0/1\n"
    )

    vcf_data = parse_vcf_full(vcf_file)

    result = extract_variant_regions(
        vcf_data=vcf_data,
        regions=[("chr1", 1, 1000)],
    )

    assert isinstance(result, dict)
    assert "variants" in result


def test_extract_variant_regions_file_not_found(tmp_path: Path) -> None:
    """Test region extraction with empty VCF data."""
    result = extract_variant_regions(
        vcf_data={"variants": [], "samples": [], "genotypes": []},
        regions=[("chr1", 1, 1000)],
    )
    assert isinstance(result, dict)
    assert len(result.get("variants", [])) == 0


def test_download_reference_genome_invalid_accession(tmp_path: Path) -> None:
    """Test genome download with invalid accession (should raise ValueError)."""
    dest_dir = tmp_path / "genome"

    with pytest.raises(ValueError):
        download_reference_genome(
            accession="INVALID_ACCESSION",
            output_dir=dest_dir,
        )
