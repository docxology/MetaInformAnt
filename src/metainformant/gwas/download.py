"""Variant and genome data download functionality."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from ..core.io import dump_json, ensure_directory
from ..dna.ncbi import download_genome_package_best_effort

logger = logging.getLogger(__name__)


def download_reference_genome(
    accession: str,
    dest_dir: str | Path,
    *,
    include: list[str] | None = None,
    ftp_url: str | None = None,
) -> dict[str, Any]:
    """Download reference genome using NCBI infrastructure.

    Args:
        accession: NCBI assembly accession (e.g., GCF_000001405.40)
        dest_dir: Destination directory for genome files
        include: List of file types to download (genome, gff3, etc.)
        ftp_url: Optional FTP URL for direct download

    Returns:
        Dictionary with download metadata and paths
    """
    logger.info(f"download_reference_genome: Downloading genome {accession} to {dest_dir}")
    out_dir = ensure_directory(dest_dir)

    include_vals = include or ["genome", "gff3"]
    record = download_genome_package_best_effort(
        accession=accession,
        dest_dir=out_dir,
        include=include_vals,
        ftp_url=ftp_url,
    )

    result = {
        "accession": accession,
        "dest_dir": str(out_dir),
        "method": record.get("method", "unknown"),
        "return_code": record.get("return_code", 1),
        "extracted_dir": record.get("extracted_dir", ""),
    }

    if record.get("return_code") == 0:
        logger.info(f"download_reference_genome: Successfully downloaded genome {accession}")
        result["status"] = "success"
    else:
        logger.warning(f"download_reference_genome: Download failed for {accession}")
        result["status"] = "failed"
        result["error"] = record.get("error", "Unknown error")

    # Write download record
    dump_json(result, out_dir / "genome_download_record.json", indent=2)
    return result


def download_variant_data(
    source: str,
    accession: str | None = None,
    region: str | None = None,
    dest_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Download variant data from public databases.

    Supported sources:
    - dbSNP: NCBI dbSNP database
    - 1000genomes: 1000 Genomes Project data
    - custom: Custom VCF URL

    Args:
        source: Data source (dbSNP, 1000genomes, custom)
        accession: Genome/assembly accession for dbSNP queries
        region: Optional genomic region filter (chr:start-end)
        dest_dir: Destination directory for variant files

    Returns:
        Dictionary with download metadata and file paths

    Note:
        This is a placeholder for future implementation. Currently
        returns a status indicating the feature is not yet implemented.
    """
    logger.info(f"download_variant_data: Downloading variants from {source}")

    if dest_dir:
        out_dir = ensure_directory(dest_dir)
    else:
        raise ValueError("dest_dir is required for variant downloads")

    result = {
        "source": source,
        "accession": accession,
        "region": region,
        "dest_dir": str(out_dir),
        "status": "pending",
        "message": "Variant download from public databases not yet fully implemented",
    }

    if source == "dbSNP":
        # TODO: Implement dbSNP download via NCBI API
        logger.warning("dbSNP download not yet implemented")
        result["message"] = "dbSNP download requires NCBI API integration"
    elif source == "1000genomes":
        # TODO: Implement 1000 Genomes download
        logger.warning("1000 Genomes download not yet implemented")
        result["message"] = "1000 Genomes download requires IGSR API integration"
    elif source == "custom":
        # TODO: Support custom VCF URL downloads
        logger.warning("Custom VCF download not yet implemented")
        result["message"] = "Custom VCF download not yet implemented"
    else:
        raise ValueError(f"Unsupported variant source: {source}")

    dump_json(result, out_dir / "variant_download_record.json", indent=2)
    return result


def extract_variant_regions(
    vcf_path: str | Path,
    regions: list[str],
    output_vcf: str | Path,
) -> dict[str, Any]:
    """Extract specific genomic regions from VCF file.

    Args:
        vcf_path: Path to input VCF file
        regions: List of genomic regions (format: "chr:start-end")
        output_vcf: Path to output VCF file

    Returns:
        Dictionary with extraction metadata

    Note:
        This requires bcftools to be available on PATH.
        Returns status indicating if extraction succeeded.
    """
    import subprocess

    logger.info(f"extract_variant_regions: Extracting {len(regions)} regions from {vcf_path}")

    vcf_path_obj = Path(vcf_path)
    output_vcf_obj = Path(output_vcf)

    if not vcf_path_obj.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path_obj}")

    # Check if bcftools is available
    try:
        subprocess.run(["bcftools", "--version"], capture_output=True, check=True, timeout=10)
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        logger.warning("bcftools not available, cannot extract regions")
        return {
            "status": "failed",
            "error": "bcftools not found on PATH",
            "message": "bcftools is required for region extraction",
        }

    # Build bcftools view command with regions
    cmd = ["bcftools", "view", "-O", "v", str(vcf_path_obj)]
    for region in regions:
        cmd.extend(["-r", region])

    try:
        ensure_directory(output_vcf_obj.parent)
        with open(output_vcf_obj, "w") as out_fh:
            result = subprocess.run(
                cmd,
                stdout=out_fh,
                stderr=subprocess.PIPE,
                text=True,
                timeout=3600,
            )

        if result.returncode == 0:
            logger.info(f"extract_variant_regions: Successfully extracted regions to {output_vcf_obj}")
            return {
                "status": "success",
                "input_vcf": str(vcf_path_obj),
                "output_vcf": str(output_vcf_obj),
                "regions": regions,
                "return_code": 0,
            }
        else:
            logger.error(f"extract_variant_regions: bcftools failed: {result.stderr}")
            return {
                "status": "failed",
                "error": result.stderr,
                "return_code": result.returncode,
            }

    except subprocess.TimeoutExpired:
        logger.error("extract_variant_regions: bcftools command timed out")
        return {
            "status": "failed",
            "error": "Command timeout",
        }
    except Exception as exc:
        logger.error(f"extract_variant_regions: Error running bcftools: {exc}")
        return {
            "status": "failed",
            "error": str(exc),
        }

