"""Variant and genome data download functionality."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..core.io import dump_json, ensure_directory
from ..core.logging import get_logger
from ..dna.ncbi import download_genome_package_best_effort

logger = get_logger(__name__)


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
    url: str | None = None,
) -> dict[str, Any]:
    """Download variant data from public databases.

    Supported sources:
    - dbSNP: NCBI dbSNP database (via FTP)
    - custom: Custom VCF URL
    - sra: Download SRA data and call variants

    Args:
        source: Data source (dbSNP, custom, sra)
        accession: Genome/assembly accession for dbSNP, or SRA run accession
        region: Optional genomic region filter (chr:start-end)
        dest_dir: Destination directory for variant files
        url: Direct URL for custom downloads

    Returns:
        Dictionary with download metadata and file paths
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
    }

    if source == "dbSNP":
        result.update(_download_from_dbsnp(accession, region, out_dir))
    elif source == "custom":
        if not url:
            raise ValueError("URL required for custom downloads")
        result.update(_download_from_url(url, out_dir))
    elif source == "sra":
        result["message"] = "SRA download requires separate workflow (download → align → call variants)"
        result["status"] = "pending"
        result["note"] = "Use SRA toolkit to download, then use variant calling workflow"
    else:
        raise ValueError(f"Unsupported variant source: {source}")

    dump_json(result, out_dir / "variant_download_record.json", indent=2)
    return result


def _download_from_dbsnp(
    accession: str | None,
    region: str | None,
    dest_dir: Path,
) -> dict[str, Any]:
    """Download variants from NCBI dbSNP via FTP.
    
    Args:
        accession: Assembly accession (e.g., GCF_003254395.2)
        region: Optional region filter
        dest_dir: Destination directory
    
    Returns:
        Status dictionary with file paths
    """
    import subprocess
    from urllib.parse import urlparse
    
    if not accession:
        return {
            "status": "failed",
            "error": "Assembly accession required for dbSNP downloads"
        }
    
    logger.info(f"Attempting dbSNP download for {accession}")
    
    # For Apis mellifera GCF_003254395.2, construct dbSNP FTP path
    # Format: ftp://ftp.ncbi.nih.gov/snp/organisms/<organism>/<assembly>/VCF/
    # However, not all assemblies have dbSNP data
    
    # Try to construct likely FTP paths
    assembly_parts = accession.split("_")
    if len(assembly_parts) >= 2:
        # Example: GCF_003254395.2 -> 003/254/395/GCF_003254395.2/
        gcf_num = assembly_parts[1].split(".")[0]
        if len(gcf_num) >= 9:
            path_parts = [gcf_num[i:i+3] for i in range(0, 9, 3)]
            ftp_path = f"ftp://ftp.ncbi.nih.gov/genomes/all/GCF/{'/'.join(path_parts)}/{accession}/"
            
            # Check for VCF files
            vcf_file = dest_dir / f"{accession}_dbsnp.vcf.gz"
            
            logger.info(f"Checking NCBI FTP for dbSNP VCF at {ftp_path}")
            
            # Try to download with wget or curl
            try:
                # Use wget to list directory
                list_cmd = ["wget", "-q", "-O", "-", "--spider", f"{ftp_path}"]
                proc = subprocess.run(list_cmd, capture_output=True, timeout=30, text=True)
                
                # Real-world note: Most assemblies don't have dbSNP VCF files directly
                # Users typically need to download from dbSNP FTP or use variant calling
                
                return {
                    "status": "failed",
                    "error": "dbSNP VCF files not directly available for this assembly",
                    "message": "Consider downloading SRA data and calling variants, or check dbSNP FTP manually",
                    "ftp_path_checked": ftp_path,
                    "note": "dbSNP data for non-model organisms often requires variant calling from raw data"
                }
                
            except (subprocess.TimeoutExpired, subprocess.CalledProcessError, FileNotFoundError) as e:
                logger.warning(f"Failed to check NCBI FTP: {e}")
                return {
                    "status": "failed",
                    "error": str(e),
                    "message": "Could not access NCBI FTP"
                }
    
    return {
        "status": "failed",
        "error": "Invalid accession format",
        "message": f"Could not parse accession: {accession}"
    }


def _download_from_url(url: str, dest_dir: Path) -> dict[str, Any]:
    """Download VCF file from a direct URL.
    
    Args:
        url: Direct URL to VCF file
        dest_dir: Destination directory
    
    Returns:
        Status dictionary with file paths
    """
    import subprocess
    from urllib.parse import urlparse
    
    logger.info(f"Downloading VCF from {url}")
    
    # Extract filename from URL
    parsed = urlparse(url)
    filename = Path(parsed.path).name
    if not filename:
        filename = "downloaded.vcf.gz"
    
    output_file = dest_dir / filename
    
    # Try wget first, then curl
    for cmd_template in [
        ["wget", "-O", str(output_file), url],
        ["curl", "-L", "-o", str(output_file), url],
    ]:
        try:
            logger.info(f"Attempting download with {cmd_template[0]}")
            proc = subprocess.run(
                cmd_template,
                capture_output=True,
                text=True,
                timeout=300,  # 5 minutes
                check=True,
            )
            
            if output_file.exists() and output_file.stat().st_size > 0:
                logger.info(f"Successfully downloaded {output_file.stat().st_size} bytes")
                return {
                    "status": "success",
                    "vcf_file": str(output_file),
                    "size_bytes": output_file.stat().st_size,
                    "url": url,
                }
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.warning(f"{cmd_template[0]} failed: {e}")
            continue
    
    return {
        "status": "failed",
        "error": "Could not download file with wget or curl",
        "message": "Ensure wget or curl is installed",
        "url": url,
    }


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

