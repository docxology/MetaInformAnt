"""Variant calling from BAM/CRAM files."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Any

from ..core.io import dump_json, ensure_directory

logger = logging.getLogger(__name__)


def check_bcftools_available() -> bool:
    """Check if bcftools is available on PATH."""
    try:
        result = subprocess.run(
            ["bcftools", "--version"],
            capture_output=True,
            check=True,
            timeout=10,
        )
        return result.returncode == 0
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return False


def check_gatk_available() -> bool:
    """Check if GATK is available on PATH or as gatk command."""
    try:
        result = subprocess.run(
            ["gatk", "--version"],
            capture_output=True,
            check=True,
            timeout=10,
        )
        return result.returncode == 0
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return False


def call_variants_bcftools(
    bam_files: list[str | Path],
    reference: str | Path,
    output_vcf: str | Path,
    *,
    threads: int = 1,
    region: str | None = None,
) -> dict[str, Any]:
    """Call variants using bcftools mpileup and call.

    Args:
        bam_files: List of BAM/CRAM file paths
        reference: Path to reference genome FASTA
        output_vcf: Path to output VCF file
        threads: Number of threads to use
        region: Optional genomic region (chr:start-end)

    Returns:
        Dictionary with calling metadata and result status
    """
    logger.info(f"call_variants_bcftools: Calling variants from {len(bam_files)} BAM file(s)")

    if not check_bcftools_available():
        return {
            "status": "failed",
            "error": "bcftools not found on PATH",
            "message": "bcftools is required for variant calling",
        }

    reference_obj = Path(reference)
    if not reference_obj.exists():
        return {
            "status": "failed",
            "error": f"Reference genome not found: {reference_obj}",
        }

    output_vcf_obj = Path(output_vcf)
    ensure_directory(output_vcf_obj.parent)

    # Step 1: mpileup to generate BCF
    bcf_intermediate = output_vcf_obj.with_suffix(".bcf")
    mpileup_cmd = [
        "bcftools",
        "mpileup",
        "-f",
        str(reference_obj),
        "-Ou",
    ]
    if region:
        mpileup_cmd.extend(["-r", region])
    mpileup_cmd.extend(["-o", str(bcf_intermediate)])
    mpileup_cmd.extend([str(Path(bam)) for bam in bam_files])

    # Step 2: call variants
    call_cmd = [
        "bcftools",
        "call",
        "-mv",
        "-Oz",
        "-o",
        str(output_vcf_obj),
        str(bcf_intermediate),
    ]

    try:
        # Run mpileup
        logger.info("Running bcftools mpileup...")
        mpileup_result = subprocess.run(
            mpileup_cmd,
            capture_output=True,
            text=True,
            timeout=7200,  # 2 hours
        )

        if mpileup_result.returncode != 0:
            return {
                "status": "failed",
                "error": mpileup_result.stderr,
                "return_code": mpileup_result.returncode,
            }

        # Run call
        logger.info("Running bcftools call...")
        call_result = subprocess.run(
            call_cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour
        )

        # Clean up intermediate BCF
        if bcf_intermediate.exists():
            bcf_intermediate.unlink()

        if call_result.returncode == 0:
            logger.info(f"call_variants_bcftools: Successfully called variants to {output_vcf_obj}")
            return {
                "status": "success",
                "method": "bcftools",
                "input_bams": [str(Path(bam)) for bam in bam_files],
                "reference": str(reference_obj),
                "output_vcf": str(output_vcf_obj),
                "return_code": 0,
            }
        else:
            return {
                "status": "failed",
                "error": call_result.stderr,
                "return_code": call_result.returncode,
            }

    except subprocess.TimeoutExpired:
        return {
            "status": "failed",
            "error": "Command timeout",
        }
    except Exception as exc:
        logger.error(f"call_variants_bcftools: Error: {exc}")
        return {
            "status": "failed",
            "error": str(exc),
        }


def call_variants_gatk(
    bam_files: list[str | Path],
    reference: str | Path,
    output_vcf: str | Path,
    *,
    intervals: str | None = None,
    threads: int = 1,
) -> dict[str, Any]:
    """Call variants using GATK HaplotypeCaller.

    Args:
        bam_files: List of BAM file paths
        reference: Path to reference genome FASTA
        output_vcf: Path to output VCF file
        intervals: Optional intervals file for targeted calling
        threads: Number of threads to use

    Returns:
        Dictionary with calling metadata and result status

    Note:
        GATK requires the reference to have a .dict and .fai index.
    """
    logger.info(f"call_variants_gatk: Calling variants from {len(bam_files)} BAM file(s)")

    if not check_gatk_available():
        return {
            "status": "failed",
            "error": "GATK not found on PATH",
            "message": "GATK (gatk command) is required for variant calling",
        }

    reference_obj = Path(reference)
    if not reference_obj.exists():
        return {
            "status": "failed",
            "error": f"Reference genome not found: {reference_obj}",
        }

    output_vcf_obj = Path(output_vcf)
    ensure_directory(output_vcf_obj.parent)

    # GATK HaplotypeCaller command
    gatk_cmd = [
        "gatk",
        "HaplotypeCaller",
        "-R",
        str(reference_obj),
        "-O",
        str(output_vcf_obj),
        "--native-pair-hmm-threads",
        str(threads),
    ]

    if intervals:
        gatk_cmd.extend(["-L", intervals])

    # Add input BAM files
    for bam in bam_files:
        bam_path = Path(bam)
        if not bam_path.exists():
            return {
                "status": "failed",
                "error": f"BAM file not found: {bam_path}",
            }
        gatk_cmd.extend(["-I", str(bam_path)])

    try:
        logger.info("Running GATK HaplotypeCaller...")
        result = subprocess.run(
            gatk_cmd,
            capture_output=True,
            text=True,
            timeout=14400,  # 4 hours
        )

        if result.returncode == 0:
            logger.info(f"call_variants_gatk: Successfully called variants to {output_vcf_obj}")
            return {
                "status": "success",
                "method": "gatk",
                "input_bams": [str(Path(bam)) for bam in bam_files],
                "reference": str(reference_obj),
                "output_vcf": str(output_vcf_obj),
                "return_code": 0,
            }
        else:
            return {
                "status": "failed",
                "error": result.stderr,
                "return_code": result.returncode,
            }

    except subprocess.TimeoutExpired:
        return {
            "status": "failed",
            "error": "Command timeout",
        }
    except Exception as exc:
        logger.error(f"call_variants_gatk: Error: {exc}")
        return {
            "status": "failed",
            "error": str(exc),
        }


def merge_vcf_files(
    vcf_list: list[str | Path],
    output_vcf: str | Path,
) -> dict[str, Any]:
    """Merge multiple VCF files into a single VCF.

    Args:
        vcf_list: List of VCF file paths to merge
        output_vcf: Path to output merged VCF file

    Returns:
        Dictionary with merge metadata and result status

    Note:
        Requires bcftools to be available. Files should have same sample IDs
        or be indexed appropriately.
    """
    logger.info(f"merge_vcf_files: Merging {len(vcf_list)} VCF file(s)")

    if not check_bcftools_available():
        return {
            "status": "failed",
            "error": "bcftools not found on PATH",
            "message": "bcftools is required for VCF merging",
        }

    if not vcf_list:
        return {
            "status": "failed",
            "error": "No VCF files provided",
        }

    output_vcf_obj = Path(output_vcf)
    ensure_directory(output_vcf_obj.parent)

    # Check all input files exist
    for vcf_path in vcf_list:
        if not Path(vcf_path).exists():
            return {
                "status": "failed",
                "error": f"VCF file not found: {vcf_path}",
            }

    # Use bcftools merge
    merge_cmd = [
        "bcftools",
        "merge",
        "-O",
        "z",
        "-o",
        str(output_vcf_obj),
    ]
    merge_cmd.extend([str(Path(vcf)) for vcf in vcf_list])

    try:
        logger.info("Running bcftools merge...")
        result = subprocess.run(
            merge_cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour
        )

        if result.returncode == 0:
            logger.info(f"merge_vcf_files: Successfully merged VCFs to {output_vcf_obj}")
            return {
                "status": "success",
                "input_vcfs": [str(Path(vcf)) for vcf in vcf_list],
                "output_vcf": str(output_vcf_obj),
                "return_code": 0,
            }
        else:
            return {
                "status": "failed",
                "error": result.stderr,
                "return_code": result.returncode,
            }

    except subprocess.TimeoutExpired:
        return {
            "status": "failed",
            "error": "Command timeout",
        }
    except Exception as exc:
        logger.error(f"merge_vcf_files: Error: {exc}")
        return {
            "status": "failed",
            "error": str(exc),
        }

