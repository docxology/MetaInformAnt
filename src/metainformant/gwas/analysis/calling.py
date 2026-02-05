"""Variant calling integration for GWAS analysis.

This module provides tools for calling genetic variants from sequencing data
using bcftools and GATK, with support for various calling strategies and
quality control.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)


def call_variants_bcftools(
    bam_files: List[str | Path], reference_fasta: str | Path, output_vcf: str | Path, threads: int = 1
) -> subprocess.CompletedProcess:
    """Call variants using bcftools mpileup and call.

    Args:
        bam_files: List of BAM files to call variants from
        reference_fasta: Reference genome FASTA file
        output_vcf: Output VCF file path
        threads: Number of threads to use

    Returns:
        CompletedProcess with command execution results

    Raises:
        FileNotFoundError: If bcftools is not available
        subprocess.CalledProcessError: If variant calling fails

    Example:
        >>> # Assuming BAM files and reference exist
        >>> result = call_variants_bcftools(["sample.bam"], "ref.fa", "output.vcf")
        >>> result.returncode == 0
        True
    """
    # Check if bcftools is available
    try:
        subprocess.run(["bcftools", "--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise FileNotFoundError("bcftools not found on PATH")

    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    bam_files_str = [str(f) for f in bam_files]

    # Build bcftools mpileup command
    mpileup_cmd = [
        "bcftools",
        "mpileup",
        "-f",
        str(reference_fasta),
        "-Ou",  # output uncompressed BCF
        "--threads",
        str(threads),
        "--max-depth",
        "500",  # prevent excessive memory usage
        "--min-BQ",
        "20",  # minimum base quality
    ]

    # Add BAM files
    mpileup_cmd.extend(bam_files_str)

    # Build bcftools call command
    call_cmd = [
        "bcftools",
        "call",
        "-mv",  # multiallelic caller, output variants only
        "-Ou",  # output uncompressed BCF
        "--threads",
        str(threads),
        "-o",
        str(output_vcf.with_suffix(".bcf")),
    ]

    # Build bcftools view command to convert BCF to VCF
    view_cmd = [
        "bcftools",
        "view",
        "-o",
        str(output_vcf),
        "-O",
        "v",  # output uncompressed VCF
        str(output_vcf.with_suffix(".bcf")),
    ]

    logger.info(f"Running bcftools variant calling: {' '.join(mpileup_cmd[:5])}...")

    try:
        # Pipe mpileup output into call (no shell=True needed)
        mpileup_proc = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        result = subprocess.run(call_cmd, stdin=mpileup_proc.stdout, capture_output=True, text=True)
        mpileup_proc.stdout.close()
        mpileup_proc.wait()
        result.check_returncode()

        # Convert BCF to VCF
        view_result = subprocess.run(view_cmd, capture_output=True, text=True)
        view_result.check_returncode()

        logger.info(f"bcftools variant calling completed: {output_vcf}")

        # Clean up intermediate BCF file
        bcf_file = output_vcf.with_suffix(".bcf")
        if bcf_file.exists():
            bcf_file.unlink()

        return view_result

    except subprocess.CalledProcessError as e:
        logger.error(f"bcftools variant calling failed: {e}")
        logger.error(f"stderr: {e.stderr}")
        raise


def call_variants_gatk(
    bam_files: List[str | Path], reference_fasta: str | Path, output_vcf: str | Path, gatk_path: Optional[str] = None
) -> subprocess.CompletedProcess:
    """Call variants using GATK HaplotypeCaller.

    Args:
        bam_files: List of BAM files to call variants from
        reference_fasta: Reference genome FASTA file
        output_vcf: Output VCF file path
        gatk_path: Path to GATK jar file (if not on PATH)

    Returns:
        CompletedProcess with command execution results

    Raises:
        FileNotFoundError: If GATK is not available
        subprocess.CalledProcessError: If variant calling fails

    Example:
        >>> # Assuming BAM files and reference exist
        >>> result = call_variants_gatk(["sample.bam"], "ref.fa", "output.vcf")
        >>> result.returncode == 0
        True
    """
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    # Check if GATK is available
    gatk_cmd = [gatk_path] if gatk_path else ["gatk"]
    gatk_cmd.extend(["--version"])

    try:
        subprocess.run(gatk_cmd, capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise FileNotFoundError("GATK not found on PATH")

    # Build GATK command
    cmd = gatk_cmd[:-1] + [  # Remove --version
        "HaplotypeCaller",
        "-R",
        str(reference_fasta),
        "-O",
        str(output_vcf),
        "--native-pair-hmm-threads",
        "4",  # Use 4 threads for pair HMM
        "--standard-min-confidence-threshold-for-calling",
        "20.0",
        "--emit-ref-confidence",
        "GVCF",  # Generate GVCF
    ]

    # Add BAM files
    for bam_file in bam_files:
        cmd.extend(["-I", str(bam_file)])

    logger.info(f"Running GATK HaplotypeCaller: {' '.join(cmd[:5])}...")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        result.check_returncode()

        logger.info(f"GATK variant calling completed: {output_vcf}")
        return result

    except subprocess.CalledProcessError as e:
        logger.error(f"GATK variant calling failed: {e}")
        logger.error(f"stderr: {e.stderr}")
        raise


def filter_variants_bcftools(
    vcf_path: str | Path, output_vcf: str | Path, filters: Dict[str, Any]
) -> subprocess.CompletedProcess:
    """Filter variants using bcftools.

    Args:
        vcf_path: Input VCF file
        output_vcf: Output filtered VCF file
        filters: Dictionary of filter parameters

    Returns:
        CompletedProcess with filtering results

    Example:
        >>> filters = {'min_qual': 20, 'max_missing': 0.1}
        >>> result = filter_variants_bcftools("input.vcf", "filtered.vcf", filters)
    """
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    cmd = ["bcftools", "filter", "-o", str(output_vcf), "-O", "v"]  # uncompressed VCF

    # Build filter expressions
    filter_expr = []

    if "min_qual" in filters:
        filter_expr.append(f'QUAL>={filters["min_qual"]}')

    if "max_missing" in filters:
        max_missing = filters["max_missing"]
        filter_expr.append(f"F_MISSING<={max_missing}")

    if "min_depth" in filters:
        min_depth = filters["min_depth"]
        filter_expr.append(f"DP>={min_depth}")

    if filter_expr:
        cmd.extend(["-e", " || ".join(f"({expr})" for expr in filter_expr)])

    cmd.append(str(vcf_path))

    logger.info(f"Running bcftools filter: {' '.join(cmd[:5])}...")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        result.check_returncode()

        logger.info(f"bcftools filtering completed: {output_vcf}")
        return result

    except subprocess.CalledProcessError as e:
        logger.error(f"bcftools filtering failed: {e}")
        logger.error(f"stderr: {e.stderr}")
        raise


def call_variants_freebayes(
    bam_files: List[str | Path], reference_fasta: str | Path, output_vcf: str | Path, ploidy: int = 2
) -> subprocess.CompletedProcess:
    """Call variants using freebayes.

    Args:
        bam_files: List of BAM files
        reference_fasta: Reference genome FASTA
        output_vcf: Output VCF file
        ploidy: Sample ploidy

    Returns:
        CompletedProcess with calling results
    """
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    # Check if freebayes is available
    try:
        subprocess.run(["freebayes", "--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise FileNotFoundError("freebayes not found on PATH")

    cmd = [
        "freebayes",
        "-f",
        str(reference_fasta),
        "-p",
        str(ploidy),
        "--min-alternate-count",
        "2",
        "--min-alternate-fraction",
        "0.2",
        "--haplotype-length",
        "3",
    ]

    # Add BAM files
    for bam_file in bam_files:
        cmd.extend(["-b", str(bam_file)])

    logger.info(f"Running freebayes: {' '.join(cmd[:5])}...")

    try:
        # Write stdout to output file instead of using shell redirection
        with open(output_vcf, "w") as vcf_out:
            result = subprocess.run(cmd, stdout=vcf_out, stderr=subprocess.PIPE, text=True)
        result.check_returncode()

        logger.info(f"freebayes variant calling completed: {output_vcf}")
        return result

    except subprocess.CalledProcessError as e:
        logger.error(f"freebayes variant calling failed: {e}")
        logger.error(f"stderr: {e.stderr}")
        raise


def merge_vcfs(vcf_files: List[str | Path], output_vcf: str | Path) -> subprocess.CompletedProcess:
    """Merge multiple VCF files.

    Args:
        vcf_files: List of VCF files to merge
        output_vcf: Output merged VCF file

    Returns:
        CompletedProcess with merge results
    """
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    if len(vcf_files) < 2:
        raise ValueError("Need at least 2 VCF files to merge")

    cmd = ["bcftools", "merge", "-o", str(output_vcf), "-O", "v"]  # uncompressed VCF

    # Add input files
    for vcf_file in vcf_files:
        cmd.append(str(vcf_file))

    logger.info(f"Merging {len(vcf_files)} VCF files")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        result.check_returncode()

        logger.info(f"VCF merge completed: {output_vcf}")
        return result

    except subprocess.CalledProcessError as e:
        logger.error(f"VCF merge failed: {e}")
        logger.error(f"stderr: {e.stderr}")
        raise


def index_vcf(vcf_path: str | Path) -> subprocess.CompletedProcess:
    """Create tabix index for VCF file.

    Args:
        vcf_path: VCF file to index

    Returns:
        CompletedProcess with indexing results
    """
    vcf_path = Path(vcf_path)

    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    # Compress if not already compressed
    if vcf_path.suffix != ".gz":
        compressed_vcf = vcf_path.with_suffix(".vcf.gz")

        # Compress with bgzip (write stdout to file instead of shell redirection)
        compress_cmd = ["bgzip", "-c", str(vcf_path)]
        with open(compressed_vcf, "wb") as gz_out:
            subprocess.run(compress_cmd, stdout=gz_out, check=True)

        vcf_path = compressed_vcf

    # Index with tabix
    cmd = ["tabix", "-p", "vcf", str(vcf_path)]

    logger.info(f"Indexing VCF: {vcf_path}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        result.check_returncode()

        logger.info(f"VCF indexing completed: {vcf_path}.tbi")
        return result

    except subprocess.CalledProcessError as e:
        logger.error(f"VCF indexing failed: {e}")
        logger.error(f"stderr: {e.stderr}")
        raise


def validate_vcf(vcf_path: str | Path) -> Dict[str, Any]:
    """Validate VCF file format and content.

    Args:
        vcf_path: VCF file to validate

    Returns:
        Validation results dictionary
    """
    vcf_path = Path(vcf_path)

    if not vcf_path.exists():
        return {"valid": False, "error": "File not found"}

    validation = {
        "valid": True,
        "errors": [],
        "warnings": [],
        "stats": {"total_variants": 0, "samples": 0, "contigs": set()},
    }

    try:
        with open(vcf_path, "r") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                if not line:
                    continue

                if line.startswith("##"):
                    # Header line
                    continue
                elif line.startswith("#"):
                    # Column header
                    fields = line[1:].split("\t")
                    if len(fields) >= 8:  # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
                        validation["stats"]["samples"] = max(0, len(fields) - 9)  # Subtract fixed columns
                else:
                    # Variant line
                    fields = line.split("\t")
                    if len(fields) < 8:
                        validation["errors"].append(f"Line {line_num}: Insufficient fields")
                        validation["valid"] = False
                        continue

                    # Basic validation
                    chrom, pos, id_field, ref, alt, qual, filter_field, info = fields[:8]

                    # Validate position
                    try:
                        pos_int = int(pos)
                        if pos_int <= 0:
                            validation["errors"].append(f"Line {line_num}: Invalid position {pos}")
                            validation["valid"] = False
                    except ValueError:
                        validation["errors"].append(f"Line {line_num}: Non-numeric position {pos}")
                        validation["valid"] = False

                    # Track contigs
                    validation["stats"]["contigs"].add(chrom)

                    # Validate alleles
                    if not ref or ref == ".":
                        validation["warnings"].append(f"Line {line_num}: Missing reference allele")

                    if not alt or alt == ".":
                        validation["warnings"].append(f"Line {line_num}: No alternate alleles")

                    validation["stats"]["total_variants"] += 1

                    # Limit validation to first 1000 variants for performance
                    if validation["stats"]["total_variants"] >= 1000:
                        break

        validation["stats"]["contigs"] = list(validation["stats"]["contigs"])

    except Exception as e:
        validation["valid"] = False
        validation["errors"].append(f"Validation failed: {str(e)}")

    return validation


def check_bcftools_available() -> bool:
    """Check if bcftools is available on the system.

    Returns:
        True if bcftools is available and executable
    """
    try:
        result = subprocess.run(["bcftools", "--version"], capture_output=True, text=True, timeout=10)
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
        return False


def check_gatk_available() -> bool:
    """Check if GATK is available on the system.

    Returns:
        True if GATK is available and executable
    """
    try:
        # Try common GATK commands
        result = subprocess.run(["gatk", "--version"], capture_output=True, text=True, timeout=10)
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
        try:
            # Try java -jar gatk.jar
            result = subprocess.run(
                ["java", "-jar", "gatk.jar", "--version"],
                capture_output=True,
                text=True,
                timeout=10,
                cwd=".",  # Look in current directory
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
            return False


def merge_vcf_files(
    vcf_files: List[str | Path], output_vcf: str | Path, threads: int = 1
) -> subprocess.CompletedProcess:
    """Merge multiple VCF files using bcftools.

    Args:
        vcf_files: List of VCF files to merge
        output_vcf: Output merged VCF file
        threads: Number of threads to use

    Returns:
        CompletedProcess with bcftools merge results

    Raises:
        FileNotFoundError: If bcftools is not available or input files don't exist
        ValueError: If no input files provided
    """
    if not vcf_files:
        raise ValueError("Must provide at least one VCF file to merge")

    # Check that all input files exist
    for vcf_file in vcf_files:
        if not Path(vcf_file).exists():
            raise FileNotFoundError(f"VCF file not found: {vcf_file}")

    # Check bcftools availability
    if not check_bcftools_available():
        raise FileNotFoundError("bcftools not available for VCF merging")

    # Build bcftools merge command
    cmd = ["bcftools", "merge", "--threads", str(threads), "-o", str(output_vcf)] + [str(f) for f in vcf_files]

    logger.info(f"Merging {len(vcf_files)} VCF files with bcftools")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)  # 5 minute timeout for merging

        if result.returncode == 0:
            logger.info(f"Successfully merged VCF files to {output_vcf}")
        else:
            logger.error(f"VCF merging failed: {result.stderr}")

        return result

    except subprocess.TimeoutExpired:
        logger.error("VCF merging timed out")
        raise RuntimeError("VCF merging operation timed out")
    except Exception as e:
        logger.error(f"VCF merging failed with error: {e}")
        raise
