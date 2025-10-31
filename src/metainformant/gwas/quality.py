"""Variant quality control and filtering."""

from __future__ import annotations

import logging
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from ..core.io import ensure_directory
from ..dna.population import allele_frequencies, observed_heterozygosity
from ..math.popgen import hardy_weinberg_genotype_freqs

logger = logging.getLogger(__name__)

# Try importing scipy for statistical tests
try:
    from scipy.stats import chi2

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    chi2 = None


def parse_vcf_full(path: str | Path) -> dict[str, Any]:
    """Parse VCF file and extract full genotype matrix.

    Extends dna.variants.parse_vcf to extract genotypes per sample.

    Args:
        path: Path to VCF file (supports .gz)

    Returns:
        Dictionary with:
        - samples: List of sample IDs
        - variants: List of variant records (CHROM, POS, ID, REF, ALT, QUAL, FILTER)
        - genotypes: 2D array (samples x variants) with genotypes encoded as:
          0 = homozygous REF (0/0)
          1 = heterozygous (0/1 or 1/0)
          2 = homozygous ALT (1/1)
          -1 = missing (./.)
    """
    logger.info(f"parse_vcf_full: Parsing VCF file {path}")

    path_obj = Path(path)
    if not path_obj.exists():
        raise FileNotFoundError(f"VCF file not found: {path_obj}")

    # Handle gzipped files
    if path_obj.suffix == ".gz":
        import gzip

        open_func = gzip.open
        mode = "rt"
    else:
        open_func = open
        mode = "r"

    samples: list[str] = []
    variants: list[dict[str, Any]] = []
    genotypes: list[list[int]] = []  # Will be transposed later

    with open_func(path_obj, mode) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                parts = line.split("\t")
                samples = parts[9:] if len(parts) > 9 else []
                continue

            # Parse variant line
            parts = line.split("\t")
            if len(parts) < 9:
                continue

            chrom = parts[0]
            pos = parts[1]
            var_id = parts[2]
            ref = parts[3]
            alt = parts[4]
            qual = parts[5]
            filt = parts[6]
            info = parts[7]
            format_field = parts[8]

            variant = {
                "CHROM": chrom,
                "POS": int(pos) if pos.isdigit() else 0,
                "ID": var_id,
                "REF": ref,
                "ALT": alt,
                "QUAL": qual,
                "FILTER": filt,
                "INFO": info,
                "FORMAT": format_field,
            }

            # Parse genotypes for this variant
            variant_genotypes: list[int] = []
            for sample_idx in range(9, len(parts)):
                sample_genotype = parts[sample_idx] if sample_idx < len(parts) else "./."
                # Parse GT field (assuming GT is first field)
                gt_part = sample_genotype.split(":")[0] if ":" in sample_genotype else sample_genotype

                if gt_part in ["./.", ".|.", "."]:
                    variant_genotypes.append(-1)  # Missing
                else:
                    try:
                        alleles = gt_part.split("/") if "/" in gt_part else gt_part.split("|")
                        if len(alleles) >= 2:
                            a1, a2 = int(alleles[0]), int(alleles[1])
                            if a1 == -1 or a2 == -1:
                                variant_genotypes.append(-1)  # Missing
                            else:
                                # Encode as 0=homREF, 1=het, 2=homALT
                                genotype_code = a1 + a2
                                variant_genotypes.append(genotype_code)
                        else:
                            variant_genotypes.append(-1)
                    except (ValueError, IndexError):
                        variant_genotypes.append(-1)  # Missing

            genotypes.append(variant_genotypes)
            variants.append(variant)

    # Transpose genotypes: samples x variants
    if genotypes:
        num_samples = len(samples)
        num_variants = len(genotypes)
        genotype_matrix: list[list[int]] = [[-1] * num_variants for _ in range(num_samples)]

        for var_idx, var_genotypes in enumerate(genotypes):
            for sample_idx, genotype in enumerate(var_genotypes):
                if sample_idx < num_samples:
                    genotype_matrix[sample_idx][var_idx] = genotype
    else:
        genotype_matrix = []

    logger.info(f"parse_vcf_full: Parsed {len(variants)} variants across {len(samples)} samples")

    return {
        "samples": samples,
        "variants": variants,
        "genotypes": genotype_matrix,  # samples x variants
    }


def filter_by_maf(
    genotypes: list[list[int]],
    maf_threshold: float = 0.01,
) -> tuple[list[int], float]:
    """Filter variants by minor allele frequency.

    Args:
        genotypes: Genotype matrix (samples x variants), encoded as 0/1/2/-1
        maf_threshold: Minimum minor allele frequency

    Returns:
        Tuple of (variant_indices_to_keep, maf_value)
        variant_indices_to_keep: List of variant indices passing MAF filter
        maf_value: The calculated MAF (for logging)
    """
    if not genotypes or not genotypes[0]:
        return ([], 0.0)

    num_samples = len(genotypes)
    num_variants = len(genotypes[0])
    passing_indices: list[int] = []

    for var_idx in range(num_variants):
        # Extract genotypes for this variant
        var_genotypes = [genotypes[sample_idx][var_idx] for sample_idx in range(num_samples)]

        # Calculate allele counts (exclude missing -1)
        allele_count = 0
        total_alleles = 0

        for gt in var_genotypes:
            if gt == -1:  # Missing
                continue
            # Genotype encoding: 0=0/0, 1=0/1 or 1/0, 2=1/1
            # Count ALT alleles (each genotype contributes 0, 1, or 2 ALT alleles)
            allele_count += gt
            total_alleles += 2

        if total_alleles == 0:
            continue  # All missing

        # Calculate MAF
        alt_freq = allele_count / total_alleles
        maf = min(alt_freq, 1.0 - alt_freq)

        if maf >= maf_threshold:
            passing_indices.append(var_idx)

    return (passing_indices, maf)


def filter_by_missing(
    genotypes: list[list[int]],
    max_missing: float = 0.05,
) -> list[int]:
    """Filter variants by missing data rate.

    Args:
        genotypes: Genotype matrix (samples x variants), encoded as 0/1/2/-1
        max_missing: Maximum allowed missing genotype rate

    Returns:
        List of variant indices passing missing data filter
    """
    if not genotypes or not genotypes[0]:
        return []

    num_samples = len(genotypes)
    num_variants = len(genotypes[0])
    passing_indices: list[int] = []

    for var_idx in range(num_variants):
        var_genotypes = [genotypes[sample_idx][var_idx] for sample_idx in range(num_samples)]
        missing_count = sum(1 for gt in var_genotypes if gt == -1)
        missing_rate = missing_count / num_samples if num_samples > 0 else 1.0

        if missing_rate <= max_missing:
            passing_indices.append(var_idx)

    return passing_indices


def test_hwe(
    genotypes: list[int],
) -> tuple[float, float]:
    """Test Hardy-Weinberg equilibrium for a variant.

    Args:
        genotypes: List of genotypes for one variant (0=homREF, 1=het, 2=homALT, -1=missing)

    Returns:
        Tuple of (chi2_statistic, p_value)
    """
    # Remove missing genotypes
    valid_genotypes = [gt for gt in genotypes if gt != -1]
    if len(valid_genotypes) < 3:
        return (0.0, 1.0)  # Not enough data

    # Count genotypes
    n_hom_ref = sum(1 for gt in valid_genotypes if gt == 0)
    n_het = sum(1 for gt in valid_genotypes if gt == 1)
    n_hom_alt = sum(1 for gt in valid_genotypes if gt == 2)
    n_total = len(valid_genotypes)

    # Calculate allele frequencies
    n_ref_alleles = 2 * n_hom_ref + n_het
    n_alt_alleles = 2 * n_hom_alt + n_het
    n_total_alleles = n_ref_alleles + n_alt_alleles

    if n_total_alleles == 0:
        return (0.0, 1.0)

    p = n_ref_alleles / n_total_alleles  # Allele A frequency
    q = 1.0 - p  # Allele a frequency

    # Expected genotype frequencies under HWE
    expected_hom_ref = n_total * p * p
    expected_het = n_total * 2 * p * q
    expected_hom_alt = n_total * q * q

    # Chi-square test
    chi2_stat = 0.0
    if expected_hom_ref > 0:
        chi2_stat += ((n_hom_ref - expected_hom_ref) ** 2) / expected_hom_ref
    if expected_het > 0:
        chi2_stat += ((n_het - expected_het) ** 2) / expected_het
    if expected_hom_alt > 0:
        chi2_stat += ((n_hom_alt - expected_hom_alt) ** 2) / expected_hom_alt

    # Calculate p-value (1 degree of freedom)
    if SCIPY_AVAILABLE and chi2:
        p_value = 1.0 - chi2.cdf(chi2_stat, df=1)
    else:
        # Simple approximation for chi2 p-value with df=1
        # P(X > x) ≈ exp(-x/2) for large x
        if chi2_stat > 0:
            p_value = math.exp(-chi2_stat / 2.0)
        else:
            p_value = 1.0

    return (chi2_stat, p_value)


def apply_qc_filters(
    vcf_path: str | Path,
    config: dict[str, Any],
    output_vcf: str | Path | None = None,
) -> dict[str, Any]:
    """Apply all quality control filters to VCF file.

    Args:
        vcf_path: Path to input VCF file
        config: QC configuration dictionary with thresholds
        output_vcf: Optional path to output filtered VCF

    Returns:
        Dictionary with QC results and statistics
    """
    logger.info(f"apply_qc_filters: Applying QC filters to {vcf_path}")

    min_maf = config.get("min_maf", 0.01)
    max_missing = config.get("max_missing", 0.05)
    min_call_rate = config.get("min_call_rate", 0.95)
    hwe_pval = config.get("hwe_pval", 1e-6)
    exclude_indels = config.get("exclude_indels", True)
    min_qual = config.get("min_qual", 30.0)

    # Parse VCF
    vcf_data = parse_vcf_full(vcf_path)
    samples = vcf_data["samples"]
    variants = vcf_data["variants"]
    genotypes = vcf_data["genotypes"]

    if not variants:
        return {
            "status": "failed",
            "error": "No variants found in VCF",
        }

    num_variants_before = len(variants)
    logger.info(f"apply_qc_filters: Starting with {num_variants_before} variants")

    # Filter by quality score
    quality_passing: list[int] = []
    for var_idx, variant in enumerate(variants):
        qual_str = variant.get("QUAL", ".")
        if qual_str != ".":
            try:
                qual_val = float(qual_str)
                if qual_val >= min_qual:
                    quality_passing.append(var_idx)
            except ValueError:
                # Missing quality, include
                quality_passing.append(var_idx)
        else:
            # Missing quality, include by default
            quality_passing.append(var_idx)

    # Filter by indels
    indel_passing: list[int] = []
    if exclude_indels:
        for var_idx, variant in enumerate(variants):
            ref = variant.get("REF", "")
            alt = variant.get("ALT", "")
            # Simple check: if REF or ALT length > 1, it's an indel
            if len(ref) == 1 and len(alt.split(",")[0]) == 1:  # Only check first ALT
                indel_passing.append(var_idx)
    else:
        indel_passing = list(range(len(variants)))

    # Filter by missing data
    missing_passing = filter_by_missing(genotypes, max_missing=max_missing)

    # Filter by MAF
    maf_passing_indices, _ = filter_by_maf(genotypes, maf_threshold=min_maf)

    # Filter by HWE
    hwe_passing: list[int] = []
    # Ensure hwe_pval is numeric
    if isinstance(hwe_pval, str):
        hwe_pval = float(hwe_pval)
    for var_idx in range(len(variants)):
        var_genotypes = [genotypes[sample_idx][var_idx] for sample_idx in range(len(samples))]
        _, p_value = test_hwe(var_genotypes)
        if isinstance(p_value, (int, float)) and not np.isnan(p_value) and p_value >= hwe_pval:
            hwe_passing.append(var_idx)

    # Combine all filters (intersection)
    all_passing = set(range(len(variants)))
    all_passing &= set(quality_passing)
    all_passing &= set(indel_passing)
    all_passing &= set(missing_passing)
    all_passing &= set(maf_passing_indices)
    all_passing &= set(hwe_passing)

    passing_indices = sorted(list(all_passing))
    num_variants_after = len(passing_indices)

    result = {
        "status": "success",
        "input_vcf": str(Path(vcf_path)),
        "num_variants_before": num_variants_before,
        "num_variants_after": num_variants_after,
        "num_variants_removed": num_variants_before - num_variants_after,
        "filters_applied": {
            "min_maf": min_maf,
            "max_missing": max_missing,
            "min_call_rate": min_call_rate,
            "hwe_pval": hwe_pval,
            "exclude_indels": exclude_indels,
            "min_qual": min_qual,
        },
        "filter_counts": {
            "quality_passing": len(quality_passing),
            "indel_passing": len(indel_passing),
            "missing_passing": len(missing_passing),
            "maf_passing": len(maf_passing_indices),
            "hwe_passing": len(hwe_passing),
            "all_passing": len(passing_indices),
        },
        "passing_variant_indices": passing_indices,
    }

    logger.info(
        f"apply_qc_filters: QC complete - {num_variants_before} -> {num_variants_after} variants "
        f"({num_variants_before - num_variants_after} removed)"
    )

    # Write filtered VCF if output path provided
    if output_vcf:
        logger.info(f"apply_qc_filters: Writing filtered VCF to {output_vcf}")
        # TODO: Implement VCF writing (requires proper VCF format handling)
        result["output_vcf"] = str(Path(output_vcf))
        result["message"] = "Filtered VCF writing not yet fully implemented"

    return result

