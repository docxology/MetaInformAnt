"""GWAS quality control and data filtering utilities.

This module provides functions for quality control of GWAS data,
including VCF parsing, genotype filtering, and data validation.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)


def parse_vcf_full(vcf_path: Union[str, Path]) -> Dict[str, Any]:
    """Parse a complete VCF file.

    Args:
        vcf_path: Path to VCF file

    Returns:
        Dictionary containing parsed VCF data

    Raises:
        FileNotFoundError: If VCF file doesn't exist
        ValueError: If VCF format is invalid
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    logger.info(f"Parsing VCF file: {vcf_path}")

    data = {
        'samples': [],
        'variants': [],
        'genotypes': [],  # Will be numpy array in full implementation
        'metadata': {}
    }

    try:
        # Open file (handle gzipped files)
        if str(vcf_path).endswith('.gz'):
            opener = gzip.open
            mode = 'rt'
        else:
            opener = open
            mode = 'r'

        with opener(vcf_path, mode) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                if not line:
                    continue

                if line.startswith('##'):
                    # Metadata lines
                    if '=' in line:
                        key, value = line[2:].split('=', 1)
                        data['metadata'][key] = value
                    continue

                if line.startswith('#'):
                    # Header line with sample names
                    parts = line[1:].split('\t')
                    if len(parts) > 9:  # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, then samples
                        data['samples'] = parts[9:]
                    continue

                # Variant line
                parts = line.split('\t')
                if len(parts) < 10:
                    continue

                try:
                    chrom = parts[0]
                    pos = int(parts[1])
                    var_id = parts[2]
                    ref = parts[3]
                    alt = parts[4]
                    qual = float(parts[5]) if parts[5] != '.' else None
                    filter_field = parts[6]
                    info = parts[7]
                    format_field = parts[8]

                    variant = {
                        'chrom': chrom,
                        'pos': pos,
                        'id': var_id,
                        'ref': ref,
                        'alt': alt.split(','),
                        'qual': qual,
                        'filter': filter_field,
                        'info': info,
                        'format': format_field
                    }

                    data['variants'].append(variant)

                    # Parse genotypes for each sample
                    genotype_row = []
                    for sample_gt in parts[9:]:
                        # Simplified genotype parsing
                        gt_value = 0  # Default to reference
                        if sample_gt and sample_gt != './.':
                            # Extract first allele from GT field
                            gt_parts = sample_gt.split(':')[0].split('/')
                            if len(gt_parts) >= 1 and gt_parts[0].isdigit():
                                gt_value = int(gt_parts[0])

                        genotype_row.append(gt_value)

                    data['genotypes'].append(genotype_row)

                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing variant at line {line_num}: {e}")
                    continue

    except Exception as e:
        raise ValueError(f"Error reading VCF file {vcf_path}: {e}")

    logger.info(f"Parsed {len(data['variants'])} variants for {len(data['samples'])} samples")

    return data


def filter_by_maf(genotypes: List[List[int]], maf_threshold: float = 0.01) -> Tuple[List[List[int]], List[int]]:
    """Filter variants by minor allele frequency.

    Args:
        genotypes: Genotype matrix (variants x samples)
        maf_threshold: Minimum minor allele frequency

    Returns:
        Tuple of (filtered_genotypes, kept_indices)
    """
    if not genotypes:
        return [], []

    kept_indices = []

    for i, variant_genotypes in enumerate(genotypes):
        # Calculate allele frequencies
        allele_counts = {}
        total_alleles = 0

        for genotype in variant_genotypes:
            if genotype == 0:  # Homozygous reference
                allele_counts[0] = allele_counts.get(0, 0) + 2
                total_alleles += 2
            elif genotype == 1:  # Heterozygous
                allele_counts[0] = allele_counts.get(0, 0) + 1
                allele_counts[1] = allele_counts.get(1, 0) + 1
                total_alleles += 2
            elif genotype == 2:  # Homozygous alternate
                allele_counts[1] = allele_counts.get(1, 0) + 2
                total_alleles += 2

        if total_alleles == 0:
            continue

        # Calculate MAF
        allele_freqs = {allele: count / total_alleles for allele, count in allele_counts.items()}
        maf = min(allele_freqs.values()) if len(allele_freqs) > 1 else 0.0

        if maf >= maf_threshold:
            kept_indices.append(i)

    # Extract kept genotypes
    filtered_genotypes = [genotypes[i] for i in kept_indices]

    logger.info(f"MAF filter: kept {len(kept_indices)}/{len(genotypes)} variants")

    return filtered_genotypes, kept_indices


def filter_by_missing(genotypes: List[List[int]], missing_threshold: float = 0.1) -> Tuple[List[List[int]], List[int]]:
    """Filter variants by missing data proportion.

    Args:
        genotypes: Genotype matrix (variants x samples)
        missing_threshold: Maximum proportion of missing data

    Returns:
        Tuple of (filtered_genotypes, kept_indices)
    """
    if not genotypes:
        return [], []

    kept_indices = []

    for i, variant_genotypes in enumerate(genotypes):
        # Count missing genotypes (represented as -1 or similar)
        missing_count = sum(1 for gt in variant_genotypes if gt < 0)
        missing_prop = missing_count / len(variant_genotypes)

        if missing_prop <= missing_threshold:
            kept_indices.append(i)

    # Extract kept genotypes
    filtered_genotypes = [genotypes[i] for i in kept_indices]

    logger.info(f"Missing data filter: kept {len(kept_indices)}/{len(genotypes)} variants")

    return filtered_genotypes, kept_indices


def test_hwe(genotypes: List[List[int]], alpha: float = 0.05) -> List[float]:
    """Test variants for Hardy-Weinberg equilibrium.

    Args:
        genotypes: Genotype matrix (variants x samples)
        alpha: Significance threshold

    Returns:
        List of p-values for HWE test
    """
    p_values = []

    for variant_genotypes in genotypes:
        # Count genotypes
        aa_count = sum(1 for gt in variant_genotypes if gt == 0)  # Homozygous reference
        ab_count = sum(1 for gt in variant_genotypes if gt == 1)  # Heterozygous
        bb_count = sum(1 for gt in variant_genotypes if gt == 2)  # Homozygous alternate

        total = aa_count + ab_count + bb_count
        if total == 0:
            p_values.append(1.0)
            continue

        # Calculate allele frequencies
        p = (2 * aa_count + ab_count) / (2 * total)  # Frequency of reference allele
        q = 1 - p  # Frequency of alternate allele

        # Expected genotype frequencies under HWE
        expected_aa = p * p * total
        expected_ab = 2 * p * q * total
        expected_bb = q * q * total

        # Chi-square test
        try:
            chi_square = (
                ((aa_count - expected_aa) ** 2) / expected_aa if expected_aa > 0 else 0 +
                ((ab_count - expected_ab) ** 2) / expected_ab if expected_ab > 0 else 0 +
                ((bb_count - expected_bb) ** 2) / expected_bb if expected_bb > 0 else 0
            )

            # Simplified p-value calculation (would use scipy.stats.chi2.sf in practice)
            # For now, use a basic approximation
            if chi_square > 5.99:  # Chi-square critical value for df=2, alpha=0.05
                p_value = 0.05
            else:
                p_value = 0.95

            p_values.append(p_value)

        except ZeroDivisionError:
            p_values.append(1.0)

    return p_values


def apply_qc_filters(vcf_data: Dict[str, Any], min_maf: float = 0.01,
                    max_missing: float = 0.1, min_hwe_p: float = 1e-6) -> Dict[str, Any]:
    """Apply comprehensive quality control filters to VCF data.

    Args:
        vcf_data: Parsed VCF data dictionary
        min_maf: Minimum minor allele frequency
        max_missing: Maximum missing data proportion
        min_hwe_p: Minimum HWE p-value threshold

    Returns:
        Filtered VCF data dictionary
    """
    logger.info("Applying QC filters to VCF data")

    genotypes = vcf_data.get('genotypes', [])
    if not genotypes:
        logger.warning("No genotype data found in VCF")
        return vcf_data

    # Apply filters
    filtered_genotypes, maf_indices = filter_by_maf(genotypes, min_maf)
    if not filtered_genotypes:
        logger.warning("No variants passed MAF filter")
        return {**vcf_data, 'genotypes': [], 'variants': []}

    # Re-filter by missing data on MAF-filtered set
    filtered_genotypes, missing_indices = filter_by_missing(filtered_genotypes, max_missing)
    if not filtered_genotypes:
        logger.warning("No variants passed missing data filter")
        return {**vcf_data, 'genotypes': [], 'variants': []}

    # Map back to original indices
    final_indices = [maf_indices[i] for i in missing_indices]

    # Apply HWE filter
    hwe_p_values = test_hwe(filtered_genotypes, alpha=min_hwe_p)
    hwe_kept = [i for i, p_val in enumerate(hwe_p_values) if p_val >= min_hwe_p]
    final_genotypes = [filtered_genotypes[i] for i in hwe_kept]
    final_indices = [final_indices[i] for i in hwe_kept]

    # Filter variants accordingly
    filtered_variants = [vcf_data['variants'][i] for i in final_indices]

    logger.info(f"QC filtering complete: {len(final_genotypes)}/{len(genotypes)} variants passed")

    return {
        **vcf_data,
        'genotypes': final_genotypes,
        'variants': filtered_variants,
        'qc_filters': {
            'min_maf': min_maf,
            'max_missing': max_missing,
            'min_hwe_p': min_hwe_p,
            'variants_before': len(genotypes),
            'variants_after': len(final_genotypes)
        }
    }


def write_filtered_vcf(filtered_data: Dict[str, Any], output_path: Union[str, Path]) -> None:
    """Write filtered VCF data to a new file.

    Args:
        filtered_data: Filtered VCF data dictionary
        output_path: Output VCF file path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Writing filtered VCF to {output_path}")

    # This is a placeholder implementation
    # Full implementation would reconstruct VCF format from filtered data
    with open(output_path, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for sample in filtered_data.get('samples', []):
            f.write(f"\t{sample}")
        f.write("\n")

        # Write filtered variants (simplified)
        for i, variant in enumerate(filtered_data.get('variants', [])):
            f.write(f"{variant['chrom']}\t{variant['pos']}\t{variant['id']}\t{variant['ref']}\t{','.join(variant['alt'])}\t")
            f.write(f"{variant.get('qual', '.')}\t{variant.get('filter', 'PASS')}\t{variant.get('info', '.')}\t{variant.get('format', 'GT')}")
            # Write genotypes (simplified)
            genotypes = filtered_data.get('genotypes', [[]] * len(filtered_data.get('variants', [])))
            if i < len(genotypes):
                for gt in genotypes[i]:
                    f.write(f"\t{gt}|{gt}")
            f.write("\n")

    logger.info(f"Filtered VCF written with {len(filtered_data.get('variants', []))} variants")


def extract_variant_regions(vcf_data: Dict[str, Any], regions: List[Tuple[str, int, int]],
                           padding: int = 0) -> Dict[str, Any]:
    """Extract variants within specified genomic regions.

    Args:
        vcf_data: Parsed VCF data dictionary
        regions: List of (chrom, start, end) tuples
        padding: Additional padding around regions (bp)

    Returns:
        Dictionary with variants filtered to specified regions
    """
    if not regions:
        logger.warning("No regions specified, returning all data")
        return vcf_data.copy()

    variants = vcf_data.get('variants', [])
    if not variants:
        logger.warning("No variants in VCF data")
        return vcf_data.copy()

    # Filter variants by regions
    filtered_variants = []
    filtered_indices = []

    for i, variant in enumerate(variants):
        chrom = variant.get('chrom', '')
        pos = variant.get('pos', 0)

        # Check if variant falls within any region
        for region_chrom, region_start, region_end in regions:
            if (chrom == region_chrom and
                region_start - padding <= pos <= region_end + padding):
                filtered_variants.append(variant)
                filtered_indices.append(i)
                break

    # Filter genotypes if present
    filtered_genotypes = []
    genotypes = vcf_data.get('genotypes', [])
    if genotypes and len(genotypes) == len(variants):
        filtered_genotypes = [genotypes[i] for i in filtered_indices]

    # Create result dictionary
    result = vcf_data.copy()
    result['variants'] = filtered_variants
    if filtered_genotypes:
        result['genotypes'] = filtered_genotypes

    logger.info(f"Extracted {len(filtered_variants)} variants from {len(regions)} regions")

    return result

