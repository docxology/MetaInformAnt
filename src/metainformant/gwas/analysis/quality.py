"""GWAS quality control and data filtering utilities.

This module provides functions for quality control of GWAS data,
including VCF parsing, genotype filtering, and data validation.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def parse_vcf_full(vcf_path: Union[str, Path]) -> Dict[str, Any]:
    """Parse entire VCF file into a dictionary (NumPy optimized for large cohorts).

    A sidecar ``.npz`` cache is used when it matches the source file size and
    modification time. The cache stores JSON metadata and NumPy arrays only, so
    cache reads do not require pickle.

    Args:
        vcf_path: Path to VCF file

    Returns:
        Dictionary with samples, variants, genotypes (numpy array), metadata
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    logger.info(f"Parsing VCF file (NumPy optimized): {vcf_path}")

    cache_path = vcf_path.with_suffix(vcf_path.suffix + ".npz")
    source_stat = vcf_path.stat()
    if cache_path.exists():
        logger.info(f"Found cached genotype matrix at {cache_path}. Loading natively...")
        try:
            with np.load(cache_path, allow_pickle=False) as cached:
                cache_mtime_ns = int(cached["source_mtime_ns"][()])
                cache_size = int(cached["source_size"][()])
                if cache_mtime_ns != source_stat.st_mtime_ns or cache_size != source_stat.st_size:
                    raise ValueError("cache source metadata does not match VCF")

                data = {
                    "samples": cached["samples"].tolist(),
                    "variants": json.loads(str(cached["variants_json"][()])),
                    "metadata": json.loads(str(cached["metadata_json"][()])),
                    "genotypes": cached["genotypes"],
                }
            logger.info("Successfully loaded VCF cache via optimized protocol.")
            return data
        except Exception as e:
            logger.warning(f"Failed to load cache from {cache_path}: {e}. Falling back to parsing.")

    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    mode = "rt" if str(vcf_path).endswith(".gz") else "r"

    data = {
        "samples": [],
        "variants": [],
        "metadata": {},
    }

    n_variants = 0
    try:
        with opener(vcf_path, mode) as f:
            for line in f:
                if not line.startswith("#"):
                    n_variants += 1
                    if n_variants % 2500000 == 0:
                        logger.info(f"Pre-scan: Processed {n_variants:,} variants...")
                elif line.startswith("##"):
                    if "=" in line:
                        key, val = line.strip()[2:].split("=", 1)
                        data["metadata"][key] = val
                elif line.startswith("#CHROM"):
                    data["samples"] = line.strip().split("\t")[9:]
    except Exception as e:
        raise ValueError(f"Error reading VCF file {vcf_path}: {e}")

    logger.info(f"Pre-scan detected {n_variants} variants")

    n_samples = len(data["samples"])
    genotypes = np.full((n_variants, n_samples), -1, dtype=np.int8)

    try:
        with opener(vcf_path, mode) as f:
            var_idx = 0
            for line_num, line in enumerate(f, 1):
                if line.startswith("#"):
                    continue

                if var_idx > 0 and var_idx % 1000000 == 0:
                    logger.info(
                        "Loading variants into NumPy matrix: %s/%s (%.1f%%)",
                        f"{var_idx:,}",
                        f"{n_variants:,}",
                        (var_idx / n_variants) * 100,
                    )

                parts = line.strip().split("\t")
                if len(parts) < 10:
                    continue

                try:
                    chrom = parts[0]
                    pos = int(parts[1])
                    var_id = parts[2]
                    ref = parts[3]
                    alt = parts[4]
                    qual = float(parts[5]) if parts[5] != "." else None

                    variant = {
                        "chrom": chrom,
                        "pos": pos,
                        "id": var_id,
                        "ref": ref,
                        "alt": alt.split(","),
                        "qual": qual,
                        "filter": parts[6],
                        "info": parts[7],
                        "format": parts[8],
                    }
                    data["variants"].append(variant)

                    for s_idx, sample_gt in enumerate(parts[9:]):
                        gt_field = sample_gt.split(":")[0] if sample_gt else "./."
                        if gt_field not in ("./.", ".|.", ".", ""):
                            sep = "|" if "|" in gt_field else "/"
                            gt_parts = gt_field.split(sep)
                            try:
                                genotypes[var_idx, s_idx] = sum(int(a) for a in gt_parts if a.isdigit())
                            except (ValueError, TypeError):
                                pass
                    var_idx += 1
                except (ValueError, IndexError) as e:
                    logger.warning(f"Error parsing variant at line {line_num}: {e}")
                    continue
    except Exception as e:
        raise ValueError(f"Error reading genotypes from {vcf_path}: {e}")

    logger.info(f"Successfully loaded {n_variants} variants for {n_samples} samples into NumPy matrix")

    # Transpose genotypes: from [variant][sample] to [sample][variant]
    data["genotypes"] = genotypes.T

    try:
        np.savez_compressed(
            cache_path,
            genotypes=data["genotypes"],
            samples=np.asarray(data["samples"], dtype=str),
            variants_json=json.dumps(data["variants"]),
            metadata_json=json.dumps(data["metadata"]),
            source_mtime_ns=np.asarray(source_stat.st_mtime_ns, dtype=np.int64),
            source_size=np.asarray(source_stat.st_size, dtype=np.int64),
        )
        logger.info(f"Cached genotype matrix to {cache_path}")
    except Exception as e:
        logger.warning(f"Failed to cache generated VCF arrays: {e}")

    return data


def filter_by_maf(genotypes: List[List[int]], maf_threshold: float = 0.01) -> Tuple[List[List[int]], List[int]]:
    """Filter variants by minor allele frequency.

    Args:
        genotypes: Genotype matrix (variants x samples)
        maf_threshold: Minimum minor allele frequency

    Returns:
        Tuple of (filtered_genotypes, kept_indices)
    """
    if len(genotypes) == 0:
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
    if len(genotypes) == 0:
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
    """Test variants for Hardy-Weinberg equilibrium using chi-square test.

    Args:
        genotypes: Genotype matrix (variants x samples)
        alpha: Significance threshold (unused, kept for API compatibility)

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

        # Chi-square test with 1 degree of freedom
        chi_square = 0.0
        if expected_aa > 0:
            chi_square += ((aa_count - expected_aa) ** 2) / expected_aa
        if expected_ab > 0:
            chi_square += ((ab_count - expected_ab) ** 2) / expected_ab
        if expected_bb > 0:
            chi_square += ((bb_count - expected_bb) ** 2) / expected_bb

        # Calculate p-value using chi-square distribution with df=1
        # Using incomplete gamma function approximation
        p_value = _chi2_sf(chi_square, df=1)
        p_values.append(p_value)

    return p_values


def _chi2_sf(x: float, df: int = 1) -> float:
    """Calculate survival function (1 - CDF) for chi-square distribution.

    Uses incomplete gamma function approximation for the chi-square distribution.

    Args:
        x: Chi-square statistic
        df: Degrees of freedom

    Returns:
        p-value (probability of observing chi-square >= x under null)
    """
    import math

    if x <= 0:
        return 1.0

    # For df=1, we can use the relationship with standard normal
    if df == 1:
        # chi2 with df=1 is square of standard normal
        # P(chi2 > x) = 2 * P(Z > sqrt(x)) = 2 * (1 - Phi(sqrt(x))) = erfc(sqrt(x/2))
        z = math.sqrt(x / 2)
        return math.erfc(z)

    # For other df, use incomplete gamma function approximation
    # P(chi2 > x) = 1 - P(chi2 <= x) = 1 - gammainc(df/2, x/2)
    # Using Pearson's approximation for larger df
    a = df / 2
    z = x / 2

    # Series approximation for regularized incomplete gamma function
    if z < a + 1:
        # Use series expansion
        ap = a
        sum_val = 1.0 / a
        delta = sum_val
        for _ in range(100):
            ap += 1
            delta *= z / ap
            sum_val += delta
            if abs(delta) < abs(sum_val) * 1e-10:
                break
        return 1.0 - sum_val * math.exp(-z + a * math.log(z) - math.lgamma(a))
    else:
        # Use continued fraction
        b = z + 1 - a
        c = 1.0 / 1e-30
        d = 1.0 / b
        h = d
        for i in range(1, 100):
            an = -i * (i - a)
            b += 2
            d = an * d + b
            if abs(d) < 1e-30:
                d = 1e-30
            c = b + an / c
            if abs(c) < 1e-30:
                c = 1e-30
            d = 1.0 / d
            delta = d * c
            h *= delta
            if abs(delta - 1.0) < 1e-10:
                break
        return h * math.exp(-z + a * math.log(z) - math.lgamma(a))


def check_haplodiploidy(
    vcf_data: Dict[str, Any],
    het_threshold: float = 0.05,
) -> Dict[str, Any]:
    """Detect haploid samples in a haplodiploid species (e.g., Hymenoptera).

    In haplodiploid species like Apis mellifera, males (drones) are haploid and
    will show near-zero heterozygosity. This function identifies such samples.

    Args:
        vcf_data: Parsed VCF data dictionary with 'genotypes' (samples x variants)
        het_threshold: Maximum heterozygosity rate to classify as haploid (default: 0.05)

    Returns:
        Dictionary with:
        - haploid_samples: list of sample indices identified as haploid
        - diploid_samples: list of sample indices identified as diploid
        - het_rates: heterozygosity rate per sample
        - n_haploid: count of haploid samples
        - n_diploid: count of diploid samples
    """
    genotypes = vcf_data.get("genotypes", [])
    if len(genotypes) == 0:
        return {
            "haploid_samples": [],
            "diploid_samples": [],
            "het_rates": [],
            "n_haploid": 0,
            "n_diploid": 0,
        }

    n_samples = len(genotypes)
    het_rates = []
    haploid_samples = []
    diploid_samples = []

    for s in range(n_samples):
        sample_gts = genotypes[s]
        n_het = sum(1 for g in sample_gts if g == 1)
        n_valid = sum(1 for g in sample_gts if g >= 0)
        het_rate = n_het / n_valid if n_valid > 0 else 0.0
        het_rates.append(het_rate)

        if het_rate <= het_threshold:
            haploid_samples.append(s)
        else:
            diploid_samples.append(s)

    logger.info(
        f"Haplodiploidy check: {len(haploid_samples)} haploid, "
        f"{len(diploid_samples)} diploid samples (threshold={het_threshold})"
    )

    return {
        "haploid_samples": haploid_samples,
        "diploid_samples": diploid_samples,
        "het_rates": het_rates,
        "n_haploid": len(haploid_samples),
        "n_diploid": len(diploid_samples),
    }


def apply_qc_filters(
    vcf_input: Union[str, Path, Dict[str, Any]],
    qc_config: Optional[Dict[str, Any]] = None,  # noqa: F821
) -> Dict[str, Any]:
    """Apply comprehensive quality control filters to VCF data.

    Args:
        vcf_input: VCF file path or pre-parsed VCF data dictionary
        qc_config: QC configuration dictionary with keys:
            - min_maf: Minimum minor allele frequency (default 0.01)
            - max_missing: Maximum missing data proportion (default 0.1)
            - min_hwe_p: Minimum HWE p-value threshold (default 1e-6)
            - min_qual: Minimum variant quality (optional)

    Returns:
        Dictionary with status and QC results
    """
    logger.info("Applying QC filters to VCF data")

    # Parse VCF if path is provided
    if isinstance(vcf_input, (str, Path)):
        vcf_data = parse_vcf_full(vcf_input)
    else:
        vcf_data = vcf_input

    # Extract QC parameters from config
    if qc_config is None:
        qc_config = {}
    min_maf = qc_config.get("min_maf", 0.01)
    max_missing = qc_config.get("max_missing", 0.1)
    min_hwe_p = qc_config.get("min_hwe_p", 1e-6)

    # Get genotypes - need to transpose from [sample][variant] to [variant][sample] for filtering
    genotypes_by_sample = vcf_data.get("genotypes", [])
    if len(genotypes_by_sample) == 0:
        logger.warning("No genotype data found in VCF")
        return {
            "status": "success",
            "num_variants_before": 0,
            "num_variants_after": 0,
            "filtered_data": vcf_data,
        }

    # Transpose genotypes to [variant][sample] for filtering functions
    n_samples = len(genotypes_by_sample)
    n_variants = len(genotypes_by_sample[0]) if n_samples > 0 else 0
    genotypes = [[genotypes_by_sample[s][v] for s in range(n_samples)] for v in range(n_variants)]

    num_variants_before = len(genotypes)

    # Apply filters
    filtered_genotypes, maf_indices = filter_by_maf(genotypes, min_maf)
    if len(filtered_genotypes) == 0:
        logger.warning("No variants passed MAF filter")
        return {
            "status": "success",
            "num_variants_before": num_variants_before,
            "num_variants_after": 0,
            "filtered_data": {**vcf_data, "genotypes": [], "variants": []},
        }

    # Re-filter by missing data on MAF-filtered set
    filtered_genotypes, missing_indices = filter_by_missing(filtered_genotypes, max_missing)
    if len(filtered_genotypes) == 0:
        logger.warning("No variants passed missing data filter")
        return {
            "status": "success",
            "num_variants_before": num_variants_before,
            "num_variants_after": 0,
            "filtered_data": {**vcf_data, "genotypes": [], "variants": []},
        }

    # Map back to original indices
    final_indices = [maf_indices[i] for i in missing_indices]

    # Apply HWE filter
    hwe_p_values = test_hwe(filtered_genotypes, alpha=min_hwe_p)
    hwe_kept = [i for i, p_val in enumerate(hwe_p_values) if p_val >= min_hwe_p]
    final_genotypes = [filtered_genotypes[i] for i in hwe_kept]
    final_indices = [final_indices[i] for i in hwe_kept]

    # Filter variants accordingly
    filtered_variants = [vcf_data["variants"][i] for i in final_indices]

    # Transpose final genotypes back to [sample][variant]
    n_final_variants = len(final_genotypes)
    final_genotypes_by_sample = [[final_genotypes[v][s] for v in range(n_final_variants)] for s in range(n_samples)]

    logger.info(f"QC filtering complete: {len(final_genotypes)}/{num_variants_before} variants passed")

    return {
        "status": "success",
        "num_variants_before": num_variants_before,
        "num_variants_after": len(final_genotypes),
        "filtered_data": {
            **vcf_data,
            "genotypes": final_genotypes_by_sample,
            "variants": filtered_variants,
        },
        "qc_filters": {
            "min_maf": min_maf,
            "max_missing": max_missing,
            "min_hwe_p": min_hwe_p,
        },
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

    with open(output_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")

        # Write metadata headers from original VCF data
        metadata = filtered_data.get("metadata", {})
        for key, value in metadata.items():
            # Skip fileformat since we already wrote it
            if key.lower() == "fileformat":
                continue
            f.write(f"##{key}={value}\n")

        # Write standard VCF header definitions if not already present
        if "INFO" not in metadata:
            f.write('##INFO=<ID=.,Number=.,Type=String,Description="Variant information">\n')
        if "FORMAT" not in metadata:
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        if "FILTER" not in metadata:
            f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')

        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for sample in filtered_data.get("samples", []):
            f.write(f"\t{sample}")
        f.write("\n")

        # Write filtered variants (simplified)
        for i, variant in enumerate(filtered_data.get("variants", [])):
            f.write(
                f"{variant['chrom']}\t{variant['pos']}\t{variant['id']}\t{variant['ref']}\t{','.join(variant['alt'])}\t"
            )
            f.write(
                f"{variant.get('qual', '.')}\t{variant.get('filter', 'PASS')}\t{variant.get('info', '.')}\t{variant.get('format', 'GT')}"
            )
            # Write genotypes (convert numeric to VCF GT format)
            variant_count = len(filtered_data.get("variants", []))
            all_genotypes = filtered_data.get("genotypes", [])
            if not all_genotypes:
                all_genotypes = [[] for _ in range(variant_count)]
            if i < len(all_genotypes):
                for gt in all_genotypes[i]:
                    # Convert numeric genotype to VCF GT format
                    if gt < 0:
                        f.write("\t./.")  # Missing genotype
                    elif gt == 0:
                        f.write("\t0/0")  # Homozygous reference
                    elif gt == 1:
                        f.write("\t0/1")  # Heterozygous
                    elif gt == 2:
                        f.write("\t1/1")  # Homozygous alternate
                    else:
                        f.write(f"\t{gt}/{gt}")  # Generic case
            f.write("\n")

    logger.info(f"Filtered VCF written with {len(filtered_data.get('variants', []))} variants")


def subset_vcf_data(
    vcf_data: Dict[str, Any],
    sample_ids: List[str] | None = None,
    sample_list_file: Union[str, Path, None] = None,
    subset_config: Dict[str, Any] | None = None,
    metadata: Dict[str, Dict[str, Any]] | None = None,
) -> Dict[str, Any]:
    """Subset parsed VCF data to specific samples.

    Priority: sample_ids > sample_list_file > subset_config (with metadata).
    Preserves genotype alignment (samples x variants format).
    Returns a new dict — does not modify the original.

    Args:
        vcf_data: Parsed VCF data dictionary from parse_vcf_full().
        sample_ids: Explicit list of sample IDs to keep.
        sample_list_file: Path to file with one sample ID per line.
        subset_config: Dict with keys like 'subspecies', 'caste',
            'max_per_subspecies' for metadata-based filtering.
        metadata: Dict mapping sample_id -> {subspecies, caste, ...}.

    Returns:
        New VCF data dict with only the selected samples.
    """
    all_samples = vcf_data.get("samples", [])
    if not all_samples:
        logger.warning("No samples in VCF data, returning as-is")
        return dict(vcf_data)

    # Determine target sample IDs
    target_ids: List[str] | None = None

    if sample_ids is not None:
        target_ids = sample_ids
    elif sample_list_file is not None:
        sample_list_path = Path(sample_list_file)
        if sample_list_path.exists():
            with open(sample_list_path) as f:
                target_ids = [line.strip() for line in f if line.strip()]
        else:
            logger.warning(f"Sample list file not found: {sample_list_file}")
    elif subset_config and metadata:
        # Metadata-based filtering
        target_ids = list(all_samples)  # Start with all

        allowed_subspecies = subset_config.get("subspecies")
        allowed_caste = subset_config.get("caste")
        max_per_sub = subset_config.get("max_per_subspecies")

        if allowed_subspecies or allowed_caste:
            filtered = []
            for sid in target_ids:
                meta = metadata.get(sid, {})
                if allowed_subspecies:
                    sample_sub = meta.get("subspecies", meta.get("population", ""))
                    if sample_sub not in allowed_subspecies:
                        continue
                if allowed_caste:
                    sample_caste = meta.get("caste", "")
                    if sample_caste not in allowed_caste:
                        continue
                filtered.append(sid)
            target_ids = filtered

        if max_per_sub and metadata:
            # Cap per subspecies
            from collections import defaultdict

            counts: Dict[str, int] = defaultdict(int)
            capped = []
            for sid in target_ids:
                meta = metadata.get(sid, {})
                sub = meta.get("subspecies", meta.get("population", "unknown"))
                if counts[sub] < max_per_sub:
                    capped.append(sid)
                    counts[sub] += 1
            target_ids = capped

    if target_ids is None:
        logger.info("No subsetting criteria specified, returning all samples")
        return dict(vcf_data)

    # Find indices of target samples in the VCF
    sample_to_idx = {s: i for i, s in enumerate(all_samples)}
    keep_indices = []
    kept_ids = []
    for sid in target_ids:
        idx = sample_to_idx.get(sid)
        if idx is not None:
            keep_indices.append(idx)
            kept_ids.append(sid)
        else:
            logger.warning(f"Sample '{sid}' not found in VCF, skipping")

    if not keep_indices:
        logger.warning("No matching samples found after subsetting")
        return {**vcf_data, "samples": [], "genotypes": []}

    # Subset genotypes (samples x variants format)
    genotypes = vcf_data.get("genotypes", [])
    subset_genotypes = [genotypes[i] for i in keep_indices] if genotypes is not None and len(genotypes) > 0 else []

    logger.info(f"Subset VCF: {len(all_samples)} -> {len(kept_ids)} samples")

    return {
        **vcf_data,
        "samples": kept_ids,
        "genotypes": subset_genotypes,
    }


def extract_variant_regions(
    vcf_data: Dict[str, Any], regions: List[Tuple[str, int, int]], padding: int = 0
) -> Dict[str, Any]:
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

    variants = vcf_data.get("variants", [])
    if not variants:
        logger.warning("No variants in VCF data")
        return vcf_data.copy()

    # Filter variants by regions
    filtered_variants = []
    filtered_indices = []

    for i, variant in enumerate(variants):
        chrom = variant.get("chrom", "")
        pos = variant.get("pos", 0)

        # Check if variant falls within any region
        for region_chrom, region_start, region_end in regions:
            if chrom == region_chrom and region_start - padding <= pos <= region_end + padding:
                filtered_variants.append(variant)
                filtered_indices.append(i)
                break

    # Filter genotypes if present
    filtered_genotypes = []
    genotypes = vcf_data.get("genotypes", [])
    if len(genotypes) > 0 and len(genotypes) == len(variants):
        filtered_genotypes = [genotypes[i] for i in filtered_indices]

    # Create result dictionary
    result = vcf_data.copy()
    result["variants"] = filtered_variants
    if len(filtered_genotypes) > 0:
        result["genotypes"] = filtered_genotypes

    logger.info(f"Extracted {len(filtered_variants)} variants from {len(regions)} regions")

    return result
