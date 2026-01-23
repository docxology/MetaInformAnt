"""DNA variant analysis and VCF file processing.

This module provides tools for analyzing genetic variants, processing VCF files,
and performing variant calling operations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)


def parse_vcf(path: str | Path) -> Dict[str, Any]:
    """Parse VCF (Variant Call Format) file.

    Args:
        path: Path to VCF file

    Returns:
        Dictionary containing VCF data and metadata

    Raises:
        FileNotFoundError: If VCF file doesn't exist
        ValueError: If VCF format is invalid

    Example:
        >>> # Assuming test.vcf exists
        >>> vcf_data = parse_vcf("test.vcf")
        >>> "variants" in vcf_data
        True
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"VCF file not found: {path}")

    variants = []
    metadata = {}
    samples = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith("##"):
                # Metadata lines
                if "=" in line:
                    key, value = line[2:].split("=", 1)
                    metadata[key] = value
                continue

            if line.startswith("#"):
                # Header line with sample names
                fields = line[1:].split("\t")
                if len(fields) >= 9:  # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, [samples...]
                    samples = fields[9:]
                continue

            # Variant line
            fields = line.split("\t")
            if len(fields) < 8:
                continue

            try:
                variant = {
                    "chrom": fields[0],
                    "pos": int(fields[1]),
                    "id": fields[2],
                    "ref": fields[3],
                    "alt": fields[4].split(",") if fields[4] != "." else [],
                    "qual": float(fields[5]) if fields[5] != "." else None,
                    "filter": fields[6].split(";") if fields[6] != "." else [],
                    "info": _parse_info_field(fields[7]),
                    "format": fields[8] if len(fields) > 8 else None,
                    "samples": {},
                }

                # Parse sample data
                if len(fields) > 9 and variant["format"]:
                    format_fields = variant["format"].split(":")
                    for i, sample_name in enumerate(samples):
                        if i + 9 < len(fields):
                            sample_data = fields[i + 9].split(":")
                            variant["samples"][sample_name] = dict(zip(format_fields, sample_data))

                variants.append(variant)

            except (ValueError, IndexError) as e:
                logger.warning(f"Skipping malformed variant line: {e}")
                continue

    return {"metadata": metadata, "samples": samples, "variants": variants, "total_variants": len(variants)}


def _parse_info_field(info_str: str) -> Dict[str, Any]:
    """Parse VCF INFO field."""
    if info_str == ".":
        return {}

    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            # Try to convert to appropriate type
            if value.isdigit():
                info[key] = int(value)
            elif value.replace(".", "").replace("-", "").isdigit():
                try:
                    info[key] = float(value)
                except ValueError:
                    info[key] = value
            else:
                info[key] = value
        else:
            info[item] = True

    return info


def filter_variants_by_quality(vcf_data: Dict[str, Any], min_qual: float = 20.0) -> Dict[str, Any]:
    """Filter variants by quality score.

    Args:
        vcf_data: Parsed VCF data from parse_vcf()
        min_qual: Minimum quality score

    Returns:
        Filtered VCF data dictionary
    """
    filtered_variants = []

    for variant in vcf_data["variants"]:
        if variant["qual"] is not None and variant["qual"] >= min_qual:
            filtered_variants.append(variant)

    result = vcf_data.copy()
    result["variants"] = filtered_variants
    result["total_variants"] = len(filtered_variants)

    logger.info(
        f"Filtered {len(vcf_data['variants'])} variants to {len(filtered_variants)} based on quality >= {min_qual}"
    )

    return result


def filter_variants_by_maf(vcf_data: Dict[str, Any], min_maf: float = 0.01) -> Dict[str, Any]:
    """Filter variants by minor allele frequency.

    Args:
        vcf_data: Parsed VCF data
        min_maf: Minimum minor allele frequency

    Returns:
        Filtered VCF data
    """
    filtered_variants = []

    for variant in vcf_data["variants"]:
        # Calculate MAF from sample data
        allele_counts = {}

        for sample_data in variant["samples"].values():
            if "GT" in sample_data:
                gt = sample_data["GT"]
                if gt not in [".", "./."]:
                    # Parse genotype (e.g., "0/1", "1/1")
                    alleles = []
                    for allele in gt.split("/"):
                        try:
                            alleles.append(int(allele))
                        except ValueError:
                            continue

                    for allele in alleles:
                        if allele > 0:  # Not reference
                            allele_counts[allele] = allele_counts.get(allele, 0) + 1

        # Calculate MAF
        total_alleles = len(variant["samples"]) * 2  # Diploid
        if total_alleles > 0:
            maf = min(allele_counts.values()) / total_alleles if allele_counts else 0.0
            if maf >= min_maf:
                filtered_variants.append(variant)

    result = vcf_data.copy()
    result["variants"] = filtered_variants
    result["total_variants"] = len(filtered_variants)

    logger.info(f"Filtered {len(vcf_data['variants'])} variants to {len(filtered_variants)} based on MAF >= {min_maf}")

    return result


def calculate_variant_statistics(vcf_data: Dict[str, Any]) -> Dict[str, Any]:
    """Calculate comprehensive variant statistics.

    Args:
        vcf_data: Parsed VCF data

    Returns:
        Statistics dictionary
    """
    variants = vcf_data["variants"]
    if not variants:
        return {"total_variants": 0}

    stats = {
        "total_variants": len(variants),
        "snps": 0,
        "indels": 0,
        "multiallelic": 0,
        "transitions": 0,
        "transversions": 0,
        "quality_distribution": {},
        "variant_types": {},
    }

    transition_pairs = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    for variant in variants:
        ref = variant["ref"].upper()
        alts = [alt.upper() for alt in variant["alt"]]

        # Count SNPs vs indels
        is_snp = all(len(ref) == len(alt) == 1 for alt in alts)
        if is_snp:
            stats["snps"] += 1

            # Count transitions/transversions
            for alt in alts:
                if len(alt) == 1 and (ref, alt) in transition_pairs:
                    stats["transitions"] += 1
                elif len(alt) == 1:
                    stats["transversions"] += 1
        else:
            stats["indels"] += 1

        # Count multiallelic variants
        if len(alts) > 1:
            stats["multiallelic"] += 1

        # Quality distribution
        if variant["qual"] is not None:
            qual_bucket = int(variant["qual"] // 10) * 10
            stats["quality_distribution"][qual_bucket] = stats["quality_distribution"].get(qual_bucket, 0) + 1

        # Variant type classification
        for alt in alts:
            if len(ref) == 1 and len(alt) == 1:
                var_type = "SNP"
            elif len(ref) > len(alt):
                var_type = "DEL"
            elif len(ref) < len(alt):
                var_type = "INS"
            else:
                var_type = "COMPLEX"

            stats["variant_types"][var_type] = stats["variant_types"].get(var_type, 0) + 1

    # Calculate ratios
    if stats["transitions"] + stats["transversions"] > 0:
        stats["transition_transversion_ratio"] = stats["transitions"] / stats["transversions"]
    else:
        stats["transition_transversion_ratio"] = 0.0

    stats["snp_ratio"] = stats["snps"] / len(variants) if variants else 0.0

    return stats


def annotate_variants(vcf_data: Dict[str, Any], annotations: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Add annotations to variants.

    Args:
        vcf_data: Parsed VCF data
        annotations: Dictionary mapping variant IDs to annotation data

    Returns:
        Annotated VCF data
    """
    annotated_variants = []

    for variant in vcf_data["variants"]:
        variant_id = variant["id"]
        if variant_id in annotations:
            # Add annotations to INFO field
            variant["info"].update(annotations[variant_id])

        annotated_variants.append(variant)

    result = vcf_data.copy()
    result["variants"] = annotated_variants

    logger.info(f"Annotated {len([v for v in annotated_variants if v['id'] in annotations])} variants")

    return result


def write_vcf(vcf_data: Dict[str, Any], output_path: str | Path) -> None:
    """Write VCF data to file.

    Args:
        vcf_data: VCF data dictionary from parse_vcf()
        output_path: Output file path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        # Write metadata
        for key, value in vcf_data["metadata"].items():
            f.write(f"##{key}={value}\n")

        # Write header
        header_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        if vcf_data["samples"]:
            header_fields.extend(["FORMAT"] + vcf_data["samples"])

        f.write("\t".join(header_fields) + "\n")

        # Write variants
        for variant in vcf_data["variants"]:
            fields = [
                variant["chrom"],
                str(variant["pos"]),
                variant["id"],
                variant["ref"],
                ",".join(variant["alt"]) if variant["alt"] else ".",
                str(variant["qual"]) if variant["qual"] is not None else ".",
                ";".join(variant["filter"]) if variant["filter"] else ".",
                _format_info_field(variant["info"]),
            ]

            # Add sample data if present
            if variant.get("format") and variant.get("samples"):
                fields.append(variant["format"])
                for sample in vcf_data["samples"]:
                    if sample in variant["samples"]:
                        sample_data = variant["samples"][sample]
                        if isinstance(sample_data, dict):
                            format_fields = variant["format"].split(":")
                            sample_values = [sample_data.get(field, ".") for field in format_fields]
                            fields.append(":".join(sample_values))
                        else:
                            fields.append(str(sample_data))
                    else:
                        fields.append(".")

            f.write("\t".join(fields) + "\n")

    logger.info(f"Wrote {len(vcf_data['variants'])} variants to {output_path}")


def _format_info_field(info: Dict[str, Any]) -> str:
    """Format INFO field for VCF output."""
    if not info:
        return "."

    items = []
    for key, value in info.items():
        if isinstance(value, bool) and value:
            items.append(key)
        else:
            items.append(f"{key}={value}")

    return ";".join(items)


def detect_variant_type(ref: str, alt: str) -> str:
    """Determine variant type from REF and ALT alleles.

    Args:
        ref: Reference allele
        alt: Alternative allele

    Returns:
        Variant type: 'SNP', 'INS', 'DEL', 'COMPLEX'

    Example:
        >>> detect_variant_type("A", "T")
        'SNP'
        >>> detect_variant_type("A", "AT")
        'INS'
    """
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    elif len(ref) > len(alt):
        return "DEL"
    elif len(ref) < len(alt):
        return "INS"
    else:
        return "COMPLEX"


def calculate_allele_frequency(vcf_data: Dict[str, Any], variant_index: int) -> Dict[str, float]:
    """Calculate allele frequencies for a specific variant.

    Args:
        vcf_data: Parsed VCF data
        variant_index: Index of variant in vcf_data['variants']

    Returns:
        Dictionary mapping alleles to frequencies
    """
    if variant_index >= len(vcf_data["variants"]):
        raise IndexError("Variant index out of range")

    variant = vcf_data["variants"][variant_index]
    allele_counts = {}

    # Count reference allele
    ref_allele = variant["ref"]
    allele_counts[ref_allele] = 0

    # Count alternative alleles
    for alt in variant["alt"]:
        allele_counts[alt] = 0

    # Count alleles in samples
    for sample_data in variant["samples"].values():
        if "GT" in sample_data:
            gt = sample_data["GT"]
            if gt not in [".", "./."]:
                alleles = gt.split("/")
                for allele_str in alleles:
                    try:
                        allele_idx = int(allele_str)
                        if allele_idx == 0:
                            allele = ref_allele
                        else:
                            allele = variant["alt"][allele_idx - 1]

                        allele_counts[allele] += 1
                    except (ValueError, IndexError):
                        continue

    # Convert to frequencies
    total_alleles = sum(allele_counts.values())
    frequencies = {}

    if total_alleles > 0:
        for allele, count in allele_counts.items():
            frequencies[allele] = count / total_alleles

    return frequencies
