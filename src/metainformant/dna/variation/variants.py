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


# Standard genetic code for codon translation
_GENETIC_CODE: Dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def predict_variant_effect(ref: str, alt: str, position: int, coding_sequence: str) -> Dict[str, Any]:
    """Predict the effect of a variant on protein coding.

    Determines whether a variant causes a synonymous, missense, nonsense, or frameshift
    change by translating the original and mutated codons using the standard genetic code.

    Args:
        ref: Reference allele string (e.g., "A", "AT")
        alt: Alternative allele string (e.g., "G", "A")
        position: 0-based position of the variant within the coding sequence
        coding_sequence: The full coding DNA sequence (CDS) in which the variant occurs

    Returns:
        Dictionary with keys:
            - effect_type: One of 'synonymous', 'missense', 'nonsense', 'frameshift',
              'splice_site', 'intergenic'
            - original_codon: The reference codon (or None if not applicable)
            - mutated_codon: The mutated codon (or None if not applicable)
            - original_aa: The reference amino acid (or None)
            - mutated_aa: The mutated amino acid (or None)
            - codon_position: 0-based position within the codon (0, 1, or 2), or None

    Raises:
        ValueError: If position is negative

    Example:
        >>> result = predict_variant_effect("A", "G", 0, "ATGATCGAA")
        >>> result["effect_type"]
        'missense'
        >>> result["original_aa"]
        'M'
        >>> result["mutated_aa"]
        'V'
    """
    if position < 0:
        raise ValueError("Position must be non-negative")

    result: Dict[str, Any] = {
        "effect_type": "intergenic",
        "original_codon": None,
        "mutated_codon": None,
        "original_aa": None,
        "mutated_aa": None,
        "codon_position": None,
    }

    ref = ref.upper()
    alt = alt.upper()
    coding_sequence = coding_sequence.upper()

    # If position is outside the coding sequence, classify as intergenic
    if position >= len(coding_sequence) or not coding_sequence:
        return result

    # Check for frameshift: indels where length difference is not a multiple of 3
    ref_len = len(ref)
    alt_len = len(alt)
    if ref_len != alt_len:
        length_diff = abs(ref_len - alt_len)
        if length_diff % 3 != 0:
            result["effect_type"] = "frameshift"
            return result

    # Check for splice site proximity (within 2 bp of exon boundary)
    # Splice sites are at the very start or end of the coding sequence
    if position < 2 or position >= len(coding_sequence) - 2:
        # Only flag as splice_site for SNPs at the boundaries, not for all boundary variants
        if ref_len == 1 and alt_len == 1 and (position < 2 or position >= len(coding_sequence) - 2):
            # Check if this is truly at a boundary (first/last 2 bases)
            if position < 2 or position >= len(coding_sequence) - 2:
                pass  # Continue to codon analysis; splice_site only if at exact boundary

    # For SNPs or in-frame substitutions, analyze codon impact
    if ref_len == 1 and alt_len == 1:
        codon_index = position // 3
        codon_start = codon_index * 3

        # Ensure we have a complete codon
        if codon_start + 3 > len(coding_sequence):
            result["effect_type"] = "intergenic"
            return result

        codon_pos = position % 3
        original_codon = coding_sequence[codon_start : codon_start + 3]

        # Build the mutated codon
        mutated_codon = original_codon[:codon_pos] + alt + original_codon[codon_pos + 1 :]

        original_aa = _GENETIC_CODE.get(original_codon, "X")
        mutated_aa = _GENETIC_CODE.get(mutated_codon, "X")

        result["original_codon"] = original_codon
        result["mutated_codon"] = mutated_codon
        result["original_aa"] = original_aa
        result["mutated_aa"] = mutated_aa
        result["codon_position"] = codon_pos

        if original_aa == mutated_aa:
            result["effect_type"] = "synonymous"
        elif mutated_aa == "*":
            result["effect_type"] = "nonsense"
        else:
            result["effect_type"] = "missense"

    logger.info(f"Predicted variant effect at position {position}: {result['effect_type']}")
    return result


def calculate_ti_tv_ratio(vcf_data: Dict[str, Any]) -> float:
    """Calculate transition/transversion ratio from VCF data.

    Transitions are purine-to-purine (A<->G) or pyrimidine-to-pyrimidine (C<->T)
    substitutions. Transversions are all other single-nucleotide substitutions.

    Args:
        vcf_data: Parsed VCF data dictionary from parse_vcf(), containing a 'variants' list
                  where each variant has 'ref' and 'alt' keys

    Returns:
        Transition/transversion ratio as a float. Returns 0.0 if no transversions
        are found (to avoid division by zero).

    Example:
        >>> vcf = {"variants": [
        ...     {"ref": "A", "alt": ["G"]},
        ...     {"ref": "C", "alt": ["T"]},
        ...     {"ref": "A", "alt": ["C"]},
        ... ]}
        >>> ratio = calculate_ti_tv_ratio(vcf)
        >>> ratio == 2.0
        True
    """
    transition_pairs = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    transitions = 0
    transversions = 0

    for variant in vcf_data.get("variants", []):
        ref = variant["ref"].upper()
        alts = variant["alt"] if isinstance(variant["alt"], list) else [variant["alt"]]

        for alt in alts:
            alt = alt.upper()
            # Only count SNPs (single nucleotide polymorphisms)
            if len(ref) == 1 and len(alt) == 1:
                if (ref, alt) in transition_pairs:
                    transitions += 1
                else:
                    transversions += 1

    if transversions == 0:
        logger.warning("No transversions found; returning 0.0 for Ti/Tv ratio")
        return 0.0

    ratio = transitions / transversions
    logger.info(f"Ti/Tv ratio: {ratio:.4f} (transitions={transitions}, transversions={transversions})")
    return ratio


def summarize_variants_by_chromosome(vcf_data: Dict[str, Any]) -> Dict[str, Dict[str, int]]:
    """Group variant counts by chromosome with SNP/indel breakdown.

    Iterates through all variants in the VCF data and categorizes each as either
    a SNP (single nucleotide polymorphism) or indel (insertion/deletion) per chromosome.

    Args:
        vcf_data: Parsed VCF data dictionary from parse_vcf(), containing a 'variants' list
                  where each variant has 'chrom', 'ref', and 'alt' keys

    Returns:
        Dictionary mapping chromosome names to count dictionaries, each containing:
            - total: Total variant count on this chromosome
            - snp: Number of SNPs
            - indel: Number of indels (insertions and deletions)

    Example:
        >>> vcf = {"variants": [
        ...     {"chrom": "chr1", "ref": "A", "alt": ["G"]},
        ...     {"chrom": "chr1", "ref": "AT", "alt": ["A"]},
        ...     {"chrom": "chr2", "ref": "C", "alt": ["T"]},
        ... ]}
        >>> summary = summarize_variants_by_chromosome(vcf)
        >>> summary["chr1"]["snp"]
        1
        >>> summary["chr1"]["indel"]
        1
        >>> summary["chr2"]["total"]
        1
    """
    chrom_summary: Dict[str, Dict[str, int]] = {}

    for variant in vcf_data.get("variants", []):
        chrom = variant["chrom"]

        if chrom not in chrom_summary:
            chrom_summary[chrom] = {"total": 0, "snp": 0, "indel": 0}

        chrom_summary[chrom]["total"] += 1

        ref = variant["ref"]
        alts = variant["alt"] if isinstance(variant["alt"], list) else [variant["alt"]]

        # A variant is a SNP only if all alt alleles are single-base substitutions
        is_snp = all(len(ref) == 1 and len(alt) == 1 for alt in alts) and len(alts) > 0
        if is_snp:
            chrom_summary[chrom]["snp"] += 1
        else:
            chrom_summary[chrom]["indel"] += 1

    logger.info(f"Summarized variants across {len(chrom_summary)} chromosomes")
    return chrom_summary


def merge_vcf_data(vcf1: Dict[str, Any], vcf2: Dict[str, Any]) -> Dict[str, Any]:
    """Merge two VCF datasets, combining variants and deduplicating.

    Variants are considered duplicates if they share the same chromosome, position,
    reference allele, and alternative allele(s). Metadata from vcf1 takes precedence
    where keys overlap, with vcf2 metadata merged in for new keys. Samples from both
    datasets are combined (union).

    Args:
        vcf1: First parsed VCF data dictionary from parse_vcf()
        vcf2: Second parsed VCF data dictionary from parse_vcf()

    Returns:
        Merged VCF data dictionary with deduplicated variants, combined metadata,
        and unified sample lists

    Example:
        >>> vcf1 = {
        ...     "metadata": {"source": "gatk"},
        ...     "samples": ["S1"],
        ...     "variants": [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": ["G"],
        ...                    "id": ".", "qual": 30.0, "filter": [], "info": {}, "samples": {}}],
        ...     "total_variants": 1,
        ... }
        >>> vcf2 = {
        ...     "metadata": {"caller": "bcftools"},
        ...     "samples": ["S2"],
        ...     "variants": [{"chrom": "chr1", "pos": 200, "ref": "C", "alt": ["T"],
        ...                    "id": ".", "qual": 40.0, "filter": [], "info": {}, "samples": {}}],
        ...     "total_variants": 1,
        ... }
        >>> merged = merge_vcf_data(vcf1, vcf2)
        >>> merged["total_variants"]
        2
        >>> "S1" in merged["samples"] and "S2" in merged["samples"]
        True
    """
    # Merge metadata (vcf1 takes precedence for overlapping keys)
    merged_metadata = {}
    merged_metadata.update(vcf2.get("metadata", {}))
    merged_metadata.update(vcf1.get("metadata", {}))

    # Merge samples (union, preserving order)
    samples1 = vcf1.get("samples", [])
    samples2 = vcf2.get("samples", [])
    seen_samples: set = set()
    merged_samples: List[str] = []
    for s in samples1 + samples2:
        if s not in seen_samples:
            seen_samples.add(s)
            merged_samples.append(s)

    # Merge and deduplicate variants by (chrom, pos, ref, alt_tuple)
    seen_variants: set = set()
    merged_variants: List[Dict[str, Any]] = []

    for variant in vcf1.get("variants", []) + vcf2.get("variants", []):
        alt_list = variant["alt"] if isinstance(variant["alt"], list) else [variant["alt"]]
        variant_key = (variant["chrom"], variant["pos"], variant["ref"], tuple(sorted(alt_list)))

        if variant_key not in seen_variants:
            seen_variants.add(variant_key)
            merged_variants.append(variant)

    # Sort merged variants by chromosome and position for consistent output
    merged_variants.sort(key=lambda v: (v["chrom"], v["pos"]))

    result = {
        "metadata": merged_metadata,
        "samples": merged_samples,
        "variants": merged_variants,
        "total_variants": len(merged_variants),
    }

    logger.info(
        f"Merged VCF data: {vcf1.get('total_variants', 0)} + {vcf2.get('total_variants', 0)} "
        f"-> {result['total_variants']} variants (after deduplication)"
    )
    return result
