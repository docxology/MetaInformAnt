"""Apis mellifera (Amel_HAv3.1) genome mapping and annotation utilities.

This module provides chromosome name normalization, genome coordinate mapping,
and GFF3 gene annotation parsing for the Apis mellifera reference genome.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Amel_HAv3.1 NCBI RefSeq chromosome accessions -> chromosome number
# Reference: GCF_003254395.2_Amel_HAv3.1
AMEL_HAV3_CHROMOSOMES: Dict[str, int] = {
    "NC_037638.1": 1,
    "NC_037639.1": 2,
    "NC_037640.1": 3,
    "NC_037641.1": 4,
    "NC_037642.1": 5,
    "NC_037643.1": 6,
    "NC_037644.1": 7,
    "NC_037645.1": 8,
    "NC_037646.1": 9,
    "NC_037647.1": 10,
    "NC_037648.1": 11,
    "NC_037649.1": 12,
    "NC_037650.1": 13,
    "NC_037651.1": 14,
    "NC_037652.1": 15,
    "NC_037653.1": 16,
}

# Reverse mapping: chromosome number -> NCBI accession
AMEL_HAV3_ACCESSIONS: Dict[int, str] = {v: k for k, v in AMEL_HAV3_CHROMOSOMES.items()}

# Approximate chromosome sizes in bp (Amel_HAv3.1)
AMEL_HAV3_CHROM_SIZES: Dict[int, int] = {
    1: 27754200,
    2: 16089512,
    3: 13619445,
    4: 11340508,
    5: 13543907,
    6: 17789102,
    7: 12717210,
    8: 12354651,
    9: 12360052,
    10: 12718924,
    11: 16352600,
    12: 11514234,
    13: 9534514,
    14: 10021392,
    15: 9946621,
    16: 7238532,
}


def normalize_chromosome_name(chrom: str) -> Optional[int]:
    """Normalize various chromosome name formats to integer chromosome number.

    Handles:
    - NCBI accessions: NC_037638.1 -> 1
    - Group names: Group1, GroupI -> 1
    - chr prefix: chr1, chrI -> 1
    - LG prefix: LG1, LGI -> 1
    - Bare integers: 1 -> 1
    - Roman numerals: I, II, III, ... XVI -> 1-16

    Args:
        chrom: Chromosome name in any supported format

    Returns:
        Integer chromosome number (1-16) or None if not recognized
    """
    chrom = str(chrom).strip()

    # NCBI accession
    if chrom in AMEL_HAV3_CHROMOSOMES:
        return AMEL_HAV3_CHROMOSOMES[chrom]

    # Strip common prefixes
    for prefix in ["Group", "group", "chr", "Chr", "CHR", "LG", "lg", "Lg"]:
        if chrom.startswith(prefix):
            chrom = chrom[len(prefix) :]
            break

    # Try bare integer
    try:
        num = int(chrom)
        if 1 <= num <= 16:
            return num
        return None
    except ValueError:
        pass

    # Try Roman numerals
    roman_map = {
        "I": 1,
        "II": 2,
        "III": 3,
        "IV": 4,
        "V": 5,
        "VI": 6,
        "VII": 7,
        "VIII": 8,
        "IX": 9,
        "X": 10,
        "XI": 11,
        "XII": 12,
        "XIII": 13,
        "XIV": 14,
        "XV": 15,
        "XVI": 16,
    }
    upper = chrom.upper()
    if upper in roman_map:
        return roman_map[upper]

    return None


def get_chromosome_order() -> List[int]:
    """Return the standard chromosome order for Apis mellifera.

    Returns:
        List of chromosome numbers [1, 2, ..., 16]
    """
    return list(range(1, 17))


def get_genome_size() -> int:
    """Return the total genome size in bp.

    Returns:
        Total genome size
    """
    return sum(AMEL_HAV3_CHROM_SIZES.values())


def parse_gff3_genes(
    gff_path: Union[str, Path],
    feature_type: str = "gene",
) -> List[Dict[str, Any]]:
    """Parse gene locations from an Amel_HAv3.1 GFF3 annotation file.

    Args:
        gff_path: Path to GFF3 annotation file
        feature_type: GFF3 feature type to extract (default: "gene")

    Returns:
        List of gene dictionaries with keys:
        - chrom: chromosome number (int)
        - start: start position (1-based)
        - end: end position (1-based)
        - strand: '+' or '-'
        - gene_id: gene identifier
        - gene_name: gene name/symbol (if available)
        - biotype: gene biotype (if available)

    Raises:
        FileNotFoundError: If GFF3 file doesn't exist
    """
    gff_path = Path(gff_path)
    if not gff_path.exists():
        raise FileNotFoundError(f"GFF3 file not found: {gff_path}")

    genes: List[Dict[str, Any]] = []
    logger.info(f"Parsing GFF3 file: {gff_path}")

    import gzip

    opener = gzip.open if str(gff_path).endswith(".gz") else open
    mode = "rt" if str(gff_path).endswith(".gz") else "r"

    with opener(gff_path, mode) as f:  # type: ignore[call-overload]
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            if parts[2] != feature_type:
                continue

            chrom_num = normalize_chromosome_name(parts[0])
            if chrom_num is None:
                continue

            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue

            strand = parts[6] if parts[6] in ("+", "-") else "."

            # Parse attributes
            attrs = _parse_gff3_attributes(parts[8])

            gene = {
                "chrom": chrom_num,
                "start": start,
                "end": end,
                "strand": strand,
                "gene_id": attrs.get("ID", attrs.get("gene_id", "")),
                "gene_name": attrs.get("Name", attrs.get("gene", "")),
                "biotype": attrs.get("gene_biotype", attrs.get("biotype", "")),
            }
            genes.append(gene)

    logger.info(f"Parsed {len(genes)} {feature_type} features from GFF3")
    return genes


def _parse_gff3_attributes(attr_string: str) -> Dict[str, str]:
    """Parse GFF3 attribute string into dictionary.

    Args:
        attr_string: GFF3 column 9 attribute string (key=value;key=value;...)

    Returns:
        Dictionary of attribute key-value pairs
    """
    attrs: Dict[str, str] = {}
    for pair in attr_string.split(";"):
        pair = pair.strip()
        if "=" in pair:
            key, value = pair.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs
