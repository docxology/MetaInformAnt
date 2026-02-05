"""SNP-to-gene annotation for GWAS results.

This module provides functions for annotating GWAS variants with nearest genes
and classifying variant locations relative to gene features.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)


def annotate_variants_with_genes(
    results: List[Dict[str, Any]],
    gff_path: Union[str, Path],
    window_kb: int = 50,
) -> List[Dict[str, Any]]:
    """Annotate GWAS results with nearest gene information.

    For each variant, finds the nearest gene within a window and classifies
    the variant's location (intragenic, upstream, downstream, intergenic).

    Args:
        results: List of GWAS result dictionaries (must have 'chrom' and 'pos')
        gff_path: Path to GFF3 annotation file
        window_kb: Window size in kilobases for gene search (default: 50kb)

    Returns:
        Annotated results (same list with added 'annotation' field per result)
    """
    from metainformant.gwas.data.genome import normalize_chromosome_name, parse_gff3_genes

    # Parse genes from GFF3
    genes = parse_gff3_genes(gff_path)

    # Organize genes by chromosome for fast lookup
    genes_by_chrom: Dict[int, List[Dict[str, Any]]] = {}
    for gene in genes:
        chrom = gene["chrom"]
        if chrom not in genes_by_chrom:
            genes_by_chrom[chrom] = []
        genes_by_chrom[chrom].append(gene)

    # Sort genes by start position for binary search
    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort(key=lambda g: g["start"])

    window_bp = window_kb * 1000
    annotated = 0

    for result in results:
        chrom_raw = result.get("chrom", "")
        pos = result.get("pos", 0)

        chrom = normalize_chromosome_name(str(chrom_raw))
        if chrom is None or pos <= 0:
            result["annotation"] = {
                "location": "unknown",
                "nearest_gene": None,
                "distance": None,
            }
            continue

        chrom_genes = genes_by_chrom.get(chrom, [])
        annotation = classify_variant_location(chrom, pos, chrom_genes, window_bp)
        result["annotation"] = annotation
        annotated += 1

    logger.info(f"Annotated {annotated}/{len(results)} variants with gene information")
    return results


def classify_variant_location(
    chrom: int,
    pos: int,
    genes: List[Dict[str, Any]],
    window_bp: int = 50000,
) -> Dict[str, Any]:
    """Classify a variant's location relative to nearby genes.

    Args:
        chrom: Chromosome number
        pos: Variant position (1-based)
        genes: List of gene dictionaries (sorted by start) on this chromosome
        window_bp: Search window in base pairs

    Returns:
        Dictionary with:
        - location: 'intragenic', 'upstream', 'downstream', or 'intergenic'
        - nearest_gene: gene_id of nearest gene (or None)
        - nearest_gene_name: gene name/symbol (or None)
        - distance: distance to nearest gene (0 if intragenic)
        - gene_strand: strand of nearest gene
    """
    if not genes:
        return {
            "location": "intergenic",
            "nearest_gene": None,
            "nearest_gene_name": None,
            "distance": None,
            "gene_strand": None,
        }

    nearest_gene = None
    nearest_distance = float("inf")
    location = "intergenic"

    for gene in genes:
        gene_start = gene["start"]
        gene_end = gene["end"]
        gene_strand = gene.get("strand", ".")

        # Check if variant is within the gene
        if gene_start <= pos <= gene_end:
            return {
                "location": "intragenic",
                "nearest_gene": gene.get("gene_id", ""),
                "nearest_gene_name": gene.get("gene_name", ""),
                "distance": 0,
                "gene_strand": gene_strand,
            }

        # Calculate distance to gene
        if pos < gene_start:
            dist = gene_start - pos
            # Relative to gene strand
            if gene_strand == "+":
                rel_location = "upstream"
            elif gene_strand == "-":
                rel_location = "downstream"
            else:
                rel_location = "upstream"
        else:
            dist = pos - gene_end
            if gene_strand == "+":
                rel_location = "downstream"
            elif gene_strand == "-":
                rel_location = "upstream"
            else:
                rel_location = "downstream"

        if dist < nearest_distance:
            nearest_distance = dist
            nearest_gene = gene
            location = rel_location

    # Check if within window
    if nearest_gene is None or nearest_distance > window_bp:
        return {
            "location": "intergenic",
            "nearest_gene": nearest_gene.get("gene_id", "") if nearest_gene else None,
            "nearest_gene_name": nearest_gene.get("gene_name", "") if nearest_gene else None,
            "distance": int(nearest_distance) if nearest_gene else None,
            "gene_strand": nearest_gene.get("strand", None) if nearest_gene else None,
        }

    return {
        "location": location,
        "nearest_gene": nearest_gene.get("gene_id", ""),
        "nearest_gene_name": nearest_gene.get("gene_name", ""),
        "distance": int(nearest_distance),
        "gene_strand": nearest_gene.get("strand", "."),
    }
