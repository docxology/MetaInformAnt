"""Gene and regulatory element overlap annotation for structural variants.

Provides interval-based overlap detection between structural variants and
genomic features (genes, regulatory elements), nearest-gene finding, and
reciprocal overlap computation.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


@dataclass
class GenomicInterval:
    """A genomic interval.

    Attributes:
        chrom: Chromosome name.
        start: Start position (0-based, inclusive).
        end: End position (0-based, exclusive).
        name: Feature name (e.g., gene name).
        feature_type: Type of feature (e.g., 'gene', 'exon', 'promoter').
        strand: Strand ('+', '-', or '.').
        info: Additional annotation fields.
    """

    chrom: str
    start: int
    end: int
    name: str = ""
    feature_type: str = ""
    strand: str = "."
    info: dict[str, Any] = field(default_factory=dict)


@dataclass
class OverlapResult:
    """Result of an overlap analysis between a variant and a genomic feature.

    Attributes:
        variant_chrom: Chromosome of the variant.
        variant_start: Start of the variant.
        variant_end: End of the variant.
        feature: The overlapping genomic feature.
        overlap_bp: Number of overlapping base pairs.
        overlap_fraction_variant: Fraction of the variant overlapping the feature.
        overlap_fraction_feature: Fraction of the feature overlapping the variant.
        relationship: Relationship type ('overlap', 'contained_in', 'contains', 'upstream', 'downstream').
    """

    variant_chrom: str
    variant_start: int
    variant_end: int
    feature: GenomicInterval
    overlap_bp: int = 0
    overlap_fraction_variant: float = 0.0
    overlap_fraction_feature: float = 0.0
    relationship: str = "overlap"


class IntervalIndex:
    """Efficient interval index for overlap queries.

    Uses a sorted-start approach with binary search for fast overlap queries.
    Suitable for moderate-size feature databases (up to millions of intervals).
    """

    def __init__(self, intervals: list[GenomicInterval]) -> None:
        """Build an interval index.

        Args:
            intervals: List of GenomicInterval objects to index.
        """
        # Group by chromosome
        self._by_chrom: dict[str, list[GenomicInterval]] = {}
        for iv in intervals:
            if iv.chrom not in self._by_chrom:
                self._by_chrom[iv.chrom] = []
            self._by_chrom[iv.chrom].append(iv)

        # Sort each chromosome's intervals by start position
        for chrom in self._by_chrom:
            self._by_chrom[chrom].sort(key=lambda x: (x.start, x.end))

        # Build start-position arrays for binary search
        self._starts: dict[str, list[int]] = {}
        for chrom, ivs in self._by_chrom.items():
            self._starts[chrom] = [iv.start for iv in ivs]

        self._total = sum(len(v) for v in self._by_chrom.values())
        logger.debug(f"Built interval index with {self._total} intervals across {len(self._by_chrom)} chromosomes")

    def query_overlap(self, chrom: str, start: int, end: int) -> list[GenomicInterval]:
        """Find all intervals overlapping the query region.

        An interval overlaps if interval.start < end AND interval.end > start.

        Args:
            chrom: Query chromosome.
            start: Query start position (inclusive).
            end: Query end position (exclusive).

        Returns:
            List of overlapping GenomicInterval objects.
        """
        if chrom not in self._by_chrom:
            return []

        intervals = self._by_chrom[chrom]
        starts = self._starts[chrom]

        # Find the leftmost interval that could overlap
        # An interval overlaps if interval.start < end AND interval.end > start
        # Use binary search to find intervals with start < end
        right_bound = bisect.bisect_left(starts, end)

        results: list[GenomicInterval] = []
        for i in range(right_bound):
            iv = intervals[i]
            if iv.end > start:
                results.append(iv)

        return results

    def query_nearest(
        self, chrom: str, position: int, max_distance: int = 1_000_000
    ) -> tuple[GenomicInterval | None, int]:
        """Find the nearest interval to a position.

        Args:
            chrom: Query chromosome.
            position: Query position.
            max_distance: Maximum distance to search.

        Returns:
            Tuple of (nearest_interval, distance). Distance is 0 if the
            position falls within the interval. Returns (None, -1) if
            no interval is within max_distance.
        """
        if chrom not in self._by_chrom:
            return None, -1

        intervals = self._by_chrom[chrom]
        starts = self._starts[chrom]

        # Binary search for position
        idx = bisect.bisect_left(starts, position)

        best_iv: GenomicInterval | None = None
        best_dist = max_distance + 1

        # Check the interval at idx and idx-1
        candidates = []
        if idx < len(intervals):
            candidates.append(idx)
        if idx > 0:
            candidates.append(idx - 1)
        if idx + 1 < len(intervals):
            candidates.append(idx + 1)

        for i in candidates:
            iv = intervals[i]
            if iv.start <= position < iv.end:
                return iv, 0
            dist = min(abs(position - iv.start), abs(position - iv.end))
            if dist < best_dist:
                best_dist = dist
                best_iv = iv

        if best_dist <= max_distance:
            return best_iv, best_dist
        return None, -1


def annotate_gene_overlap(
    variants: list[dict[str, Any]],
    gene_db: list[dict[str, Any]] | list[GenomicInterval],
) -> list[dict[str, Any]]:
    """Annotate structural variants with overlapping genes.

    For each variant, finds all genes that overlap with the variant
    interval and adds gene annotation information.

    Args:
        variants: List of variant dictionaries, each containing:
            - 'chrom': Chromosome
            - 'start': Start position
            - 'end': End position
            - Additional fields are preserved in the output.
        gene_db: List of gene entries as dictionaries with 'chrom', 'start',
            'end', 'name' keys, or as GenomicInterval objects.

    Returns:
        List of variant dictionaries with added fields:
            - 'overlapping_genes': List of gene names that overlap
            - 'gene_overlaps': List of OverlapResult objects
            - 'n_genes_affected': Number of genes affected
    """
    # Build interval index from gene database
    gene_intervals = _normalize_intervals(gene_db, feature_type="gene")
    index = IntervalIndex(gene_intervals)

    annotated: list[dict[str, Any]] = []

    for variant in variants:
        chrom = variant.get("chrom", "")
        start = variant.get("start", 0)
        end = variant.get("end", 0)

        overlapping = index.query_overlap(chrom, start, end)

        gene_names: list[str] = []
        overlap_results: list[OverlapResult] = []

        for gene in overlapping:
            overlap_bp = _compute_overlap_bp(start, end, gene.start, gene.end)
            variant_size = max(1, end - start)
            gene_size = max(1, gene.end - gene.start)

            frac_variant = overlap_bp / variant_size
            frac_gene = overlap_bp / gene_size

            relationship = _classify_relationship(start, end, gene.start, gene.end)

            gene_names.append(gene.name)
            overlap_results.append(
                OverlapResult(
                    variant_chrom=chrom,
                    variant_start=start,
                    variant_end=end,
                    feature=gene,
                    overlap_bp=overlap_bp,
                    overlap_fraction_variant=frac_variant,
                    overlap_fraction_feature=frac_gene,
                    relationship=relationship,
                )
            )

        result = dict(variant)
        result["overlapping_genes"] = gene_names
        result["gene_overlaps"] = overlap_results
        result["n_genes_affected"] = len(gene_names)
        annotated.append(result)

    total_annotated = sum(1 for v in annotated if v["n_genes_affected"] > 0)
    logger.info(f"Annotated {len(annotated)} variants: " f"{total_annotated} overlap with genes")
    return annotated


def annotate_regulatory_overlap(
    variants: list[dict[str, Any]],
    regulatory_db: list[dict[str, Any]] | list[GenomicInterval],
) -> list[dict[str, Any]]:
    """Annotate structural variants with overlapping regulatory elements.

    For each variant, finds all regulatory elements (enhancers, promoters,
    silencers, insulators, etc.) that overlap with the variant interval.

    Args:
        variants: List of variant dictionaries (same format as annotate_gene_overlap).
        regulatory_db: List of regulatory element entries as dictionaries with
            'chrom', 'start', 'end', 'name', 'feature_type' keys, or as
            GenomicInterval objects.

    Returns:
        List of variant dictionaries with added fields:
            - 'overlapping_regulatory': List of regulatory element names
            - 'regulatory_overlaps': List of OverlapResult objects
            - 'regulatory_types': Set of affected regulatory element types
            - 'n_regulatory_affected': Count of affected elements
    """
    reg_intervals = _normalize_intervals(regulatory_db, feature_type="regulatory")
    index = IntervalIndex(reg_intervals)

    annotated: list[dict[str, Any]] = []

    for variant in variants:
        chrom = variant.get("chrom", "")
        start = variant.get("start", 0)
        end = variant.get("end", 0)

        overlapping = index.query_overlap(chrom, start, end)

        reg_names: list[str] = []
        overlap_results: list[OverlapResult] = []
        reg_types: set[str] = set()

        for reg in overlapping:
            overlap_bp = _compute_overlap_bp(start, end, reg.start, reg.end)
            variant_size = max(1, end - start)
            reg_size = max(1, reg.end - reg.start)

            frac_variant = overlap_bp / variant_size
            frac_reg = overlap_bp / reg_size

            relationship = _classify_relationship(start, end, reg.start, reg.end)

            reg_names.append(reg.name)
            reg_types.add(reg.feature_type)
            overlap_results.append(
                OverlapResult(
                    variant_chrom=chrom,
                    variant_start=start,
                    variant_end=end,
                    feature=reg,
                    overlap_bp=overlap_bp,
                    overlap_fraction_variant=frac_variant,
                    overlap_fraction_feature=frac_reg,
                    relationship=relationship,
                )
            )

        result = dict(variant)
        result["overlapping_regulatory"] = reg_names
        result["regulatory_overlaps"] = overlap_results
        result["regulatory_types"] = reg_types
        result["n_regulatory_affected"] = len(reg_names)
        annotated.append(result)

    return annotated


def calculate_overlap_fraction(
    interval_a: dict[str, Any] | tuple[int, int],
    interval_b: dict[str, Any] | tuple[int, int],
) -> tuple[float, float]:
    """Calculate reciprocal overlap fractions between two intervals.

    Computes what fraction of interval A is covered by interval B and
    vice versa. Both intervals must be on the same chromosome (caller
    is responsible for this check).

    Args:
        interval_a: First interval as dict with 'start'/'end' keys or
            (start, end) tuple.
        interval_b: Second interval in the same format.

    Returns:
        Tuple of (fraction_of_a_covered, fraction_of_b_covered).
        Each value is between 0.0 and 1.0.
    """
    if isinstance(interval_a, tuple):
        a_start, a_end = interval_a
    else:
        a_start = interval_a.get("start", 0)
        a_end = interval_a.get("end", 0)

    if isinstance(interval_b, tuple):
        b_start, b_end = interval_b
    else:
        b_start = interval_b.get("start", 0)
        b_end = interval_b.get("end", 0)

    overlap_start = max(a_start, b_start)
    overlap_end = min(a_end, b_end)
    overlap_bp = max(0, overlap_end - overlap_start)

    a_size = max(1, a_end - a_start)
    b_size = max(1, b_end - b_start)

    return overlap_bp / a_size, overlap_bp / b_size


def find_nearest_gene(
    variant: dict[str, Any],
    gene_db: list[dict[str, Any]] | list[GenomicInterval],
    max_distance: int = 100_000,
) -> dict[str, Any]:
    """Find the nearest gene to a structural variant.

    If the variant overlaps a gene, the distance is 0. Otherwise, finds
    the closest gene within max_distance base pairs of either breakpoint.

    Args:
        variant: Variant dictionary with 'chrom', 'start', 'end' keys.
        gene_db: Gene database (see annotate_gene_overlap for format).
        max_distance: Maximum distance to search for nearest gene.

    Returns:
        Dictionary with:
            - 'nearest_gene': Name of nearest gene (empty string if none found)
            - 'distance': Distance to nearest gene in bp (0 if overlapping, -1 if none found)
            - 'direction': 'upstream', 'downstream', 'overlapping', or 'none'
            - 'gene_interval': The GenomicInterval of the nearest gene, or None
    """
    gene_intervals = _normalize_intervals(gene_db, feature_type="gene")
    index = IntervalIndex(gene_intervals)

    chrom = variant.get("chrom", "")
    start = variant.get("start", 0)
    end = variant.get("end", 0)

    # First check for overlap
    overlapping = index.query_overlap(chrom, start, end)
    if overlapping:
        # Return the gene with largest overlap
        best_gene = max(
            overlapping,
            key=lambda g: _compute_overlap_bp(start, end, g.start, g.end),
        )
        return {
            "nearest_gene": best_gene.name,
            "distance": 0,
            "direction": "overlapping",
            "gene_interval": best_gene,
        }

    # Search from both breakpoints
    midpoint = (start + end) // 2
    nearest, dist = index.query_nearest(chrom, midpoint, max_distance)

    if nearest is None:
        # Try from start and end
        nearest_start, dist_start = index.query_nearest(chrom, start, max_distance)
        nearest_end, dist_end = index.query_nearest(chrom, end, max_distance)

        if nearest_start is not None and (nearest_end is None or dist_start <= dist_end):
            nearest = nearest_start
            dist = dist_start
        elif nearest_end is not None:
            nearest = nearest_end
            dist = dist_end

    if nearest is None:
        return {
            "nearest_gene": "",
            "distance": -1,
            "direction": "none",
            "gene_interval": None,
        }

    # Determine direction
    gene_mid = (nearest.start + nearest.end) // 2
    variant_mid = (start + end) // 2
    direction = "upstream" if variant_mid < gene_mid else "downstream"

    return {
        "nearest_gene": nearest.name,
        "distance": dist,
        "direction": direction,
        "gene_interval": nearest,
    }


def _normalize_intervals(
    features: list[dict[str, Any]] | list[GenomicInterval],
    feature_type: str = "",
) -> list[GenomicInterval]:
    """Normalize a list of features to GenomicInterval objects.

    Args:
        features: List of dictionaries or GenomicInterval objects.
        feature_type: Default feature type if not specified in the data.

    Returns:
        List of GenomicInterval objects.
    """
    intervals: list[GenomicInterval] = []

    for feat in features:
        if isinstance(feat, GenomicInterval):
            intervals.append(feat)
        else:
            intervals.append(
                GenomicInterval(
                    chrom=feat.get("chrom", ""),
                    start=feat.get("start", 0),
                    end=feat.get("end", 0),
                    name=feat.get("name", ""),
                    feature_type=feat.get("feature_type", feature_type),
                    strand=feat.get("strand", "."),
                    info=feat.get("info", {}),
                )
            )

    return intervals


def _compute_overlap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Compute the number of overlapping base pairs between two intervals.

    Args:
        a_start: Start of interval A.
        a_end: End of interval A.
        b_start: Start of interval B.
        b_end: End of interval B.

    Returns:
        Number of overlapping base pairs (0 if no overlap).
    """
    overlap_start = max(a_start, b_start)
    overlap_end = min(a_end, b_end)
    return max(0, overlap_end - overlap_start)


def _classify_relationship(v_start: int, v_end: int, f_start: int, f_end: int) -> str:
    """Classify the spatial relationship between a variant and a feature.

    Args:
        v_start: Variant start.
        v_end: Variant end.
        f_start: Feature start.
        f_end: Feature end.

    Returns:
        Relationship string: 'contains', 'contained_in', or 'overlap'.
    """
    if v_start <= f_start and v_end >= f_end:
        return "contains"
    elif f_start <= v_start and f_end >= v_end:
        return "contained_in"
    else:
        return "overlap"
