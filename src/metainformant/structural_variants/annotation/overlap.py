"""Gene and regulatory element overlap annotation for structural variants.

Provides interval-based overlap detection between structural variants and
genomic features (genes, regulatory elements), nearest-gene finding, and
reciprocal overlap computation.

A translocation (``TRA``/``BND``) joins two loci, so it touches each feature at
two breakpoints, ``chrom``/``start`` and ``chrom2``/``pos2``. The functions here
look at both, each against the features of its own chromosome.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from typing import Any

from metainformant.core.utils.logging import get_logger
from metainformant.structural_variants._breakends import variant_regions

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
        relationship: Relationship type ('overlap', 'contained_in', 'contains').
    """

    variant_chrom: str
    variant_start: int
    variant_end: int
    feature: GenomicInterval
    overlap_bp: int = 0
    overlap_fraction_variant: float = 0.0
    overlap_fraction_feature: float = 0.0
    relationship: str = "overlap"


# Sentinel below every coordinate: the maximum of an empty range of interval ends.
_NO_END = -(1 << 62)

# A distance no pair of coordinates reaches: ``query_nearest`` with no search limit.
_NO_LIMIT = 1 << 62


class IntervalIndex:
    """Interval index for overlap and nearest-interval queries.

    Intervals are grouped by chromosome and sorted by ``(start, end)``. Two
    structures sit on that order, so a query never scans the whole chromosome:

    * a segment tree of the maximum ``end`` over ranges of the sorted intervals
      reports the ``k`` intervals that overlap a region in ``O(k log n)``,
      however long the intervals are (a gene spanning most of a chromosome is
      found, and does not slow down the queries after it);
    * the running maximum of the ends, together with the interval that reaches
      it, gives the nearest interval to the left of a position in ``O(1)``.
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

        self._starts: dict[str, list[int]] = {}
        self._max_end_tree: dict[str, list[int]] = {}
        self._tree_size: dict[str, int] = {}
        self._running_max_end: dict[str, list[int]] = {}
        self._running_max_arg: dict[str, list[int]] = {}

        for chrom, ivs in self._by_chrom.items():
            # Sort each chromosome's intervals by start position
            ivs.sort(key=lambda x: (x.start, x.end))
            self._starts[chrom] = [iv.start for iv in ivs]
            ends = [iv.end for iv in ivs]

            # Segment tree over the ends: node i has children 2i and 2i+1, the leaves start at ``size``.
            size = 1 << (len(ends) - 1).bit_length()
            tree = [_NO_END] * (2 * size)
            tree[size : size + len(ends)] = ends
            for node in range(size - 1, 0, -1):
                left, right = tree[2 * node], tree[2 * node + 1]
                tree[node] = left if left >= right else right
            self._max_end_tree[chrom] = tree
            self._tree_size[chrom] = size

            # Running maximum of the ends and the (last) interval that reaches it.
            running_end = _NO_END
            running_arg = -1
            max_ends: list[int] = []
            max_args: list[int] = []
            for i, end in enumerate(ends):
                if end >= running_end:
                    running_end = end
                    running_arg = i
                max_ends.append(running_end)
                max_args.append(running_arg)
            self._running_max_end[chrom] = max_ends
            self._running_max_arg[chrom] = max_args

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
            List of overlapping GenomicInterval objects, sorted by (start, end).
        """
        if chrom not in self._by_chrom:
            return []

        intervals = self._by_chrom[chrom]
        tree = self._max_end_tree[chrom]

        # Only intervals that start before the query end can overlap: a prefix of the sorted list.
        limit = bisect.bisect_left(self._starts[chrom], end)

        # Descend into every tree node that covers part of that prefix and holds an end beyond
        # the query start. Left children are visited first, so hits come out in sorted order.
        results: list[GenomicInterval] = []
        stack = [(1, 0, self._tree_size[chrom])]
        while stack:
            node, low, high = stack.pop()
            if low >= limit or tree[node] <= start:
                continue
            if high - low == 1:
                results.append(intervals[low])
                continue
            middle = (low + high) >> 1
            stack.append((2 * node + 1, middle, high))
            stack.append((2 * node, low, middle))

        return results

    def query_nearest(
        self, chrom: str, position: int, max_distance: int = 1_000_000
    ) -> tuple[GenomicInterval | None, int]:
        """Find the nearest interval to a position.

        Every interval on the chromosome is considered, however far its start lies
        from the position: an interval that spans the position is always found, and a
        long interval ending just before it is not hidden by the shorter ones that start
        after it.

        Args:
            chrom: Query chromosome.
            position: Query position.
            max_distance: Maximum distance to search.

        Returns:
            Tuple of (nearest_interval, distance). Distance is 0 if the
            position falls within the interval. Returns (None, -1) if
            no interval is within max_distance. On a tie the interval that starts
            at or after the position wins.
        """
        if chrom not in self._by_chrom:
            return None, -1

        intervals = self._by_chrom[chrom]

        # The position falls within an interval: the one starting last wins.
        containing = self.query_overlap(chrom, position, position + 1)
        if containing:
            return containing[-1], 0

        # No interval spans the position. Intervals before it all end at or before it, so the
        # closest of them is the one with the largest end; after it, the one with the smallest start.
        split = bisect.bisect_left(self._starts[chrom], position)

        best_iv: GenomicInterval | None = None
        best_dist = max_distance + 1

        if split < len(intervals):
            best_iv = intervals[split]
            best_dist = best_iv.start - position
        if split > 0:
            left = intervals[self._running_max_arg[chrom][split - 1]]
            left_dist = position - left.end
            if left_dist < best_dist:
                best_iv = left
                best_dist = left_dist

        if best_iv is not None and best_dist <= max_distance:
            return best_iv, best_dist
        return None, -1

    def query_nearest_region(
        self, chrom: str, start: int, end: int, max_distance: int = 1_000_000
    ) -> tuple[GenomicInterval | None, int]:
        """Find the interval nearest to the region ``[start, end)``.

        The distance is measured between the region and the interval (as ``bedtools closest`` does), not from a
        point of the region: the number of bases between them in half-open coordinates, ``interval.start - end``
        or ``start - interval.end``. It is 0 when they overlap, and also when they abut. A region that overlaps
        several intervals gets the one with the largest overlap (the first of them on a tie).

        Args:
            chrom: Query chromosome.
            start: Region start (inclusive).
            end: Region end (exclusive); equal to ``start`` for a point such as an insertion or a breakpoint.
            max_distance: Largest distance to report.

        Returns:
            Tuple of (nearest_interval, distance). Returns (None, -1) if no interval on the chromosome is within
            max_distance.
        """
        if chrom not in self._by_chrom:
            return None, -1

        overlapping = self.query_overlap(chrom, start, end)
        if overlapping:
            best_overlap = max(overlapping, key=lambda iv: _compute_overlap_bp(start, end, iv.start, iv.end))
            return best_overlap, 0

        # No interval overlaps the region, so each lies wholly to its left (end <= start) or to its right
        # (start >= end). The nearest on the left is the nearest interval to ``start`` among those ending before
        # it, and the nearest on the right is the nearest to ``end``; ``query_nearest`` from each of the two
        # positions finds one of them or the other, and together they always include the nearest overall.
        best_iv: GenomicInterval | None = None
        best_gap = max_distance + 1
        for position in (start, end):
            candidate, _ = self.query_nearest(chrom, position, _NO_LIMIT)
            if candidate is None:
                continue
            gap = max(candidate.start - end, start - candidate.end, 0)
            if gap < best_gap:
                best_iv, best_gap = candidate, gap

        if best_iv is not None and best_gap <= max_distance:
            return best_iv, best_gap
        return None, -1

    @classmethod
    def from_features(
        cls, features: list[dict[str, Any]] | list[GenomicInterval], feature_type: str = ""
    ) -> IntervalIndex:
        """Build an index from feature dictionaries (``chrom``, ``start``, ``end``, ``name``) or intervals.

        Build the index once and pass it in place of the feature list to :func:`annotate_gene_overlap`,
        :func:`annotate_regulatory_overlap` and :func:`find_nearest_gene` to annotate many variants without
        indexing the features again for each call.

        Args:
            features: Feature dictionaries (as accepted by :func:`annotate_gene_overlap`) or GenomicInterval objects.
            feature_type: Default feature type for dictionaries that give none.
        """
        return cls(_normalize_intervals(features, feature_type=feature_type))


def annotate_gene_overlap(
    variants: list[dict[str, Any]],
    gene_db: list[dict[str, Any]] | list[GenomicInterval] | IntervalIndex,
) -> list[dict[str, Any]]:
    """Annotate structural variants with overlapping genes.

    For each variant, finds all genes that overlap with the variant
    interval and adds gene annotation information. A translocation
    (``TRA``/``BND``) has no interval: the genes it disrupts are those that
    contain either breakpoint, ``chrom``/``start`` or ``chrom2``/``pos2``
    (the genes at the second breakpoint follow those at the first), and each
    ``OverlapResult`` records the breakpoint it was found at as its variant
    chromosome and position.

    Args:
        variants: List of variant dictionaries, each containing:
            - 'chrom': Chromosome
            - 'start': Start position
            - 'end': End position
            - 'chrom2', 'pos2': Second breakpoint (translocations)
            - Additional fields are preserved in the output.
        gene_db: List of gene entries as dictionaries with 'chrom', 'start',
            'end', 'name' keys, or as GenomicInterval objects; or an
            :class:`IntervalIndex` built once for many calls
            (``IntervalIndex.from_features(genes, "gene")``).

    Returns:
        List of variant dictionaries with added fields:
            - 'overlapping_genes': List of gene names that overlap
            - 'gene_overlaps': List of OverlapResult objects
            - 'n_genes_affected': Number of genes affected
    """
    index = _as_index(gene_db, feature_type="gene")

    annotated: list[dict[str, Any]] = []

    for variant in variants:
        gene_names: list[str] = []
        overlap_results: list[OverlapResult] = []
        seen: set[int] = set()  # genes found at an earlier breakpoint of a translocation

        for chrom, start, end in variant_regions(variant):
            found = index.query_overlap(chrom, start, end)
            for gene in found:
                if id(gene) in seen:
                    continue  # both breakpoints of a translocation fall in this gene
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
            seen.update(id(gene) for gene in found)

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
    regulatory_db: list[dict[str, Any]] | list[GenomicInterval] | IntervalIndex,
) -> list[dict[str, Any]]:
    """Annotate structural variants with overlapping regulatory elements.

    For each variant, finds all regulatory elements (enhancers, promoters,
    silencers, insulators, etc.) that overlap with the variant interval. A
    translocation (``TRA``/``BND``) is checked at each of its two breakpoints
    (``chrom``/``start`` and ``chrom2``/``pos2``), as in :func:`annotate_gene_overlap`.

    Args:
        variants: List of variant dictionaries (same format as annotate_gene_overlap).
        regulatory_db: List of regulatory element entries as dictionaries with
            'chrom', 'start', 'end', 'name', 'feature_type' keys, or as
            GenomicInterval objects; or an :class:`IntervalIndex` built once
            for many calls (``IntervalIndex.from_features(elements, "regulatory")``).

    Returns:
        List of variant dictionaries with added fields:
            - 'overlapping_regulatory': List of regulatory element names
            - 'regulatory_overlaps': List of OverlapResult objects
            - 'regulatory_types': Set of affected regulatory element types
            - 'n_regulatory_affected': Count of affected elements
    """
    index = _as_index(regulatory_db, feature_type="regulatory")

    annotated: list[dict[str, Any]] = []

    for variant in variants:
        reg_names: list[str] = []
        overlap_results: list[OverlapResult] = []
        reg_types: set[str] = set()
        seen: set[int] = set()  # elements found at an earlier breakpoint of a translocation

        for chrom, start, end in variant_regions(variant):
            found = index.query_overlap(chrom, start, end)
            for reg in found:
                if id(reg) in seen:
                    continue  # both breakpoints of a translocation fall in this element
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
            seen.update(id(reg) for reg in found)

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
    gene_db: list[dict[str, Any]] | list[GenomicInterval] | IntervalIndex,
    max_distance: int = 100_000,
) -> dict[str, Any]:
    """Find the nearest gene to a structural variant.

    The distance runs from the variant *interval* to the gene interval, as
    ``bedtools closest`` measures it, and not from the variant's midpoint: the
    number of bases between them in half-open coordinates (``gene.start -
    variant.end`` or ``variant.start - gene.end``). It is 0 when the variant
    overlaps a gene, in which case the gene with the largest overlap is
    returned, and also when a gene merely abuts the variant (``'overlapping'``
    is reported only for a real overlap). A point variant (``start == end``,
    such as an insertion) overlaps a gene that contains it strictly.

    A translocation (``TRA``/``BND``) has no interval: its two breakpoints,
    ``chrom``/``start`` and ``chrom2``/``pos2``, are each looked up on their own
    chromosome and the nearer gene is returned (the first breakpoint's on a tie).

    Indexing the genes costs more than looking one variant up. To annotate many
    variants, build the index once (``IntervalIndex.from_features(genes, "gene")``)
    and pass it as ``gene_db``; a list of genes is indexed again on every call.

    Args:
        variant: Variant dictionary with 'chrom', 'start', 'end' keys.
        gene_db: Gene database (see annotate_gene_overlap for format), or an
            :class:`IntervalIndex` of the genes.
        max_distance: Maximum distance to search for nearest gene.

    Returns:
        Dictionary with:
            - 'nearest_gene': Name of nearest gene (empty string if none found)
            - 'distance': Distance to nearest gene in bp (0 if overlapping, -1 if none found)
            - 'direction': 'upstream' or 'downstream' (the variant lies before or after
              the gene in genomic coordinates, whatever the gene's strand),
              'overlapping', or 'none'
            - 'gene_interval': The GenomicInterval of the nearest gene, or None
    """
    index = _as_index(gene_db, feature_type="gene")

    nearest: GenomicInterval | None = None
    distance = -1
    region = (0, 0)
    for chrom, start, end in variant_regions(variant):
        gene, gene_distance = index.query_nearest_region(chrom, start, end, max_distance)
        if gene is not None and (nearest is None or gene_distance < distance):
            nearest, distance, region = gene, gene_distance, (start, end)

    if nearest is None:
        return {
            "nearest_gene": "",
            "distance": -1,
            "direction": "none",
            "gene_interval": None,
        }

    start, end = region
    if nearest.start < end and nearest.end > start:
        direction = "overlapping"
    else:
        direction = "upstream" if end <= nearest.start else "downstream"

    return {
        "nearest_gene": nearest.name,
        "distance": distance,
        "direction": direction,
        "gene_interval": nearest,
    }


def _as_index(
    features: list[dict[str, Any]] | list[GenomicInterval] | IntervalIndex, feature_type: str = ""
) -> IntervalIndex:
    """Return ``features`` as an interval index: a prebuilt index is used as is, a list is indexed."""
    if isinstance(features, IntervalIndex):
        return features
    return IntervalIndex(_normalize_intervals(features, feature_type=feature_type))


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
