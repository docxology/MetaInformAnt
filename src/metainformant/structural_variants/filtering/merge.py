"""Multi-caller merging and consensus for structural variants.

Implements algorithms for merging SV callsets from multiple detection
tools, computing reciprocal overlap between SVs, SURVIVOR-style
distance-based merging, and deduplication of variant calls.

Intra-chromosomal variants are matched on their interval. A translocation
(``TRA``/``BND``) is a zero-length record, so reciprocal overlap can never join
two of them: it is matched on both of its breakpoints instead (``chrom`` and
``start``, plus ``chrom2`` and ``pos2``), each pair of positions within the
merge distance, whichever breakpoint a caller wrote first. A translocation never
matches an interval variant. See :mod:`metainformant.structural_variants._breakends`
for how the second breakpoint is read.
"""

from __future__ import annotations

import math
import re
import statistics
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger
from metainformant.structural_variants._breakends import (
    TRANSLOCATION_TYPES,
    canonical_junction,
    consensus_junction,
    is_translocation,
    junction_links,
)

np: Any  # None when numpy is unavailable (optional dependency)
try:
    import numpy as np
except ImportError:
    np = None

logger = get_logger(__name__)


@dataclass
class MergedVariant:
    """A structural variant resulting from merging multiple callsets.

    Attributes:
        chrom: Chromosome.
        start: Consensus start position.
        end: Consensus end position.
        sv_type: Consensus SV type.
        size: Consensus size.
        n_callers: Number of callers that detected this variant.
        caller_names: Names of callers that detected this variant (sorted).
        support_variants: Original variant calls that were merged.
        consensus_quality: Quality based on multi-caller agreement (the rule
            differs per merge function; see :func:`merge_callsets` and
            :func:`survivor_merge`).
        genotype: Consensus genotype ('./.' when no merged call has one).
        chrom2: Chromosome of the second breakpoint of a merged translocation
            ('' for every other variant).
        pos2: Position of the second breakpoint on ``chrom2`` (None unless the
            merged variant is a translocation with a known partner). A merged
            translocation lists its breakpoints in (chromosome, position)
            order, with ``start`` and ``end`` both the first one and ``size`` 0.
    """

    chrom: str
    start: int
    end: int
    sv_type: str
    size: int = 0
    n_callers: int = 0
    caller_names: list[str] = field(default_factory=list)
    support_variants: list[dict[str, Any]] = field(default_factory=list)
    consensus_quality: float = 0.0
    genotype: str = "./."
    chrom2: str = ""
    pos2: int | None = None


@dataclass
class MergeStats:
    """Statistics from a merge operation.

    Attributes:
        n_input_callsets: Number of input callsets.
        n_input_variants: Total variants across all callsets.
        n_output_variants: Number of merged variants.
        n_multi_caller: Variants detected by multiple callers.
        n_single_caller: Variants detected by only one caller.
        caller_counts: Variants per caller.
    """

    n_input_callsets: int
    n_input_variants: int
    n_output_variants: int
    n_multi_caller: int = 0
    n_single_caller: int = 0
    caller_counts: dict[str, int] = field(default_factory=dict)


def merge_callsets(
    callsets: dict[str, list[dict[str, Any]]],
    min_overlap: float = 0.5,
    type_match: bool = True,
    max_distance: int = 1000,
) -> tuple[list[MergedVariant], MergeStats]:
    """Merge structural variant callsets from multiple callers.

    Uses reciprocal overlap to identify the same SV detected by different
    tools. For each group of matching variants across callers, produces a
    single consensus call with aggregated support: median breakpoints, the
    majority SV type, the majority genotype (see :func:`_consensus_genotype`)
    and ``consensus_quality`` = the mean ``quality`` (or ``qual``; 0 when
    absent or None) of the merged calls plus ``100 * log10(n_callers + 1)``.

    A translocation (``TRA``/``BND``) has no interval, so reciprocal overlap
    cannot match it. It is matched instead on both breakpoints: ``chrom`` and
    ``start`` of one call against the other's, and ``chrom2`` and ``pos2``
    likewise, each pair of positions within ``max_distance``. A junction
    recorded from its other side (the second breakpoint first) matches too. A
    translocation never merges with an interval variant, and ``min_overlap``
    does not apply to it. The consensus of a group of translocations takes the
    median of each breakpoint and lists them in (chromosome, position) order,
    with the second one in ``chrom2``/``pos2`` of the :class:`MergedVariant`.

    Args:
        callsets: Dictionary mapping caller names to lists of variant
            dictionaries. Each variant should have 'chrom', 'start',
            'end', 'sv_type' keys; a translocation also 'chrom2' and 'pos2'.
        min_overlap: Minimum reciprocal overlap fraction (0-1) to consider
            two variants as the same event (default 0.5 = 50%).
        type_match: Whether SV types must match for merging (default True).
            ``TRA`` and ``BND`` count as the same type.
        max_distance: Largest distance in bp between the corresponding
            breakpoints of two translocations that are the same event
            (default 1000, the SURVIVOR default).

    Returns:
        Tuple of (merged_variants, merge_statistics).
    """
    if not callsets:
        return [], MergeStats(0, 0, 0)

    # Tag each variant with its caller
    all_variants: list[tuple[str, dict[str, Any]]] = []
    caller_counts: dict[str, int] = {}

    for caller_name, variants in callsets.items():
        caller_counts[caller_name] = len(variants)
        for v in variants:
            tagged = dict(v)
            tagged["_caller"] = caller_name
            all_variants.append((caller_name, tagged))

    total_input = len(all_variants)

    # Sort by chromosome and start position
    all_variants.sort(key=lambda x: (x[1].get("chrom", ""), x[1].get("start", 0)))

    # Graph-based clustering: build adjacency lists for matching variants
    n = len(all_variants)
    merged_into: list[int] = list(range(n))  # Union-find parent array

    for i in range(n):
        v_i = all_variants[i][1]
        if is_translocation(v_i):
            continue  # matched on both breakpoints below
        chrom_i = v_i.get("chrom", "")
        end_i = v_i.get("end", 0)
        type_i = v_i.get("sv_type", "")
        if hasattr(type_i, "value"):
            type_i = type_i.value

        for j in range(i + 1, n):
            v_j = all_variants[j][1]
            chrom_j = v_j.get("chrom", "")

            if chrom_j != chrom_i:
                break  # Sorted by chrom, no more matches

            start_j = v_j.get("start", 0)

            # Quick distance check before computing overlap
            if start_j > end_i + 10_000:
                break

            if is_translocation(v_j):
                continue  # a translocation never matches an interval

            type_j = v_j.get("sv_type", "")
            if hasattr(type_j, "value"):
                type_j = type_j.value

            # Type matching
            if type_match and type_i and type_j and type_i != type_j:
                continue

            # Don't merge variants from the same caller
            if all_variants[i][0] == all_variants[j][0]:
                continue

            # Reciprocal overlap check
            ro = calculate_reciprocal_overlap(v_i, v_j)
            if ro >= min_overlap:
                # Union-find merge
                _union(merged_into, i, j)

    # Translocations: the same junction when both breakpoints agree (TRA and BND are one type)
    for i, j in _translocation_links([v for _, v in all_variants], max_distance):
        if all_variants[i][0] != all_variants[j][0]:  # Don't merge variants from the same caller
            _union(merged_into, i, j)

    # Group variants by their union-find root
    groups: dict[int, list[int]] = {}
    for i in range(n):
        root = _find(merged_into, i)
        if root not in groups:
            groups[root] = []
        groups[root].append(i)

    # Build merged variants from groups
    merged_list: list[MergedVariant] = []

    for indices in groups.values():
        group_variants = [all_variants[i] for i in indices]
        callers = sorted({caller for caller, _ in group_variants})
        variant_dicts = [v for _, v in group_variants]

        chrom, consensus_start, consensus_end, chrom2, pos2 = _group_position(variant_dicts)
        consensus_type = _majority_sv_type(variant_dicts)

        size = abs(consensus_end - consensus_start)

        # Quality from multi-caller support
        n_callers = len(callers)
        base_qual = sum(_call_quality(v) for v in variant_dicts) / len(variant_dicts)
        consensus_quality = base_qual + 10.0 * math.log10(n_callers + 1) * 10

        # Consensus genotype
        genotypes = [v.get("genotype", "./.") for v in variant_dicts]
        consensus_gt = _consensus_genotype(genotypes)

        merged = MergedVariant(
            chrom=chrom,
            start=consensus_start,
            end=consensus_end,
            sv_type=consensus_type,
            size=size,
            n_callers=n_callers,
            caller_names=callers,
            support_variants=variant_dicts,
            consensus_quality=consensus_quality,
            genotype=consensus_gt,
            chrom2=chrom2,
            pos2=pos2,
        )
        merged_list.append(merged)

    # Sort by position
    merged_list.sort(key=lambda m: (m.chrom, m.start))

    n_multi = sum(1 for m in merged_list if m.n_callers > 1)
    n_single = sum(1 for m in merged_list if m.n_callers == 1)

    stats = MergeStats(
        n_input_callsets=len(callsets),
        n_input_variants=total_input,
        n_output_variants=len(merged_list),
        n_multi_caller=n_multi,
        n_single_caller=n_single,
        caller_counts=caller_counts,
    )

    logger.info(
        f"Merged {stats.n_input_variants} variants from {stats.n_input_callsets} callers "
        f"into {stats.n_output_variants} ({n_multi} multi-caller, {n_single} single-caller)"
    )
    return merged_list, stats


def calculate_reciprocal_overlap(
    sv1: dict[str, Any],
    sv2: dict[str, Any],
) -> float:
    """Calculate the minimum reciprocal overlap between two SVs.

    Reciprocal overlap requires that the overlap constitutes at least a
    certain fraction of BOTH variants. The minimum of the two fractions
    is returned (i.e., the stricter criterion).

    Variants on different chromosomes (or without coordinates) never
    overlap and yield 0.0.

    Args:
        sv1: First variant dictionary with 'chrom', 'start', 'end' keys.
        sv2: Second variant dictionary with 'chrom', 'start', 'end' keys.

    Returns:
        Reciprocal overlap fraction between 0.0 and 1.0.
    """
    chrom1 = sv1.get("chrom", "")
    chrom2 = sv2.get("chrom", "")

    if chrom1 != chrom2:
        return 0.0

    start1: int = sv1.get("start", 0)
    end1: int = sv1.get("end", 0)
    start2: int = sv2.get("start", 0)
    end2: int = sv2.get("end", 0)

    size1 = max(1, end1 - start1)
    size2 = max(1, end2 - start2)

    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap_bp = max(0, overlap_end - overlap_start)

    if overlap_bp == 0:
        return 0.0

    frac1 = overlap_bp / size1
    frac2 = overlap_bp / size2

    return min(frac1, frac2)


def survivor_merge(
    vcf_files: list[str] | list[dict[str, Any]],
    max_distance: int = 1000,
    min_callers: int = 2,
    type_match: bool = True,
    strand_match: bool = False,
) -> tuple[list[MergedVariant], MergeStats]:
    """SURVIVOR-style merging of structural variant callsets.

    Implements the SURVIVOR merge algorithm that uses breakpoint distance
    rather than reciprocal overlap for matching variants. This is more
    appropriate for imprecise SV calls where exact breakpoints are uncertain.

    Note: If vcf_files contains strings (file paths), this function requires
    pysam or a VCF parser. If vcf_files contains dictionaries, they are
    used directly.

    Each merged variant combines the calls of its group as follows:

    * ``genotype``: majority vote among the calls' non-missing genotypes
      (allele order ignored, so '1/0' is '0/1'). A tie goes to the genotype
      whose calls have the larger summed quality, then to the one that sorts
      first. ``'./.'`` when no call has a genotype. Genotypes are read from
      each dict's ``genotype`` and from the first sample of VCF files.
    * ``consensus_quality``: the mean ``quality`` (or ``qual``; 0 when absent)
      of the merged calls, plus an agreement bonus of
      ``10 * n_callers * log10(n_callers + 1)``.

    A translocation (``TRA``/``BND``) is compared on both of its breakpoints,
    not on its interval: ``chrom`` and ``start`` of one call against the other's,
    and ``chrom2`` and ``pos2`` likewise, each pair of positions within
    ``max_distance``. A junction recorded from its other side (the second
    breakpoint first) matches too, and a translocation never merges with an
    interval variant. The consensus of a group of translocations takes the
    median of each breakpoint and lists them in (chromosome, position) order,
    the second one in ``chrom2``/``pos2`` of the :class:`MergedVariant`. VCF
    records give their second breakpoint through a BND ALT (``N[chr2:500[``) or
    ``CHR2`` with ``END``/``POS2`` in INFO.

    Args:
        vcf_files: Either file paths to VCF files or lists of variant
            dictionaries. If file paths, each file is treated as one caller.
            If dicts, each dict must have a '_caller' key naming its caller:
            a call without one cannot be assigned to a caller, and counting
            it as a caller of its own would inflate ``n_callers``.
        max_distance: Maximum distance between breakpoints to consider
            variants as the same event (default 1000bp).
        min_callers: Minimum number of callers required to keep a merged
            variant (default 2).
        type_match: Whether SV types must match (default True). ``TRA`` and
            ``BND`` count as the same type.
        strand_match: Whether strand orientations must match (default False).

    Returns:
        Tuple of (merged_variants, merge_statistics).

    Raises:
        ValueError: If a dict in ``vcf_files`` has no '_caller' key.
        TypeError: If an item of ``vcf_files`` is not a path, a dict or a list.
    """
    # Load variants from VCF files or use provided dicts
    callsets: dict[str, list[dict[str, Any]]] = {}

    for i, item in enumerate(vcf_files):
        if isinstance(item, str):
            # Load from VCF file
            caller_name = f"caller_{i}"
            variants = _load_vcf_variants(item)
            callsets[caller_name] = variants
        elif isinstance(item, dict):
            if item.get("_caller") is None:
                raise ValueError(
                    f"vcf_files[{i}] is a dict without a '_caller' key: tag each call with the name of the "
                    "caller that made it (a call counted as a caller of its own would inflate n_callers)"
                )
            caller_name = item["_caller"]
            if caller_name not in callsets:
                callsets[caller_name] = []
            callsets[caller_name].append(item)
        elif isinstance(item, list):
            caller_name = f"caller_{i}"
            callsets[caller_name] = item
        else:
            raise TypeError(
                f"vcf_files[{i}] is a {type(item).__name__}: expected a VCF path (str), a variant dict with a "
                "'_caller' key, or a list of variant dicts"
            )

    # Flatten all variants with caller tags
    all_variants: list[tuple[str, dict[str, Any]]] = []
    caller_counts: dict[str, int] = {}
    for caller_name, variants in callsets.items():
        caller_counts[caller_name] = len(variants)
        for v in variants:
            all_variants.append((caller_name, v))

    total_input = len(all_variants)

    # Sort by chromosome and position
    all_variants.sort(key=lambda x: (x[1].get("chrom", ""), x[1].get("start", 0)))

    # Distance-based merging
    n = len(all_variants)
    parent: list[int] = list(range(n))

    for i in range(n):
        v_i = all_variants[i][1]
        if is_translocation(v_i):
            continue  # matched on both breakpoints below
        chrom_i = v_i.get("chrom", "")
        start_i = v_i.get("start", 0)
        end_i = v_i.get("end", 0)
        type_i = v_i.get("sv_type", "")
        if hasattr(type_i, "value"):
            type_i = type_i.value
        strand_i = v_i.get("strand", "+")

        for j in range(i + 1, n):
            v_j = all_variants[j][1]
            chrom_j = v_j.get("chrom", "")

            if chrom_j != chrom_i:
                break

            start_j = v_j.get("start", 0)

            # Quick exit: if starts are too far apart
            if start_j - start_i > max_distance + (end_i - start_i):
                break

            if is_translocation(v_j):
                continue  # a translocation never matches an interval

            end_j = v_j.get("end", 0)
            type_j = v_j.get("sv_type", "")
            if hasattr(type_j, "value"):
                type_j = type_j.value
            strand_j = v_j.get("strand", "+")

            # Don't merge same-caller variants
            if all_variants[i][0] == all_variants[j][0]:
                continue

            # Type match check
            if type_match and type_i and type_j and type_i != type_j:
                continue

            # Strand match check
            if strand_match and strand_i != strand_j:
                continue

            # Distance check: breakpoints must be within max_distance
            dist_start = abs(start_i - start_j)
            dist_end = abs(end_i - end_j)

            if dist_start <= max_distance and dist_end <= max_distance:
                _union(parent, i, j)

    # Translocations: the same junction when both breakpoints are within max_distance
    for i, j in _translocation_links([v for _, v in all_variants], max_distance):
        if all_variants[i][0] == all_variants[j][0]:
            continue  # Don't merge same-caller variants
        if strand_match and all_variants[i][1].get("strand", "+") != all_variants[j][1].get("strand", "+"):
            continue
        _union(parent, i, j)

    # Group by root
    groups: dict[int, list[int]] = {}
    for i in range(n):
        root = _find(parent, i)
        if root not in groups:
            groups[root] = []
        groups[root].append(i)

    # Build merged variants, filtering by min_callers
    merged_list: list[MergedVariant] = []

    for indices in groups.values():
        group_variants = [all_variants[i] for i in indices]
        callers = sorted({caller for caller, _ in group_variants})

        if len(callers) < min_callers:
            continue

        variant_dicts = [v for _, v in group_variants]

        chrom, consensus_start, consensus_end, chrom2, pos2 = _group_position(variant_dicts)
        consensus_type = _majority_sv_type(variant_dicts)

        size = abs(consensus_end - consensus_start)
        n_callers = len(callers)

        # Callers' evidence: the mean call quality, plus a bonus for how many callers agree.
        call_qualities = [_call_quality(v) for v in variant_dicts]
        quality = sum(call_qualities) / len(call_qualities) + 10.0 * n_callers * math.log10(n_callers + 1)
        genotype = _consensus_genotype([v.get("genotype", "./.") for v in variant_dicts], call_qualities)

        merged = MergedVariant(
            chrom=chrom,
            start=consensus_start,
            end=consensus_end,
            sv_type=consensus_type,
            size=size,
            n_callers=n_callers,
            caller_names=callers,
            support_variants=variant_dicts,
            consensus_quality=quality,
            genotype=genotype,
            chrom2=chrom2,
            pos2=pos2,
        )
        merged_list.append(merged)

    merged_list.sort(key=lambda m: (m.chrom, m.start))

    n_multi = sum(1 for m in merged_list if m.n_callers > 1)

    stats = MergeStats(
        n_input_callsets=len(callsets),
        n_input_variants=total_input,
        n_output_variants=len(merged_list),
        n_multi_caller=n_multi,
        n_single_caller=len(merged_list) - n_multi,
        caller_counts=caller_counts,
    )

    logger.info(
        f"SURVIVOR merge: {stats.n_output_variants} variants from {stats.n_input_variants} "
        f"(min_callers={min_callers}, max_distance={max_distance}bp)"
    )
    return merged_list, stats


def deduplicate_variants(
    variants: list[dict[str, Any]],
    max_distance: int = 100,
    min_reciprocal_overlap: float = 0.8,
    type_match: bool = True,
) -> tuple[list[dict[str, Any]], int]:
    """Remove duplicate structural variant calls.

    Identifies duplicate calls based on proximity and reciprocal overlap,
    keeping the highest-quality call from each group of duplicates.

    A translocation (``TRA``/``BND``) has no interval, so reciprocal overlap
    cannot identify its duplicates. Two translocations are duplicates when both
    breakpoints agree: ``chrom`` and ``start`` of one against the other's, and
    ``chrom2`` and ``pos2`` likewise, each pair of positions within
    ``max_distance`` (``min_reciprocal_overlap`` does not apply), whichever
    breakpoint each call wrote first. A translocation is never a duplicate of
    an interval variant.

    Args:
        variants: List of variant dictionaries with 'chrom', 'start', 'end',
            'sv_type', and optionally 'quality' keys (a translocation also
            'chrom2' and 'pos2').
        max_distance: Maximum breakpoint distance for duplicate candidates.
        min_reciprocal_overlap: Minimum reciprocal overlap to classify as
            duplicate (default 0.8 = 80%).
        type_match: Whether SV types must match (default True). ``TRA`` and
            ``BND`` count as the same type.

    Returns:
        Tuple of (deduplicated_variants, n_removed).
    """
    if not variants:
        return [], 0

    # Sort by chromosome and position
    sorted_vars = sorted(variants, key=lambda v: (v.get("chrom", ""), v.get("start", 0)))
    n = len(sorted_vars)

    # Mark duplicates using union-find
    parent: list[int] = list(range(n))

    for i in range(n):
        v_i = sorted_vars[i]
        if is_translocation(v_i):
            continue  # matched on both breakpoints below
        chrom_i = v_i.get("chrom", "")
        start_i = v_i.get("start", 0)
        end_i = v_i.get("end", 0)
        type_i = v_i.get("sv_type", "")
        if hasattr(type_i, "value"):
            type_i = type_i.value

        for j in range(i + 1, n):
            v_j = sorted_vars[j]
            chrom_j = v_j.get("chrom", "")

            if chrom_j != chrom_i:
                break

            start_j = v_j.get("start", 0)
            if start_j - start_i > max_distance + (end_i - start_i):
                break

            if is_translocation(v_j):
                continue  # a translocation is never a duplicate of an interval

            type_j = v_j.get("sv_type", "")
            if hasattr(type_j, "value"):
                type_j = type_j.value

            if type_match and type_i and type_j and type_i != type_j:
                continue

            # Check distance and overlap
            end_j = v_j.get("end", 0)
            dist_start = abs(start_i - start_j)
            dist_end = abs(end_i - end_j)

            if dist_start <= max_distance and dist_end <= max_distance:
                ro = calculate_reciprocal_overlap(v_i, v_j)
                if ro >= min_reciprocal_overlap:
                    _union(parent, i, j)

    # Translocations: duplicates when both breakpoints are within max_distance
    for i, j in _translocation_links(sorted_vars, max_distance):
        _union(parent, i, j)

    # Group duplicates and keep highest quality
    groups: dict[int, list[int]] = {}
    for i in range(n):
        root = _find(parent, i)
        if root not in groups:
            groups[root] = []
        groups[root].append(i)

    deduped: list[dict[str, Any]] = []
    for indices in groups.values():
        # Pick the variant with highest quality
        best_idx = max(indices, key=lambda i: _call_quality(sorted_vars[i]))
        deduped.append(sorted_vars[best_idx])

    # Sort output
    deduped.sort(key=lambda v: (v.get("chrom", ""), v.get("start", 0)))

    n_removed = n - len(deduped)
    logger.info(f"Deduplication: {len(deduped)} unique variants ({n_removed} duplicates removed)")

    return deduped, n_removed


def _union(parent: list[int], i: int, j: int) -> None:
    """Union operation for union-find with path compression.

    Args:
        parent: Parent array.
        i: First element.
        j: Second element.
    """
    root_i = _find(parent, i)
    root_j = _find(parent, j)
    if root_i != root_j:
        parent[root_j] = root_i


def _find(parent: list[int], i: int) -> int:
    """Find operation for union-find with path compression.

    Args:
        parent: Parent array.
        i: Element to find root of.

    Returns:
        Root element.
    """
    while parent[i] != i:
        parent[i] = parent[parent[i]]  # Path compression
        i = parent[i]
    return i


def _median_position(values: list[int]) -> int:
    """Median of a list of positions (consensus breakpoint).

    Args:
        values: Positions from the variants in a merge group.

    Returns:
        Median position as an int.
    """
    if np is not None:
        return int(np.median(values))
    return int(statistics.median(values))


def _translocation_links(variants: Sequence[Mapping[str, Any]], max_distance: int) -> Iterator[tuple[int, int]]:
    """Index pairs of translocation records that are the same junction.

    Args:
        variants: Variant records (interval variants are skipped).
        max_distance: Largest distance between corresponding breakpoints.

    Yields:
        ``(i, j)`` indices into ``variants`` of two translocations whose breakpoints agree on both chromosomes
        and lie within ``max_distance``, in either orientation.
    """
    indices = [i for i, variant in enumerate(variants) if is_translocation(variant)]
    if len(indices) < 2:
        return
    junctions = [canonical_junction(variants[i]) for i in indices]
    for a, b in junction_links(junctions, max_distance):
        yield indices[a], indices[b]


def _group_position(variant_dicts: Sequence[Mapping[str, Any]]) -> tuple[str, int, int, str, int | None]:
    """Consensus position of a merge group.

    Args:
        variant_dicts: The variants of one group. A group holds either interval variants or translocations,
            never both (they are never linked).

    Returns:
        ``(chrom, start, end, chrom2, pos2)``: the median start and end of interval variants (``chrom2`` ''
        and ``pos2`` None), or the median of each breakpoint of translocations, in (chromosome, position)
        order, with ``end == start``.
    """
    if is_translocation(variant_dicts[0]):
        (chrom, start), second = consensus_junction([canonical_junction(v) for v in variant_dicts], _median_position)
        chrom2, pos2 = ("", None) if second is None else second
        return chrom, start, start, chrom2, pos2
    return (
        variant_dicts[0].get("chrom", ""),
        _median_position([v.get("start", 0) for v in variant_dicts]),
        _median_position([v.get("end", 0) for v in variant_dicts]),
        "",
        None,
    )


def _majority_sv_type(variant_dicts: list[dict[str, Any]]) -> str:
    """Consensus SV type by majority vote within a merge group.

    Args:
        variant_dicts: Variants in the merge group.

    Returns:
        Most frequent SV type (enum values are unwrapped first).
    """
    type_counts: dict[str, int] = {}
    for v in variant_dicts:
        st = v.get("sv_type", "UNKNOWN")
        if hasattr(st, "value"):
            st = st.value
        type_counts[st] = type_counts.get(st, 0) + 1
    return max(type_counts, key=type_counts.get)  # type: ignore[arg-type]


def _call_quality(variant: dict[str, Any]) -> float:
    """Quality of one call (``quality``, else ``qual``); 0 when absent or None."""
    quality = variant.get("quality", variant.get("qual", 0.0))
    return 0.0 if quality is None else float(quality)


def _consensus_genotype(genotypes: list[str], qualities: Sequence[float] | None = None) -> str:
    """Determine consensus genotype from multiple calls.

    Uses majority voting among non-missing genotypes, ignoring allele order.

    Args:
        genotypes: List of genotype strings.
        qualities: Optional quality of each call, parallel to ``genotypes``. When
            given, a tie in the vote goes to the genotype with the larger summed
            quality, then to the one that sorts first. Without it a tie goes to the
            genotype seen first.

    Returns:
        Consensus genotype string.
    """
    weights = qualities if qualities is not None else [0.0] * len(genotypes)
    counts: dict[str, int] = {}
    total_quality: dict[str, float] = {}
    for g, quality in zip(genotypes, weights):
        if not isinstance(g, str) or g in ("./.", "."):
            continue
        # Normalize genotype ordering
        normalized = "/".join(sorted(g.split("/")))
        counts[normalized] = counts.get(normalized, 0) + 1
        total_quality[normalized] = total_quality.get(normalized, 0.0) + quality

    if not counts:
        return "./."
    if qualities is None:
        return max(counts, key=counts.get)  # type: ignore[arg-type]
    return min(counts, key=lambda g: (-counts[g], -total_quality[g], g))


def _format_genotype(alleles: Sequence[int | None] | None) -> str | None:
    """Format a called allele tuple (as pysam reports ``GT``) as ``'0/1'``; None when absent."""
    if not alleles:
        return None
    return "/".join("." if allele is None else str(allele) for allele in alleles)


def _vcf_sample_genotype(fields: list[str]) -> str | None:
    """Genotype of the first sample of a VCF record line (``'0/1'``), or None when there is none.

    Phased genotypes (``0|1``) are written with ``/``; a missing call is ``'./.'``.
    """
    if len(fields) < 10:
        return None
    keys = fields[8].split(":")
    if "GT" not in keys:
        return None
    values = fields[9].split(":")
    index = keys.index("GT")
    if index >= len(values):
        return None
    genotype = values[index].replace("|", "/")
    return "./." if genotype in ("", ".") else genotype


def _load_vcf_variants(vcf_path: str) -> list[dict[str, Any]]:
    """Load structural variants from a VCF file.

    Attempts to use pysam for VCF parsing. Falls back to basic text
    parsing if pysam is not available.

    Args:
        vcf_path: Path to VCF file.

    Returns:
        List of variant dictionaries.
    """
    try:
        import pysam

        variants: list[dict[str, Any]] = []
        with pysam.VariantFile(vcf_path) as vcf:
            for record in vcf:
                sv_type = record.info.get("SVTYPE", "UNKNOWN")
                end = record.info.get("END", record.stop)
                variant = {
                    "chrom": record.chrom,
                    "start": record.start,
                    "end": end,
                    "sv_type": sv_type,
                    "quality": record.qual or 0.0,
                    "filter": list(record.filter),
                }
                if sv_type in TRANSLOCATION_TYPES:
                    alt = record.alts[0] if record.alts else None
                    chrom2 = None
                    pos2 = None
                    # ALT already has precedence in the breakpoint helper.
                    # Read fallback metadata only if it is actually needed.
                    if not (alt and _BND_MATE.search(alt)):
                        chrom2 = record.info.get("CHR2") if "CHR2" in record.info else None
                        if "CHR2" in record.info or "POS2" in record.info:
                            pos2 = record.info.get("POS2", record.info.get("END"))
                    _set_translocation_breakpoints(variant, alt, chrom2, pos2)
                genotype = _pysam_first_genotype(record)
                if genotype is not None:
                    variant["genotype"] = genotype
                variants.append(variant)
        return variants

    except ImportError:
        logger.warning("pysam not available; falling back to basic VCF parsing")
        return _parse_vcf_basic(vcf_path)


def _pysam_first_genotype(record: Any) -> str | None:
    """Genotype of the first sample of a pysam record, or None when it has no ``GT``."""
    try:
        return _format_genotype(record.samples[0]["GT"])
    except (IndexError, KeyError, TypeError, ValueError):
        return None


#: The mate of a BND ALT: ``N[chr2:500[``, ``N]chr2:500]``, ``]chr2:500]N`` or ``[chr2:500[N`` (1-based position).
_BND_MATE = re.compile(r"[\[\]]([^\[\]]+):(\d+)[\[\]]")


def _breakend_partner(alt: str | None, chrom2: str | None, position2: int | None) -> tuple[str, int] | None:
    """Second breakpoint of a translocation record, as ``(chromosome, 0-based position)``, or None.

    The mate named in a BND ALT wins; else ``CHR2`` with its 1-based position (``POS2``, else ``END``). A record
    with neither (a single breakend such as ``N.``) has no partner.
    """
    if alt:
        mate = _BND_MATE.search(alt)
        if mate:
            return mate.group(1), int(mate.group(2)) - 1
    if chrom2 and position2 is not None:
        return str(chrom2), int(position2) - 1
    return None


def _set_translocation_breakpoints(variant: dict[str, Any], alt: str | None, chrom2: Any, position2: Any) -> None:
    """Give a TRA/BND record the package's translocation shape: a zero-length junction plus ``chrom2``/``pos2``.

    ``end`` repeats ``start`` as in :class:`~metainformant.structural_variants.detection.sv_calling.StructuralVariant`
    (a VCF ``END`` of such a record is a position on ``CHR2``, not the end of an interval). The second breakpoint is
    added when the record names one.
    """
    variant["end"] = variant["start"]
    partner = _breakend_partner(alt, chrom2, int(position2) if position2 is not None else None)
    if partner is not None:
        variant["chrom2"], variant["pos2"] = partner


def _parse_vcf_basic(vcf_path: str) -> list[dict[str, Any]]:
    """Basic VCF parser without pysam dependency.

    Args:
        vcf_path: Path to VCF file.

    Returns:
        List of variant dictionaries. A record with a sample column also carries
        the ``genotype`` of its first sample. A translocation (``SVTYPE`` TRA or
        BND) is a zero-length junction (``end == start``) with its second
        breakpoint in ``chrom2``/``pos2`` when the record names one: the mate of a
        BND ALT, or ``CHR2`` with ``POS2``/``END`` in INFO.
    """
    variants: list[dict[str, Any]] = []

    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = int(fields[1]) - 1  # Convert to 0-based
            qual = float(fields[5]) if fields[5] != "." else 0.0
            info = fields[7]

            # Parse INFO field
            info_dict: dict[str, str] = {}
            for item in info.split(";"):
                if "=" in item:
                    key, val = item.split("=", 1)
                    info_dict[key] = val
                else:
                    info_dict[item] = ""

            sv_type = info_dict.get("SVTYPE", "UNKNOWN")
            end = int(info_dict.get("END", str(pos + 1)))

            variant: dict[str, Any] = {
                "chrom": chrom,
                "start": pos,
                "end": end,
                "sv_type": sv_type,
                "quality": qual,
            }
            if sv_type in TRANSLOCATION_TYPES:
                position2 = info_dict.get("POS2", info_dict.get("END"))
                _set_translocation_breakpoints(
                    variant, fields[4], info_dict.get("CHR2"), int(position2) if position2 else None
                )
            genotype = _vcf_sample_genotype(fields)
            if genotype is not None:
                variant["genotype"] = genotype
            variants.append(variant)

    return variants
