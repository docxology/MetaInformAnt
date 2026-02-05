"""Multi-caller merging and consensus for structural variants.

Implements algorithms for merging SV callsets from multiple detection
tools, computing reciprocal overlap between SVs, SURVIVOR-style
distance-based merging, and deduplication of variant calls.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]

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
        caller_names: Names of callers that detected this variant.
        support_variants: Original variant calls that were merged.
        consensus_quality: Quality based on multi-caller agreement.
        genotype: Consensus genotype.
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
) -> tuple[list[MergedVariant], MergeStats]:
    """Merge structural variant callsets from multiple callers.

    Uses reciprocal overlap to identify the same SV detected by different
    tools. For each group of matching variants across callers, produces a
    single consensus call with aggregated support.

    Args:
        callsets: Dictionary mapping caller names to lists of variant
            dictionaries. Each variant should have 'chrom', 'start',
            'end', 'sv_type' keys.
        min_overlap: Minimum reciprocal overlap fraction (0-1) to consider
            two variants as the same event (default 0.5 = 50%).
        type_match: Whether SV types must match for merging (default True).

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
        chrom_i = v_i.get("chrom", "")
        start_i = v_i.get("start", 0)
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

            end_j = v_j.get("end", 0)
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
        callers = list(set(caller for caller, _ in group_variants))
        variant_dicts = [v for _, v in group_variants]

        # Consensus position: median of starts and ends
        starts = [v.get("start", 0) for v in variant_dicts]
        ends = [v.get("end", 0) for v in variant_dicts]

        if np is not None:
            consensus_start = int(np.median(starts))
            consensus_end = int(np.median(ends))
        else:
            starts.sort()
            ends.sort()
            consensus_start = starts[len(starts) // 2]
            consensus_end = ends[len(ends) // 2]

        # Consensus SV type: majority vote
        type_counts: dict[str, int] = {}
        for v in variant_dicts:
            st = v.get("sv_type", "UNKNOWN")
            if hasattr(st, "value"):
                st = st.value
            type_counts[st] = type_counts.get(st, 0) + 1
        consensus_type = max(type_counts, key=type_counts.get)  # type: ignore[arg-type]

        size = abs(consensus_end - consensus_start)

        # Quality from multi-caller support
        n_callers = len(callers)
        base_qual = sum(v.get("quality", v.get("qual", 0.0)) for v in variant_dicts) / len(variant_dicts)
        consensus_quality = base_qual + 10.0 * math.log10(n_callers + 1) * 10

        # Consensus genotype
        genotypes = [v.get("genotype", "./.") for v in variant_dicts]
        consensus_gt = _consensus_genotype(genotypes)

        merged = MergedVariant(
            chrom=variant_dicts[0].get("chrom", ""),
            start=consensus_start,
            end=consensus_end,
            sv_type=consensus_type,
            size=size,
            n_callers=n_callers,
            caller_names=callers,
            support_variants=variant_dicts,
            consensus_quality=consensus_quality,
            genotype=consensus_gt,
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

    For inter-chromosomal SVs or SVs without clear coordinates,
    uses breakpoint distance as a proxy.

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

    start1 = sv1.get("start", 0)
    end1 = sv1.get("end", 0)
    start2 = sv2.get("start", 0)
    end2 = sv2.get("end", 0)

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

    Args:
        vcf_files: Either file paths to VCF files or lists of variant
            dictionaries. If file paths, each file is treated as one caller.
            If dicts, each dict should have a '_caller' key.
        max_distance: Maximum distance between breakpoints to consider
            variants as the same event (default 1000bp).
        min_callers: Minimum number of callers required to keep a merged
            variant (default 2).
        type_match: Whether SV types must match (default True).
        strand_match: Whether strand orientations must match (default False).

    Returns:
        Tuple of (merged_variants, merge_statistics).
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
            caller_name = item.get("_caller", f"caller_{i}")
            if caller_name not in callsets:
                callsets[caller_name] = []
            callsets[caller_name].append(item)
        elif isinstance(item, list):
            caller_name = f"caller_{i}"
            callsets[caller_name] = item

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
        callers = list(set(caller for caller, _ in group_variants))

        if len(callers) < min_callers:
            continue

        variant_dicts = [v for _, v in group_variants]

        # Consensus position
        starts = [v.get("start", 0) for v in variant_dicts]
        ends = [v.get("end", 0) for v in variant_dicts]

        if np is not None:
            consensus_start = int(np.median(starts))
            consensus_end = int(np.median(ends))
        else:
            starts.sort()
            ends.sort()
            consensus_start = starts[len(starts) // 2]
            consensus_end = ends[len(ends) // 2]

        # Consensus type
        type_counts: dict[str, int] = {}
        for v in variant_dicts:
            st = v.get("sv_type", "UNKNOWN")
            if hasattr(st, "value"):
                st = st.value
            type_counts[st] = type_counts.get(st, 0) + 1
        consensus_type = max(type_counts, key=type_counts.get)  # type: ignore[arg-type]

        size = abs(consensus_end - consensus_start)
        n_callers = len(callers)

        quality = 10.0 * n_callers * math.log10(n_callers + 1)

        merged = MergedVariant(
            chrom=variant_dicts[0].get("chrom", ""),
            start=consensus_start,
            end=consensus_end,
            sv_type=consensus_type,
            size=size,
            n_callers=n_callers,
            caller_names=callers,
            support_variants=variant_dicts,
            consensus_quality=quality,
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

    Args:
        variants: List of variant dictionaries with 'chrom', 'start', 'end',
            'sv_type', and optionally 'quality' keys.
        max_distance: Maximum breakpoint distance for duplicate candidates.
        min_reciprocal_overlap: Minimum reciprocal overlap to classify as
            duplicate (default 0.8 = 80%).
        type_match: Whether SV types must match (default True).

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
        best_idx = max(
            indices,
            key=lambda i: sorted_vars[i].get("quality", sorted_vars[i].get("qual", 0.0)),
        )
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


def _consensus_genotype(genotypes: list[str]) -> str:
    """Determine consensus genotype from multiple calls.

    Uses majority voting among non-missing genotypes.

    Args:
        genotypes: List of genotype strings.

    Returns:
        Consensus genotype string.
    """
    valid = [g for g in genotypes if g not in ("./.", ".")]
    if not valid:
        return "./."

    counts: dict[str, int] = {}
    for g in valid:
        # Normalize genotype ordering
        normalized = "/".join(sorted(g.split("/")))
        counts[normalized] = counts.get(normalized, 0) + 1

    return max(counts, key=counts.get)  # type: ignore[arg-type]


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
                variants.append(
                    {
                        "chrom": record.chrom,
                        "start": record.start,
                        "end": end,
                        "sv_type": sv_type,
                        "quality": record.qual or 0.0,
                        "filter": list(record.filter),
                    }
                )
        return variants

    except ImportError:
        logger.warning("pysam not available; falling back to basic VCF parsing")
        return _parse_vcf_basic(vcf_path)


def _parse_vcf_basic(vcf_path: str) -> list[dict[str, Any]]:
    """Basic VCF parser without pysam dependency.

    Args:
        vcf_path: Path to VCF file.

    Returns:
        List of variant dictionaries.
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

            variants.append(
                {
                    "chrom": chrom,
                    "start": pos,
                    "end": end,
                    "sv_type": sv_type,
                    "quality": qual,
                }
            )

    return variants
