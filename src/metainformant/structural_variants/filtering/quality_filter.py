"""Quality-based filtering for structural variants.

Provides filters based on variant quality scores, read support, size
constraints, population allele frequency, and blacklist genomic regions.
"""

from __future__ import annotations

import math
import numbers
import os
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, NamedTuple

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "SV_"


@dataclass
class FilterStats:
    """Statistics from a filtering operation.

    Attributes:
        input_count: Number of variants before filtering.
        output_count: Number of variants passing filter.
        filtered_count: Number of variants removed.
        filter_name: Name of the filter applied.
        parameters: Filter parameters used.
    """

    input_count: int
    output_count: int
    filtered_count: int
    filter_name: str
    parameters: dict[str, Any]

    @property
    def pass_rate(self) -> float:
        """Fraction of variants passing the filter."""
        if self.input_count == 0:
            return 0.0
        return self.output_count / self.input_count


def filter_by_quality(
    variants: list[dict[str, Any]],
    min_qual: float = 20.0,
    min_support: int = 3,
    min_mapq: float = 0.0,
    min_gq: float = 0.0,
) -> tuple[list[dict[str, Any]], FilterStats]:
    """Filter structural variants by quality metrics.

    Removes variants that do not meet minimum quality thresholds.
    Multiple quality metrics can be applied simultaneously.

    Args:
        variants: List of variant dictionaries. Expected keys:
            - 'quality' or 'qual': Quality score (float)
            - 'support' or 'n_support': Number of supporting reads (int)
            - 'mapq': Mean mapping quality of supporting reads (float, optional)
            - 'gq': Genotype quality (float, optional)
        min_qual: Minimum quality score (Phred-scaled, default 20).
        min_support: Minimum number of supporting reads (default 3).
        min_mapq: Minimum mean mapping quality (default 0, no filter).
        min_gq: Minimum genotype quality (default 0, no filter).

    Returns:
        Tuple of (filtered_variants, filter_statistics).
    """
    # Environment overrides
    env_qual = os.environ.get(f"{_ENV_PREFIX}MIN_QUAL")
    env_support = os.environ.get(f"{_ENV_PREFIX}MIN_SUPPORT")
    if env_qual is not None:
        min_qual = float(env_qual)
    if env_support is not None:
        min_support = int(env_support)

    input_count = len(variants)
    passed: list[dict[str, Any]] = []

    for v in variants:
        qual = v.get("quality", v.get("qual", 0.0))
        support = v.get("support", v.get("n_support", 0))

        # If variant has an evidence object, extract support from it
        evidence = v.get("evidence")
        if evidence is not None and support == 0:
            if hasattr(evidence, "split_reads"):
                support = evidence.split_reads + evidence.discordant_pairs
            elif isinstance(evidence, dict):
                support = evidence.get("split_reads", 0) + evidence.get("discordant_pairs", 0)

        mapq = v.get("mapq", 60.0)
        gq = v.get("gq", 99.0)

        if qual < min_qual:
            continue
        if support < min_support:
            continue
        if mapq < min_mapq:
            continue
        if gq < min_gq:
            continue

        passed.append(v)

    stats = FilterStats(
        input_count=input_count,
        output_count=len(passed),
        filtered_count=input_count - len(passed),
        filter_name="quality_filter",
        parameters={
            "min_qual": min_qual,
            "min_support": min_support,
            "min_mapq": min_mapq,
            "min_gq": min_gq,
        },
    )

    logger.info(
        f"Quality filter: {stats.output_count}/{stats.input_count} passed "
        f"(removed {stats.filtered_count}, pass rate {stats.pass_rate:.1%})"
    )
    return passed, stats


def filter_by_size(
    variants: list[dict[str, Any]],
    min_size: int = 50,
    max_size: int | None = None,
) -> tuple[list[dict[str, Any]], FilterStats]:
    """Filter structural variants by size.

    Removes variants outside the specified size range. This is commonly
    used to separate SVs (>=50bp) from indels (<50bp) and to exclude
    extremely large artifacts.

    Args:
        variants: List of variant dictionaries. Expected keys:
            - 'size': Variant size in base pairs (int)
            - Or 'start' and 'end' to compute size
            - 'sv_type': Optional; translocations (TRA, BND) are always kept
              since they don't have a meaningful "size", and so are insertions
              (INS) without a positive 'size': the inserted length of a call from
              split reads is not measured (``size`` 0), and the start-to-end span
              of such a call is its junction, not that length. An insertion with a
              positive 'size' is filtered on it like any other variant.
        min_size: Minimum size in base pairs (default 50, the standard
            SV size threshold).
        max_size: Maximum size in base pairs (default None, no upper limit).

    Returns:
        Tuple of (filtered_variants, filter_statistics).
    """
    input_count = len(variants)
    passed: list[dict[str, Any]] = []

    for v in variants:
        sv_type = v.get("sv_type", "")
        if hasattr(sv_type, "value"):
            sv_type = sv_type.value

        # Translocations don't have a standard size
        if sv_type in ("TRA", "BND"):
            passed.append(v)
            continue

        size = v.get("size", 0)

        # An insertion whose length was not measured has no size to filter on
        if sv_type == "INS" and not size:
            passed.append(v)
            continue

        if size == 0:
            start = v.get("start", 0)
            end = v.get("end", 0)
            size = abs(end - start)

        if size < min_size:
            continue
        if max_size is not None and size > max_size:
            continue

        passed.append(v)

    stats = FilterStats(
        input_count=input_count,
        output_count=len(passed),
        filtered_count=input_count - len(passed),
        filter_name="size_filter",
        parameters={"min_size": min_size, "max_size": max_size},
    )

    logger.info(
        f"Size filter: {stats.output_count}/{stats.input_count} passed "
        f"(range: {min_size}-{max_size if max_size else 'inf'} bp)"
    )
    return passed, stats


def filter_by_frequency(
    variants: list[dict[str, Any]],
    population_db: dict[str, float] | list[dict[str, Any]] | None = None,
    max_af: float = 0.01,
    match_window: int = 200,
    min_reciprocal_overlap: float = 0.5,
) -> tuple[list[dict[str, Any]], FilterStats]:
    """Filter structural variants by population allele frequency.

    Removes common variants that are likely benign based on their
    frequency in population databases (e.g., gnomAD-SV, DGV).

    Variant matching uses both position proximity and reciprocal overlap
    to find matching database entries.

    Args:
        variants: List of variant dictionaries with 'chrom', 'start', 'end',
            'sv_type' keys. May also have 'af' or 'allele_frequency' pre-computed.
        population_db: Population frequency database. Can be:
            - dict mapping 'chrom:start-end' to allele frequency
            - list of dicts with 'chrom', 'start', 'end', 'af' keys
            - None (no filtering applied, all variants pass)
        max_af: Maximum population allele frequency to keep (default 0.01 = 1%).
        match_window: Window in bp around breakpoints for matching database entries.
        min_reciprocal_overlap: Minimum reciprocal overlap fraction for matching.

    Returns:
        Tuple of (filtered_variants, filter_statistics).

    Raises:
        ValueError: If ``population_db`` is given but none of its entries is in a
            supported layout (so the filter would silently do nothing). Entries
            of a mostly usable database that cannot be read are skipped with a
            logged warning that counts them.
    """
    input_count = len(variants)

    # Check if any variants have pre-computed AF values
    has_precomputed_af = any(v.get("af") is not None or v.get("allele_frequency") is not None for v in variants)

    if population_db is None and not has_precomputed_af:
        stats = FilterStats(
            input_count=input_count,
            output_count=input_count,
            filtered_count=0,
            filter_name="frequency_filter",
            parameters={"max_af": max_af, "population_db": "none"},
        )
        logger.info("No population database provided; skipping frequency filter")
        return list(variants), stats

    # Normalize population database to a lookup structure
    pop_lookup: dict[str, list[tuple[int, int, float, str]]] = {}
    if population_db is not None:
        parsed = _parse_population_db(population_db)
        if parsed.n_records and not parsed.lookup:
            raise ValueError(
                f"population_db layout is not supported: none of its {parsed.n_records} entries could be used "
                f"(for example {parsed.examples[0]}). {_POPULATION_DB_LAYOUTS}"
            )
        if parsed.n_skipped:
            logger.warning(
                f"Frequency filter: skipped {parsed.n_skipped} of {parsed.n_records} population_db entries "
                f"in an unsupported layout (for example {parsed.examples[0]}). {_POPULATION_DB_LAYOUTS}"
            )
        pop_lookup = parsed.lookup

    passed: list[dict[str, Any]] = []

    for v in variants:
        # Check if AF is pre-computed on the variant
        af = v.get("af", v.get("allele_frequency", None))

        if af is None and pop_lookup:
            # Look up in population database
            chrom = v.get("chrom", "")
            start = v.get("start", 0)
            end = v.get("end", 0)
            sv_type = v.get("sv_type", "")
            if hasattr(sv_type, "value"):
                sv_type = sv_type.value

            af = _lookup_population_frequency(
                pop_lookup, chrom, start, end, sv_type, match_window, min_reciprocal_overlap
            )

        if af is not None and af > max_af:
            continue

        # Add the looked-up AF to the variant for downstream use
        result = dict(v)
        if af is not None:
            result["population_af"] = af
        passed.append(result)

    stats = FilterStats(
        input_count=input_count,
        output_count=len(passed),
        filtered_count=input_count - len(passed),
        filter_name="frequency_filter",
        parameters={"max_af": max_af, "match_window": match_window},
    )

    logger.info(f"Frequency filter (AF <= {max_af}): {stats.output_count}/{stats.input_count} passed")
    return passed, stats


def apply_blacklist(
    variants: list[dict[str, Any]],
    blacklist_regions: list[dict[str, Any]] | list[tuple[str, int, int]],
    min_overlap_fraction: float = 0.0,
) -> tuple[list[dict[str, Any]], FilterStats]:
    """Remove variants overlapping blacklist regions.

    Blacklist regions are genomic intervals known to produce false-positive
    SV calls, such as centromeres, telomeres, segmental duplications,
    low-complexity regions, and assembly gaps.

    Args:
        variants: List of variant dictionaries with 'chrom', 'start', 'end'.
        blacklist_regions: Blacklist as a list of dicts with 'chrom', 'start',
            'end' keys, or as (chrom, start, end) tuples.
        min_overlap_fraction: Minimum fraction of variant overlapping a blacklist
            region to trigger removal (default 0.0, any overlap removes the variant).
            A zero-length variant (``start == end``, such as an insertion point)
            has no bases to overlap: it counts as fully overlapping (fraction 1.0)
            when it lies strictly inside a region, and not at all on a region's
            edge, the same as a variant that only abuts a region.

    Returns:
        Tuple of (filtered_variants, filter_statistics).
    """
    input_count = len(variants)

    if not blacklist_regions:
        stats = FilterStats(
            input_count=input_count,
            output_count=input_count,
            filtered_count=0,
            filter_name="blacklist_filter",
            parameters={"n_blacklist_regions": 0},
        )
        return list(variants), stats

    # Normalize blacklist to dicts and group by chromosome
    bl_by_chrom: dict[str, list[tuple[int, int]]] = {}
    for region in blacklist_regions:
        if isinstance(region, tuple):
            chrom, start, end = region
        else:
            chrom = region.get("chrom", "")
            start = region.get("start", 0)
            end = region.get("end", 0)

        if chrom not in bl_by_chrom:
            bl_by_chrom[chrom] = []
        bl_by_chrom[chrom].append((start, end))

    # Sort blacklist regions for each chromosome
    for chrom in bl_by_chrom:
        bl_by_chrom[chrom].sort()

    passed: list[dict[str, Any]] = []

    for v in variants:
        chrom = v.get("chrom", "")
        start = v.get("start", 0)
        end = v.get("end", 0)
        variant_size = max(1, end - start)

        if chrom not in bl_by_chrom:
            passed.append(v)
            continue

        # Check overlap with any blacklist region
        max_overlap_frac = 0.0
        for bl_start, bl_end in bl_by_chrom[chrom]:
            if bl_start >= end:
                break  # No more possible overlaps (sorted)
            if bl_end <= start:
                continue

            if end == start:
                # Zero-length variant strictly inside the region (the two checks above
                # ruled out a region that begins at or after it, or ends at or before it).
                max_overlap_frac = 1.0
                break

            overlap_start = max(start, bl_start)
            overlap_end = min(end, bl_end)
            overlap_bp = max(0, overlap_end - overlap_start)
            overlap_frac = overlap_bp / variant_size
            max_overlap_frac = max(max_overlap_frac, overlap_frac)

        if max_overlap_frac > min_overlap_fraction:
            continue

        passed.append(v)

    stats = FilterStats(
        input_count=input_count,
        output_count=len(passed),
        filtered_count=input_count - len(passed),
        filter_name="blacklist_filter",
        parameters={
            "n_blacklist_regions": sum(len(v) for v in bl_by_chrom.values()),
            "min_overlap_fraction": min_overlap_fraction,
        },
    )

    logger.info(
        f"Blacklist filter: {stats.output_count}/{stats.input_count} passed "
        f"({stats.filtered_count} removed in blacklist regions)"
    )
    return passed, stats


_POPULATION_DB_LAYOUTS = (
    "Supported layouts: a dict mapping 'chrom:start-end' (or 'chrom:start-end-svtype') to an allele "
    "frequency, or a list of dicts with 'chrom', 'start', 'end' and 'af' (or 'allele_frequency') "
    "and optionally 'sv_type'."
)


class _PopulationParse(NamedTuple):
    """Outcome of reading a population database (internal).

    Attributes:
        lookup: Chromosome -> sorted ``(start, end, af, sv_type)`` entries that could be read.
        n_records: Number of entries in the database.
        n_skipped: Entries that could not be read.
        examples: Descriptions of the first few skipped entries.
    """

    lookup: dict[str, list[tuple[int, int, float, str]]]
    n_records: int
    n_skipped: int
    examples: list[str]


def _as_position(value: Any) -> int | None:
    """Return ``value`` as an int coordinate, or None when it is not an integral number."""
    if isinstance(value, bool):
        return None
    if isinstance(value, numbers.Integral):
        return int(value)
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return None


def _as_frequency(value: Any) -> float | None:
    """Return ``value`` as an allele frequency in [0, 1], or None when it is not one."""
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        return None
    frequency = float(value)
    if not math.isfinite(frequency) or not 0.0 <= frequency <= 1.0:
        return None
    return frequency


def _parse_population_db(population_db: Any) -> _PopulationParse:
    """Read a population database, counting the entries that are not in a supported layout.

    Args:
        population_db: A dict ``{'chrom:start-end[-svtype]': af}`` or an iterable of dicts
            with ``chrom``, ``start``, ``end``, ``af`` (or ``allele_frequency``) and
            optionally ``sv_type``.

    Returns:
        The usable entries by chromosome, plus how many entries were skipped and why.

    Raises:
        ValueError: If ``population_db`` is neither a dict nor an iterable.
    """
    lookup: dict[str, list[tuple[int, int, float, str]]] = {}
    n_records = 0
    n_skipped = 0
    examples: list[str] = []

    def skip(description: str) -> None:
        nonlocal n_skipped
        n_skipped += 1
        if len(examples) < 3:
            examples.append(description)

    if isinstance(population_db, Mapping):
        for key, af in population_db.items():
            n_records += 1
            parts = key.replace(":", "-").split("-") if isinstance(key, str) else []
            start = end = None
            if len(parts) >= 3 and parts[0]:
                try:
                    start, end = int(parts[1]), int(parts[2])
                except ValueError:
                    pass
            frequency = _as_frequency(af)
            if start is None or end is None or frequency is None:
                skip(f"key {key!r}: {af!r}")
                continue
            lookup.setdefault(parts[0], []).append((start, end, frequency, parts[3] if len(parts) > 3 else ""))
    else:
        try:
            entries = list(population_db)
        except TypeError:
            raise ValueError(
                f"population_db must be a dict or a list of dicts, not {type(population_db).__name__}. "
                f"{_POPULATION_DB_LAYOUTS}"
            ) from None
        for index, entry in enumerate(entries):
            n_records += 1
            if not isinstance(entry, Mapping):
                skip(f"entry {index}: {type(entry).__name__}, not a dict")
                continue
            chrom = entry.get("chrom")
            start = _as_position(entry.get("start"))
            end = _as_position(entry.get("end"))
            frequency = _as_frequency(entry.get("af", entry.get("allele_frequency")))
            if not isinstance(chrom, str) or not chrom or start is None or end is None or frequency is None:
                skip(f"entry {index}: {sorted(str(k) for k in entry)}")
                continue
            sv_type = entry.get("sv_type", "")
            if hasattr(sv_type, "value"):
                sv_type = sv_type.value
            lookup.setdefault(chrom, []).append((start, end, frequency, sv_type))

    # Sort each chromosome
    for chrom in lookup:
        lookup[chrom].sort()

    return _PopulationParse(lookup, n_records, n_skipped, examples)


def _build_population_lookup(
    population_db: dict[str, float] | list[dict[str, Any]],
) -> dict[str, list[tuple[int, int, float, str]]]:
    """Build a chromosome-keyed lookup from population database.

    Entries that are not in a supported layout are left out (use
    :func:`_parse_population_db` to count them).

    Args:
        population_db: Population database in either dict or list format.

    Returns:
        Dictionary mapping chromosome to sorted list of (start, end, af, sv_type) tuples.
    """
    return _parse_population_db(population_db).lookup


def _lookup_population_frequency(
    pop_lookup: dict[str, list[tuple[int, int, float, str]]],
    chrom: str,
    start: int,
    end: int,
    sv_type: str,
    match_window: int,
    min_ro: float,
) -> float | None:
    """Look up population frequency for a variant.

    Args:
        pop_lookup: Population lookup structure.
        chrom: Variant chromosome.
        start: Variant start.
        end: Variant end.
        sv_type: Variant SV type.
        match_window: Position matching window.
        min_ro: Minimum reciprocal overlap.

    Returns:
        Allele frequency if a match is found, None otherwise.
    """
    if chrom not in pop_lookup:
        return None

    variant_size = max(1, end - start)
    best_af: float | None = None

    for p_start, p_end, p_af, p_type in pop_lookup[chrom]:
        if p_start > end + match_window:
            break
        if p_end < start - match_window:
            continue

        # Check SV type match (if specified)
        if p_type and sv_type and p_type != sv_type:
            continue

        # Check reciprocal overlap
        overlap_start = max(start, p_start)
        overlap_end = min(end, p_end)
        overlap_bp = max(0, overlap_end - overlap_start)

        pop_size = max(1, p_end - p_start)
        ro_variant = overlap_bp / variant_size
        ro_pop = overlap_bp / pop_size

        if ro_variant >= min_ro and ro_pop >= min_ro:
            if best_af is None or p_af > best_af:
                best_af = p_af

    return best_af
