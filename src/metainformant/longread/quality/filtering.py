"""Read filtering and adapter trimming for long-read sequencing data.

Provides length-based filtering, quality filtering, adapter detection and
trimming, and chimeric read detection/splitting. All algorithms are real
implementations with proper bioinformatics logic.

Adapter matching model
----------------------
An adapter is located by semi-global alignment: the whole adapter must align to
a stretch of the read, the stretch may start and end anywhere, and every
substitution, insertion or deletion costs one edit. A placement is accepted when
it needs at most ``max_errors`` edits, the smallest of

* ``floor(len(adapter) * (1 - min_identity))``, the divergence the caller tolerates;
* the largest budget whose false-positive probability stays within
  ``max_false_positive_rate`` (:func:`adapter_false_positive_bound`);
* ``len(adapter) // 7 - 1``, which keeps the exact seeds that propose candidate
  placements at least 7 bp long. Every placement within the budget contains one
  of the ``max_errors + 1`` seeds exactly (pigeonhole), so the seed step loses
  nothing the alignment model would accept, and adapters are found at any offset.

False-positive bound: for i.i.d. uniform bases, the probability that a given read
offset ends a spurious placement of a ``length``-bp adapter within ``max_errors``
edits is at most :func:`adapter_false_positive_bound` (a union bound over edit
scripts, so it is conservative). A read of ``n`` adapter-free bases therefore has
at most ``2 * n * bound`` expected spurious hits per adapter, counting both
orientations. With the default rate of 1e-9 that is at most 2e-5 hits for a 10 kb
read, and the shortest adapter that can be searched is 15 bp
(:func:`minimum_adapter_length`). The bound assumes uniform base composition;
GC- or AT-rich sequence raises the true rate, so lower ``max_false_positive_rate``
for such genomes. An adapter shorter than the minimum is searched by exact match
only and a warning names it.

Optional dependencies:
    - numpy: For efficient numerical computation
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Sequence, TypeVar, cast

from metainformant.core.sequence import reverse_complement
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ReadT = TypeVar("_ReadT")
np: Any  # None when numpy is unavailable (optional dependency)
try:
    import numpy as np
except ImportError:
    np = None

#: Default upper bound on the probability that one read offset starts or ends a
#: spurious adapter placement in adapter-free random sequence.
DEFAULT_MAX_FALSE_POSITIVE_RATE = 1e-9

_MIN_SEED_LENGTH = 7
_DETECT_SEARCH_LENGTH = 300


# Standard ONT adapter sequences. Every entry must be long enough to be searched with
# a bounded false-positive rate (see minimum_adapter_length). The 8 bp ONT_LSK110_top
# and ONT_LSK110_bottom entries (CAGCACCT and AGGTGCTG, reverse complements of each
# other) were removed: they matched random sequence at about one offset in 240.
ONT_ADAPTERS: dict[str, str] = {
    "ONT_LSK109_top": "AATGTACTTCGTTCAGTTACGTATTGCT",
    "ONT_LSK109_bottom": "GCAATACGTAACTGAACGAAGT",
    "ONT_RNA_adapter": "GGCGTCTGCTTGGGTGTTTAACCT",
    "ONT_rapid_top": "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
    "ONT_rapid_barcode": "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA",
}

# Standard PacBio adapter sequences
PACBIO_ADAPTERS: dict[str, str] = {
    "PacBio_SMRTbell": "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT",
    "PacBio_barcoded": "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
}

ALL_KNOWN_ADAPTERS: dict[str, str] = {**ONT_ADAPTERS, **PACBIO_ADAPTERS}


@dataclass
class ReadRecord:
    """A simple read record for filtering operations.

    Attributes:
        read_id: Unique read identifier.
        sequence: Nucleotide sequence.
        quality_string: ASCII-encoded quality string (Phred+33).
        metadata: Additional key-value metadata.
    """

    read_id: str
    sequence: str
    quality_string: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class AdapterMatch:
    """Describes a detected adapter in a read sequence.

    Attributes:
        adapter_name: Name of the matched adapter.
        adapter_sequence: Adapter sequence that matched.
        position: Start position in the read (0-based).
        end_position: End position in the read (0-based, exclusive).
        score: Alignment score of the match (adapter length minus edits).
        identity: Fraction of the adapter that matched (1 - edit_distance / adapter length).
        location: Where the adapter was found ('start', 'end', 'internal').
        edit_distance: Substitutions, insertions and deletions in the placement.
        false_positive_bound: Upper bound on the probability that one offset of an
            adapter-free random read yields a placement with this adapter length
            and edit budget (see :func:`adapter_false_positive_bound`).
    """

    adapter_name: str
    adapter_sequence: str
    position: int
    end_position: int
    score: float
    identity: float
    location: str
    edit_distance: int = 0
    false_positive_bound: float = 0.0


def filter_by_length(
    reads: Sequence[_ReadT],
    min_length: int = 1000,
    max_length: int | None = None,
) -> list[_ReadT]:
    """Filter reads by sequence length.

    Args:
        reads: Sequence of read dictionaries (with 'sequence' key) or ReadRecord objects.
        min_length: Minimum read length to keep (inclusive). Default 1000.
        max_length: Maximum read length to keep (inclusive). None for no upper limit.

    Returns:
        List of reads passing the length filter.
    """
    filtered: list[_ReadT] = []
    total = 0
    passed = 0

    for read in reads:
        total += 1
        seq = _get_sequence(read)
        if seq is None:
            continue

        length = len(seq)
        if length < min_length:
            continue
        if max_length is not None and length > max_length:
            continue

        filtered.append(read)
        passed += 1

    logger.info(
        "Length filter: %d/%d reads passed (min=%d, max=%s)",
        passed,
        total,
        min_length,
        max_length or "None",
    )
    return filtered


def filter_by_quality(
    reads: Sequence[_ReadT],
    min_q: float = 7.0,
) -> list[_ReadT]:
    """Filter reads by mean Phred quality score.

    Calculates the mean quality for each read and retains only those
    with mean quality >= min_q.

    Args:
        reads: Sequence of read dictionaries or ReadRecord objects.
        min_q: Minimum mean Phred quality score to keep. Default 7.0.

    Returns:
        List of reads passing the quality filter.
    """
    filtered: list[_ReadT] = []
    total = 0
    passed = 0

    for read in reads:
        total += 1
        qual_str = _get_quality(read)
        if not qual_str:
            # If no quality info, include the read
            filtered.append(read)
            passed += 1
            continue

        mean_q = _mean_phred(qual_str)
        if mean_q >= min_q:
            filtered.append(read)
            passed += 1

    logger.info("Quality filter: %d/%d reads passed (min_q=%.1f)", passed, total, min_q)
    return filtered


@lru_cache(maxsize=None)
def adapter_false_positive_bound(adapter_length: int, max_errors: int) -> float:
    """Upper bound on the chance that a random read offset ends a spurious adapter placement.

    For an i.i.d. uniform ACGT read, the probability that the bases ending at a
    given offset are within ``max_errors`` edits (substitutions, insertions,
    deletions) of a fixed adapter of ``adapter_length`` bp. It is a union bound
    over edit scripts, summed over the lengths a placement can span::

        sum over (d, i, s) with d + i + s <= max_errors of
            C(L, d) * C(L - d, s) * 3**s * C(L - d + i, i) * 4**i / 4**(L - d + i)

    where d, i and s count deletions, insertions and substitutions. Different
    scripts can yield the same string, so the true probability is at most this value.

    Args:
        adapter_length: Adapter length L in bp (at least 1).
        max_errors: Edit budget (at least 0).

    Returns:
        A probability in (0, 1]. ``4 ** -adapter_length`` for ``max_errors == 0``.

    Raises:
        ValueError: If ``adapter_length`` is below 1 or ``max_errors`` is negative.
    """
    if adapter_length < 1:
        raise ValueError(f"adapter_length must be at least 1, got {adapter_length}")
    if max_errors < 0:
        raise ValueError(f"max_errors must be non-negative, got {max_errors}")

    total = 0.0
    for deletions in range(max_errors + 1):
        kept = adapter_length - deletions
        for insertions in range(max_errors + 1 - deletions):
            span = kept + insertions
            if span < 1:
                continue
            substitutions = sum(math.comb(kept, s) * 3**s for s in range(max_errors - deletions - insertions + 1))
            scripts = math.comb(adapter_length, deletions) * math.comb(span, insertions) * 4**insertions * substitutions
            total += scripts / 4**span
    return min(1.0, total)


def minimum_adapter_length(max_false_positive_rate: float = DEFAULT_MAX_FALSE_POSITIVE_RATE) -> int:
    """Shortest adapter, in bp, whose exact match has a false-positive probability within the rate.

    The smallest L with ``4 ** -L <= max_false_positive_rate`` (15 bp for the default 1e-9).

    Raises:
        ValueError: If ``max_false_positive_rate`` is outside (0, 1).
    """
    _validate_search_parameters(1.0, max_false_positive_rate)
    length = 1
    while 4.0**-length > max_false_positive_rate:
        length += 1
    return length


def adapter_error_budget(
    adapter_length: int,
    min_identity: float = 0.75,
    max_false_positive_rate: float = DEFAULT_MAX_FALSE_POSITIVE_RATE,
) -> int:
    """Edits tolerated when searching for an adapter of ``adapter_length`` bp.

    The smallest of the divergence implied by ``min_identity``, the largest budget
    whose :func:`adapter_false_positive_bound` is within ``max_false_positive_rate``,
    and ``adapter_length // 7 - 1`` (keeps the exact search seeds at least 7 bp).

    Returns:
        The edit budget, or -1 when even an exact match of this length exceeds the
        false-positive rate (the adapter is shorter than :func:`minimum_adapter_length`).

    Raises:
        ValueError: If ``adapter_length`` is below 1, ``min_identity`` is outside 0-1,
            or ``max_false_positive_rate`` is outside (0, 1).
    """
    _validate_search_parameters(min_identity, max_false_positive_rate)
    if adapter_length < 1:
        raise ValueError(f"adapter_length must be at least 1, got {adapter_length}")

    tolerated = math.floor(adapter_length * (1.0 - min_identity) + 1e-9)
    seed_limit = max(0, adapter_length // _MIN_SEED_LENGTH - 1)

    budget = -1
    for errors in range(min(tolerated, seed_limit) + 1):
        if adapter_false_positive_bound(adapter_length, errors) > max_false_positive_rate:
            break
        budget = errors
    return budget


def trim_adapters(
    reads: Sequence[dict[str, Any] | ReadRecord],
    adapter_sequences: dict[str, str] | None = None,
    max_search_length: int = 200,
    min_identity: float = 0.75,
    max_false_positive_rate: float = DEFAULT_MAX_FALSE_POSITIVE_RATE,
) -> list[ReadRecord]:
    """Trim adapter sequences from the ends of reads.

    Searches the first and last ``max_search_length`` bases of each read for
    adapters by semi-global alignment (see the module docstring for the error
    model) and trims everything up to the end of an adapter found at the start,
    or from the start of an adapter found at the end. Adapter-free reads are
    returned unchanged: the false-positive probability per searched offset is
    at most ``max_false_positive_rate`` for adapters of at least
    :func:`minimum_adapter_length` bases.

    Args:
        reads: Sequence of read dictionaries or ReadRecord objects.
        adapter_sequences: Dictionary of adapter_name -> sequence. If None,
            uses all known ONT and PacBio adapters. Adapters shorter than
            :func:`minimum_adapter_length` are matched exactly only, with a warning.
        max_search_length: Maximum number of bases to search at each end.
        min_identity: Minimum fraction of the adapter that must match (0-1). It
            caps the tolerated edits; the false-positive bound can lower the cap.
        max_false_positive_rate: Upper bound on the per-offset probability of a
            spurious placement in adapter-free i.i.d. uniform sequence.

    Returns:
        List of ReadRecord objects with adapters trimmed.

    Raises:
        ValueError: If ``min_identity`` is outside 0-1 or ``max_false_positive_rate``
            is outside (0, 1).
    """
    if adapter_sequences is None:
        adapter_sequences = ALL_KNOWN_ADAPTERS

    plan = _plan_adapters(adapter_sequences, min_identity, max_false_positive_rate)
    start_patterns = [pattern for pattern in plan if not pattern.is_reverse]

    trimmed: list[ReadRecord] = []

    for read in reads:
        seq = _get_sequence(read)
        qual_str = _get_quality(read) or ""
        read_id = _get_read_id(read)

        if seq is None:
            continue

        trim_start = 0
        trim_end = len(seq)

        # Check for adapters at the start (forward orientation only)
        search_start = seq[:max_search_length].upper()
        for pattern in start_patterns:
            match = _best_placement(search_start, pattern)
            if match is not None:
                adapter_end = match[1]
                trim_start = max(trim_start, adapter_end)

        # Check for adapters at the end, in both orientations
        search_end = seq[-max_search_length:].upper() if len(seq) > max_search_length else seq.upper()
        offset = len(seq) - len(search_end)
        for pattern in plan:
            match = _best_placement(search_end, pattern)
            if match is not None:
                adapter_start = offset + match[0]
                trim_end = min(trim_end, adapter_start)

        # Apply trimming
        trimmed_seq = seq[trim_start:trim_end]
        trimmed_qual = qual_str[trim_start:trim_end] if qual_str else ""

        if len(trimmed_seq) > 0:
            metadata = {}
            if isinstance(read, dict):
                metadata = {
                    k: v for k, v in read.items() if k not in ("sequence", "quality", "quality_string", "read_id")
                }
            elif hasattr(read, "metadata"):
                metadata = dict(read.metadata)

            metadata["original_length"] = len(seq)
            metadata["trimmed_start"] = trim_start
            metadata["trimmed_end"] = len(seq) - trim_end

            trimmed.append(
                ReadRecord(
                    read_id=read_id,
                    sequence=trimmed_seq,
                    quality_string=trimmed_qual,
                    metadata=metadata,
                )
            )

    logger.info("Adapter trimming: processed %d reads, %d with sequence remaining", len(reads), len(trimmed))
    return trimmed


def detect_adapters(
    sequence: str,
    known_adapters: dict[str, str] | None = None,
    max_search_length: int = _DETECT_SEARCH_LENGTH,
    min_identity: float = 0.70,
    max_false_positive_rate: float = DEFAULT_MAX_FALSE_POSITIVE_RATE,
) -> list[AdapterMatch]:
    """Detect adapter sequences within a read.

    Searches for known adapter sequences throughout the read, at every offset
    and in both orientations, to identify internal adapters that may indicate
    chimeric reads. An adapter with up to ``max_errors`` substitutions,
    insertions or deletions is found wherever it sits (see the module docstring
    for how ``max_errors`` follows from ``min_identity`` and
    ``max_false_positive_rate``). On an adapter-free i.i.d. uniform read of
    ``n`` bases the expected number of spurious matches is at most
    ``2 * n * false_positive_bound`` per adapter, where the bound is reported on
    every :class:`AdapterMatch`.

    Args:
        sequence: Nucleotide sequence to search.
        known_adapters: Dictionary of adapter_name -> sequence. If None,
            uses all known adapters. Adapters shorter than
            :func:`minimum_adapter_length` are matched exactly only, with a warning.
        max_search_length: Length of the start/end regions used to label a
            match ``'start'`` or ``'end'`` rather than ``'internal'``.
        min_identity: Minimum fraction of the adapter that must match (0-1). It
            caps the tolerated edits; the false-positive bound can lower the cap.
        max_false_positive_rate: Upper bound on the per-offset probability of a
            spurious placement in adapter-free i.i.d. uniform sequence.

    Returns:
        List of AdapterMatch objects for all detected adapters, sorted by
        position, with overlapping matches reduced to the best-scoring one.

    Raises:
        ValueError: If ``min_identity`` is outside 0-1 or ``max_false_positive_rate``
            is outside (0, 1).
    """
    if known_adapters is None:
        known_adapters = ALL_KNOWN_ADAPTERS

    plan = _plan_adapters(known_adapters, min_identity, max_false_positive_rate)

    if not sequence:
        return []

    return _detect_with_plan(sequence.upper(), plan, max_search_length)


def split_chimeric_reads(
    reads: Sequence[dict[str, Any] | ReadRecord],
    known_adapters: dict[str, str] | None = None,
    min_fragment_length: int = 500,
    min_adapter_identity: float = 0.75,
    max_false_positive_rate: float = DEFAULT_MAX_FALSE_POSITIVE_RATE,
) -> list[ReadRecord]:
    """Detect and split chimeric reads at internal adapter sequences.

    Chimeric reads are artifacts where two separate DNA molecules are
    joined during library preparation. They contain internal adapter
    sequences. This function finds internal adapters (see
    :func:`detect_adapters` for the matching model and its false-positive
    bound) and splits the read into separate fragments. Adapter-free reads are
    not split: each adapter has at most ``max_false_positive_rate`` probability
    per read offset of a spurious internal hit.

    Args:
        reads: Sequence of read dictionaries or ReadRecord objects.
        known_adapters: Adapter sequences to search for. If None, uses defaults.
            Adapters shorter than :func:`minimum_adapter_length` are matched
            exactly only, with a warning.
        min_fragment_length: Minimum length of resulting fragments to keep.
        min_adapter_identity: Minimum fraction of the adapter that must match (0-1).
        max_false_positive_rate: Upper bound on the per-offset probability of a
            spurious placement in adapter-free i.i.d. uniform sequence.

    Returns:
        List of ReadRecord objects, with chimeric reads split into fragments.
        Non-chimeric reads are returned unchanged.

    Raises:
        ValueError: If ``min_adapter_identity`` is outside 0-1 or
            ``max_false_positive_rate`` is outside (0, 1).
    """
    if known_adapters is None:
        known_adapters = ALL_KNOWN_ADAPTERS

    plan = _plan_adapters(known_adapters, min_adapter_identity, max_false_positive_rate)

    result: list[ReadRecord] = []
    chimeric_count = 0

    for read in reads:
        seq = _get_sequence(read)
        qual_str = _get_quality(read) or ""
        read_id = _get_read_id(read)

        if seq is None:
            logger.warning(f"split_chimeric_reads: read {read_id} has no sequence; skipping")
            continue
        if len(seq) < min_fragment_length:
            # Too short to be chimeric: pass through unchanged
            result.append(
                ReadRecord(
                    read_id=read_id,
                    sequence=seq,
                    quality_string=qual_str,
                    metadata={"chimeric": False},
                )
            )
            continue

        # Detect internal adapters
        adapter_matches = _detect_with_plan(seq.upper(), plan, _DETECT_SEARCH_LENGTH)
        internal_matches = [m for m in adapter_matches if m.location == "internal"]

        if not internal_matches:
            # Not chimeric, return as-is
            result.append(
                ReadRecord(
                    read_id=read_id,
                    sequence=seq,
                    quality_string=qual_str,
                    metadata={"chimeric": False},
                )
            )
            continue

        # Split at internal adapters
        chimeric_count += 1
        split_points: list[tuple[int, AdapterMatch | None]] = [(0, None)]
        for match in internal_matches:
            split_points.append((match.position, match))
            split_points.append((match.end_position, None))
        split_points.append((len(seq), None))

        # Generate fragments between adapters
        fragment_idx = 0
        i = 0
        while i < len(split_points) - 1:
            start = split_points[i][0]
            # Find the next split point that is not within an adapter
            end = split_points[i + 1][0]
            adapter = split_points[i][1]

            if adapter is not None:
                # This is the start of an adapter, skip to after it
                i += 1
                continue

            # This is a sequence region between adapters
            frag_seq = seq[start:end]
            frag_qual = qual_str[start:end] if qual_str else ""

            if len(frag_seq) >= min_fragment_length:
                result.append(
                    ReadRecord(
                        read_id=f"{read_id}_fragment_{fragment_idx}",
                        sequence=frag_seq,
                        quality_string=frag_qual,
                        metadata={
                            "chimeric": True,
                            "parent_read": read_id,
                            "fragment_index": fragment_idx,
                            "fragment_start": start,
                            "fragment_end": end,
                        },
                    )
                )
                fragment_idx += 1

            i += 1

    logger.info(
        "Chimeric detection: %d/%d reads were chimeric, produced %d total fragments",
        chimeric_count,
        len(reads),
        len(result),
    )
    return result


# --- Internal helper functions ---


def _get_sequence(read: Any) -> str | None:
    """Extract sequence from various read representations."""
    if isinstance(read, dict):
        return read.get("sequence")
    elif hasattr(read, "sequence"):
        return cast(str, read.sequence)
    elif hasattr(read, "query_sequence"):
        return cast(str, read.query_sequence)
    return None


def _get_quality(read: Any) -> str | None:
    """Extract quality string from various read representations."""
    if isinstance(read, dict):
        return read.get("quality") or read.get("quality_string")
    elif hasattr(read, "quality_string"):
        return cast(str, read.quality_string)
    return None


def _get_read_id(read: Any) -> str:
    """Extract read ID from various read representations."""
    if isinstance(read, dict):
        return str(read.get("read_id", "unknown"))
    elif hasattr(read, "read_id"):
        return str(read.read_id)
    elif hasattr(read, "read_name"):
        return str(read.read_name)
    return "unknown"


def _mean_phred(quality_string: str) -> float:
    """Calculate mean Phred quality from ASCII quality string."""
    if not quality_string:
        return 0.0
    scores = [ord(c) - 33 for c in quality_string]
    return sum(scores) / len(scores)


@dataclass(frozen=True)
class _AdapterPattern:
    """One oriented adapter sequence with its edit budget and exact seeds."""

    name: str
    sequence: str
    is_reverse: bool
    max_errors: int
    seeds: tuple[tuple[int, str], ...]
    false_positive_bound: float


def _validate_search_parameters(min_identity: float, max_false_positive_rate: float) -> None:
    """Reject identity and false-positive settings outside their valid ranges."""
    if not 0.0 <= min_identity <= 1.0:
        raise ValueError(f"min_identity must be between 0 and 1, got {min_identity}")
    if not 0.0 < max_false_positive_rate < 1.0:
        raise ValueError(f"max_false_positive_rate must be between 0 and 1 (exclusive), got {max_false_positive_rate}")


def _split_seeds(sequence: str, pieces: int) -> tuple[tuple[int, str], ...]:
    """Cut ``sequence`` into ``pieces`` contiguous seeds of near-equal length as (offset, seed)."""
    base, extra = divmod(len(sequence), pieces)
    seeds: list[tuple[int, str]] = []
    offset = 0
    for index in range(pieces):
        size = base + (1 if index < extra else 0)
        seeds.append((offset, sequence[offset : offset + size]))
        offset += size
    return tuple(seeds)


def _plan_adapters(
    adapters: dict[str, str],
    min_identity: float,
    max_false_positive_rate: float,
) -> list[_AdapterPattern]:
    """Turn an adapter table into oriented search patterns, each with its edit budget.

    Adapters too short for any false-positive-bounded search keep an exact-match
    budget and are reported once through a warning.

    Raises:
        ValueError: If ``min_identity`` or ``max_false_positive_rate`` is out of
            range, or an adapter contains characters that are not IUPAC DNA.
    """
    _validate_search_parameters(min_identity, max_false_positive_rate)

    plan: list[_AdapterPattern] = []
    for name, raw_sequence in adapters.items():
        adapter = raw_sequence.upper()
        if not adapter:
            continue

        max_errors = adapter_error_budget(len(adapter), min_identity, max_false_positive_rate)
        if max_errors < 0:
            logger.warning(
                "Adapter %s is %d bp but %d bp are needed to keep spurious matches below %.1e per read offset; "
                "searching exact matches only (spurious-match probability %.1e per offset)",
                name,
                len(adapter),
                minimum_adapter_length(max_false_positive_rate),
                max_false_positive_rate,
                adapter_false_positive_bound(len(adapter), 0),
            )
            max_errors = 0

        bound = adapter_false_positive_bound(len(adapter), max_errors)
        oriented = [(name, adapter, False)]
        reverse = reverse_complement(adapter)
        if reverse != adapter:
            oriented.append((f"{name}_rc", reverse, True))

        for pattern_name, sequence, is_reverse in oriented:
            plan.append(
                _AdapterPattern(
                    name=pattern_name,
                    sequence=sequence,
                    is_reverse=is_reverse,
                    max_errors=max_errors,
                    seeds=_split_seeds(sequence, max_errors + 1),
                    false_positive_bound=bound,
                )
            )

    return plan


def _align_adapter(adapter: str, text: str, max_errors: int) -> tuple[int, int, int] | None:
    """Place the whole ``adapter`` inside ``text`` with the fewest edits.

    Semi-global (Sellers) alignment: the adapter is aligned end to end, the
    stretch of ``text`` it covers may start and end anywhere, and substitutions,
    insertions and deletions each cost one edit.

    Returns:
        ``(start, end, edits)`` of the best placement (leftmost end on ties), or
        None when every placement needs more than ``max_errors`` edits.
    """
    rows: list[list[int]] = [[0] * (len(text) + 1)]
    for row_index, base in enumerate(adapter, start=1):
        previous = rows[-1]
        current = [row_index]
        for column, other in enumerate(text, start=1):
            current.append(min(previous[column - 1] + (base != other), previous[column] + 1, current[column - 1] + 1))
        rows.append(current)

    last_row = rows[-1]
    edits = min(last_row)
    if edits > max_errors:
        return None

    end = last_row.index(edits)
    row_index, column = len(adapter), end
    while row_index > 0:
        here = rows[row_index][column]
        if column > 0 and here == rows[row_index - 1][column - 1] + (adapter[row_index - 1] != text[column - 1]):
            row_index, column = row_index - 1, column - 1
        elif here == rows[row_index - 1][column] + 1:
            row_index -= 1
        else:
            column -= 1
    return column, end, edits


def _find_adapter_placements(text: str, pattern: _AdapterPattern) -> list[tuple[int, int, int]]:
    """Find every placement of ``pattern`` in upper-case ``text`` within its edit budget.

    Exact seeds propose candidate regions and the alignment confirms them. A
    placement with at most ``max_errors`` edits always contains one of the
    ``max_errors + 1`` seeds unchanged, so no accepted placement is missed.

    Returns:
        ``(start, end, edits)`` tuples; overlapping placements of one adapter
        occurrence may appear more than once and are reduced by the caller.
    """
    adapter = pattern.sequence
    length = len(adapter)
    errors = pattern.max_errors

    if errors == 0:
        exact: list[tuple[int, int, int]] = []
        position = text.find(adapter)
        while position != -1:
            exact.append((position, position + length, 0))
            position = text.find(adapter, position + 1)
        return exact

    placements: dict[tuple[int, int], int] = {}
    checked: set[tuple[int, int]] = set()
    for seed_offset, seed in pattern.seeds:
        position = text.find(seed)
        while position != -1:
            anchor = position - seed_offset
            region = (max(0, anchor - errors), min(len(text), anchor + length + errors))
            if region not in checked:
                checked.add(region)
                hit = _align_adapter(adapter, text[region[0] : region[1]], errors)
                if hit is not None:
                    placements[(region[0] + hit[0], region[0] + hit[1])] = hit[2]
            position = text.find(seed, position + 1)

    return [(start, end, edits) for (start, end), edits in placements.items()]


def _best_placement(text: str, pattern: _AdapterPattern) -> tuple[int, int] | None:
    """Return (start, end) of the placement with the fewest edits (leftmost on ties), or None."""
    placements = _find_adapter_placements(text, pattern)
    if not placements:
        return None
    start, end, _ = min(placements, key=lambda placement: (placement[2], placement[0], placement[1]))
    return start, end


def _detect_with_plan(sequence: str, plan: list[_AdapterPattern], max_search_length: int) -> list[AdapterMatch]:
    """Find adapters anywhere in the upper-case ``sequence`` for every planned pattern."""
    matches: list[AdapterMatch] = []
    for pattern in plan:
        length = len(pattern.sequence)
        for start, end, edits in _find_adapter_placements(sequence, pattern):
            if start < max_search_length:
                location = "start"
            elif start > len(sequence) - max_search_length:
                location = "end"
            else:
                location = "internal"

            matches.append(
                AdapterMatch(
                    adapter_name=pattern.name,
                    adapter_sequence=pattern.sequence,
                    position=start,
                    end_position=end,
                    score=float(length - edits),
                    identity=(length - edits) / length,
                    location=location,
                    edit_distance=edits,
                    false_positive_bound=pattern.false_positive_bound,
                )
            )

    # Sort by position and deduplicate overlapping matches
    matches.sort(key=lambda m: m.position)
    return _deduplicate_adapter_matches(matches)


def _deduplicate_adapter_matches(matches: list[AdapterMatch]) -> list[AdapterMatch]:
    """Remove overlapping adapter matches, keeping the highest-scoring one."""
    if not matches:
        return []

    # Sort by score descending; ties go to the leftmost, shortest match so boundaries do not
    # depend on which seed proposed the placement
    sorted_matches = sorted(matches, key=lambda m: (-m.score, m.position, m.end_position))
    kept: list[AdapterMatch] = []

    for match in sorted_matches:
        overlap_found = False
        for kept_match in kept:
            # Check for overlap
            if match.position < kept_match.end_position and match.end_position > kept_match.position:
                overlap_found = True
                break
        if not overlap_found:
            kept.append(match)

    # Re-sort by position
    kept.sort(key=lambda m: m.position)
    return kept
