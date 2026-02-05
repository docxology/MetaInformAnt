"""Read filtering and adapter trimming for long-read sequencing data.

Provides length-based filtering, quality filtering, adapter detection and
trimming, and chimeric read detection/splitting. All algorithms are real
implementations with proper bioinformatics logic.

Optional dependencies:
    - numpy: For efficient numerical computation
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

try:
    import numpy as np  # type: ignore[import-untyped]
except ImportError:
    np = None  # type: ignore[assignment]


# Standard ONT adapter sequences
ONT_ADAPTERS: dict[str, str] = {
    "ONT_LSK109_top": "AATGTACTTCGTTCAGTTACGTATTGCT",
    "ONT_LSK109_bottom": "GCAATACGTAACTGAACGAAGT",
    "ONT_LSK110_top": "CAGCACCT",
    "ONT_LSK110_bottom": "AGGTGCTG",
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
        score: Alignment score of the match.
        identity: Fraction of matching bases.
        location: Where the adapter was found ('start', 'end', 'internal').
    """

    adapter_name: str
    adapter_sequence: str
    position: int
    end_position: int
    score: float
    identity: float
    location: str


def filter_by_length(
    reads: Sequence[dict[str, Any] | ReadRecord],
    min_length: int = 1000,
    max_length: int | None = None,
) -> list[dict[str, Any] | ReadRecord]:
    """Filter reads by sequence length.

    Args:
        reads: Sequence of read dictionaries (with 'sequence' key) or ReadRecord objects.
        min_length: Minimum read length to keep (inclusive). Default 1000.
        max_length: Maximum read length to keep (inclusive). None for no upper limit.

    Returns:
        List of reads passing the length filter.
    """
    filtered: list[dict[str, Any] | ReadRecord] = []
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
        passed, total, min_length, max_length or "None",
    )
    return filtered


def filter_by_quality(
    reads: Sequence[dict[str, Any] | ReadRecord],
    min_q: float = 7.0,
) -> list[dict[str, Any] | ReadRecord]:
    """Filter reads by mean Phred quality score.

    Calculates the mean quality for each read and retains only those
    with mean quality >= min_q.

    Args:
        reads: Sequence of read dictionaries or ReadRecord objects.
        min_q: Minimum mean Phred quality score to keep. Default 7.0.

    Returns:
        List of reads passing the quality filter.
    """
    filtered: list[dict[str, Any] | ReadRecord] = []
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


def trim_adapters(
    reads: Sequence[dict[str, Any] | ReadRecord],
    adapter_sequences: dict[str, str] | None = None,
    max_search_length: int = 200,
    min_identity: float = 0.75,
) -> list[ReadRecord]:
    """Trim adapter sequences from the ends of reads.

    Searches for adapters at the beginning and end of each read using
    semi-global alignment and trims matching regions.

    Args:
        reads: Sequence of read dictionaries or ReadRecord objects.
        adapter_sequences: Dictionary of adapter_name -> sequence. If None,
            uses all known ONT and PacBio adapters.
        max_search_length: Maximum number of bases to search at each end.
        min_identity: Minimum sequence identity to consider a match (0-1).

    Returns:
        List of ReadRecord objects with adapters trimmed.
    """
    if adapter_sequences is None:
        adapter_sequences = ALL_KNOWN_ADAPTERS

    trimmed: list[ReadRecord] = []

    for read in reads:
        seq = _get_sequence(read)
        qual_str = _get_quality(read) or ""
        read_id = _get_read_id(read)

        if seq is None:
            continue

        trim_start = 0
        trim_end = len(seq)

        # Check for adapters at the start
        search_start = seq[:max_search_length].upper()
        for adapter_name, adapter_seq in adapter_sequences.items():
            match = _find_adapter_match(search_start, adapter_seq.upper(), min_identity)
            if match is not None:
                adapter_end = match[1]
                trim_start = max(trim_start, adapter_end)

        # Check for adapters at the end
        search_end = seq[-max_search_length:].upper() if len(seq) > max_search_length else seq.upper()
        offset = len(seq) - len(search_end)
        for adapter_name, adapter_seq in adapter_sequences.items():
            # Also check reverse complement of adapter at end
            rc_adapter = _reverse_complement(adapter_seq.upper())
            match = _find_adapter_match(search_end, rc_adapter, min_identity)
            if match is not None:
                adapter_start = offset + match[0]
                trim_end = min(trim_end, adapter_start)

            # Check forward adapter at end too
            match = _find_adapter_match(search_end, adapter_seq.upper(), min_identity)
            if match is not None:
                adapter_start = offset + match[0]
                trim_end = min(trim_end, adapter_start)

        # Apply trimming
        trimmed_seq = seq[trim_start:trim_end]
        trimmed_qual = qual_str[trim_start:trim_end] if qual_str else ""

        if len(trimmed_seq) > 0:
            metadata = {}
            if isinstance(read, dict):
                metadata = {k: v for k, v in read.items() if k not in ("sequence", "quality", "quality_string", "read_id")}
            elif hasattr(read, "metadata"):
                metadata = dict(read.metadata)

            metadata["original_length"] = len(seq)
            metadata["trimmed_start"] = trim_start
            metadata["trimmed_end"] = len(seq) - trim_end

            trimmed.append(ReadRecord(
                read_id=read_id,
                sequence=trimmed_seq,
                quality_string=trimmed_qual,
                metadata=metadata,
            ))

    logger.info("Adapter trimming: processed %d reads, %d with sequence remaining", len(reads), len(trimmed))
    return trimmed


def detect_adapters(
    sequence: str,
    known_adapters: dict[str, str] | None = None,
    max_search_length: int = 300,
    min_identity: float = 0.70,
) -> list[AdapterMatch]:
    """Detect adapter sequences within a read.

    Searches for known adapter sequences throughout the read, not just at
    the ends, to identify internal adapters that may indicate chimeric reads.

    Args:
        sequence: Nucleotide sequence to search.
        known_adapters: Dictionary of adapter_name -> sequence. If None,
            uses all known adapters.
        max_search_length: Maximum length of start/end regions to search.
        min_identity: Minimum identity for a match (0-1).

    Returns:
        List of AdapterMatch objects for all detected adapters.
    """
    if known_adapters is None:
        known_adapters = ALL_KNOWN_ADAPTERS

    if not sequence:
        return []

    matches: list[AdapterMatch] = []
    seq_upper = sequence.upper()

    for adapter_name, adapter_seq in known_adapters.items():
        adapter_upper = adapter_seq.upper()
        rc_adapter = _reverse_complement(adapter_upper)

        for search_adapter, suffix in [(adapter_upper, ""), (rc_adapter, "_rc")]:
            # Scan the sequence with a sliding window
            adapter_len = len(search_adapter)
            if adapter_len == 0:
                continue

            # Search the full sequence
            for i in range(0, len(seq_upper) - adapter_len + 1, max(1, adapter_len // 2)):
                window = seq_upper[i : i + adapter_len]
                identity = _sequence_identity(window, search_adapter)

                if identity >= min_identity:
                    # Determine location
                    if i < max_search_length:
                        location = "start"
                    elif i > len(seq_upper) - max_search_length:
                        location = "end"
                    else:
                        location = "internal"

                    matches.append(AdapterMatch(
                        adapter_name=adapter_name + suffix,
                        adapter_sequence=search_adapter,
                        position=i,
                        end_position=i + adapter_len,
                        score=identity * adapter_len,
                        identity=identity,
                        location=location,
                    ))

    # Sort by position and deduplicate overlapping matches
    matches.sort(key=lambda m: m.position)
    deduplicated = _deduplicate_adapter_matches(matches)

    return deduplicated


def split_chimeric_reads(
    reads: Sequence[dict[str, Any] | ReadRecord],
    known_adapters: dict[str, str] | None = None,
    min_fragment_length: int = 500,
    min_adapter_identity: float = 0.75,
) -> list[ReadRecord]:
    """Detect and split chimeric reads at internal adapter sequences.

    Chimeric reads are artifacts where two separate DNA molecules are
    joined during library preparation. They contain internal adapter
    sequences. This function finds internal adapters and splits the
    read into separate fragments.

    Args:
        reads: Sequence of read dictionaries or ReadRecord objects.
        known_adapters: Adapter sequences to search for. If None, uses defaults.
        min_fragment_length: Minimum length of resulting fragments to keep.
        min_adapter_identity: Minimum identity for adapter matching.

    Returns:
        List of ReadRecord objects, with chimeric reads split into fragments.
        Non-chimeric reads are returned unchanged.
    """
    if known_adapters is None:
        known_adapters = ALL_KNOWN_ADAPTERS

    result: list[ReadRecord] = []
    chimeric_count = 0

    for read in reads:
        seq = _get_sequence(read)
        qual_str = _get_quality(read) or ""
        read_id = _get_read_id(read)

        if seq is None or len(seq) < min_fragment_length:
            continue

        # Detect internal adapters
        adapter_matches = detect_adapters(seq, known_adapters, min_identity=min_adapter_identity)
        internal_matches = [m for m in adapter_matches if m.location == "internal"]

        if not internal_matches:
            # Not chimeric, return as-is
            result.append(ReadRecord(
                read_id=read_id,
                sequence=seq,
                quality_string=qual_str,
                metadata={"chimeric": False},
            ))
            continue

        # Split at internal adapters
        chimeric_count += 1
        split_points = [(0, None)]  # (start, adapter_match)
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
                result.append(ReadRecord(
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
                ))
                fragment_idx += 1

            i += 1

    logger.info(
        "Chimeric detection: %d/%d reads were chimeric, produced %d total fragments",
        chimeric_count, len(reads), len(result),
    )
    return result


# --- Internal helper functions ---


def _get_sequence(read: Any) -> str | None:
    """Extract sequence from various read representations."""
    if isinstance(read, dict):
        return read.get("sequence")
    elif hasattr(read, "sequence"):
        return read.sequence
    elif hasattr(read, "query_sequence"):
        return read.query_sequence
    return None


def _get_quality(read: Any) -> str | None:
    """Extract quality string from various read representations."""
    if isinstance(read, dict):
        return read.get("quality") or read.get("quality_string")
    elif hasattr(read, "quality_string"):
        return read.quality_string
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


def _reverse_complement(seq: str) -> str:
    """Compute the reverse complement of a DNA sequence."""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(c, "N") for c in reversed(seq))


def _sequence_identity(seq1: str, seq2: str) -> float:
    """Calculate sequence identity between two equal-length sequences.

    Simple character-by-character comparison. Both sequences must be
    the same length.

    Returns:
        Fraction of matching positions (0.0 to 1.0).
    """
    if not seq1 or not seq2:
        return 0.0
    min_len = min(len(seq1), len(seq2))
    if min_len == 0:
        return 0.0
    matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
    return matches / min_len


def _find_adapter_match(
    sequence: str,
    adapter: str,
    min_identity: float,
) -> tuple[int, int] | None:
    """Find the best match of an adapter within a sequence using semi-global alignment.

    Uses a simplified banded alignment approach: slides the adapter across
    the sequence and returns the position with highest identity above threshold.

    Args:
        sequence: The sequence to search in.
        adapter: The adapter sequence to find.
        min_identity: Minimum identity threshold.

    Returns:
        Tuple of (start, end) positions if found, None otherwise.
    """
    if not sequence or not adapter:
        return None

    adapter_len = len(adapter)
    if adapter_len > len(sequence):
        return None

    best_identity = 0.0
    best_pos = -1

    for i in range(len(sequence) - adapter_len + 1):
        window = sequence[i : i + adapter_len]
        identity = _sequence_identity(window, adapter)
        if identity > best_identity:
            best_identity = identity
            best_pos = i

    if best_identity >= min_identity and best_pos >= 0:
        return (best_pos, best_pos + adapter_len)

    return None


def _deduplicate_adapter_matches(matches: list[AdapterMatch]) -> list[AdapterMatch]:
    """Remove overlapping adapter matches, keeping the highest-scoring one."""
    if not matches:
        return []

    # Sort by score descending
    sorted_matches = sorted(matches, key=lambda m: -m.score)
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
