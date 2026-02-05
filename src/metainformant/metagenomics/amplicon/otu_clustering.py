"""OTU clustering for amplicon metagenomics.

Implements greedy centroid-based OTU clustering (VSEARCH/UCLUST-style),
pairwise sequence identity calculation, and de novo chimera detection
using the UCHIME algorithm approach.

The clustering algorithm sorts sequences by abundance (or length), then
iteratively assigns each sequence to the most similar existing centroid
if the identity exceeds the threshold, or creates a new OTU centroid.
"""

from __future__ import annotations

import os
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Config prefix for environment overrides
_ENV_PREFIX = "META_"


@dataclass
class OTU:
    """Represents an Operational Taxonomic Unit.

    Attributes:
        centroid_id: Identifier of the centroid sequence.
        centroid_sequence: The representative centroid sequence.
        member_ids: List of sequence identifiers assigned to this OTU.
        member_sequences: List of member sequences.
        size: Number of sequences in this OTU.
    """

    centroid_id: str
    centroid_sequence: str
    member_ids: list[str] = field(default_factory=list)
    member_sequences: list[str] = field(default_factory=list)
    size: int = 1

    def add_member(self, seq_id: str, sequence: str) -> None:
        """Add a member sequence to this OTU."""
        self.member_ids.append(seq_id)
        self.member_sequences.append(sequence)
        self.size += 1


@dataclass
class ClusteringResult:
    """Result of OTU clustering.

    Attributes:
        otus: List of OTU objects.
        threshold: Identity threshold used for clustering.
        total_sequences: Total number of input sequences.
        num_otus: Number of OTUs formed.
        otu_table: Mapping of OTU centroid IDs to member counts.
    """

    otus: list[OTU]
    threshold: float
    total_sequences: int
    num_otus: int
    otu_table: dict[str, int] = field(default_factory=dict)


def calculate_identity(seq1: str, seq2: str) -> float:
    """Calculate pairwise sequence identity using global alignment.

    Uses a banded Needleman-Wunsch global alignment with affine gap penalties
    to compute the fraction of identical positions in the alignment.

    Args:
        seq1: First nucleotide sequence (uppercase ACGT expected).
        seq2: Second nucleotide sequence (uppercase ACGT expected).

    Returns:
        Fraction of identical positions (0.0 to 1.0).

    Raises:
        ValueError: If either sequence is empty.

    Examples:
        >>> calculate_identity("ATCGATCG", "ATCGATCG")
        1.0
        >>> calculate_identity("ATCGATCG", "ATCGTTCG")
        0.875
    """
    if not seq1 or not seq2:
        raise ValueError("Both sequences must be non-empty")

    s1 = seq1.upper().replace("-", "")
    s2 = seq2.upper().replace("-", "")

    if s1 == s2:
        return 1.0

    # Needleman-Wunsch global alignment with affine gap penalties
    match_score = 2
    mismatch_penalty = -1
    gap_open = -5
    gap_extend = -1

    n, m = len(s1), len(s2)

    # Three matrices: M (match/mismatch), X (gap in seq2), Y (gap in seq1)
    NEG_INF = float("-inf")

    # For memory efficiency, use two-row approach
    prev_m = [NEG_INF] * (m + 1)
    prev_x = [NEG_INF] * (m + 1)
    prev_y = [NEG_INF] * (m + 1)

    curr_m = [NEG_INF] * (m + 1)
    curr_x = [NEG_INF] * (m + 1)
    curr_y = [NEG_INF] * (m + 1)

    # Traceback matrices (we need full matrices for traceback)
    trace_m = [[0] * (m + 1) for _ in range(n + 1)]
    trace_x = [[0] * (m + 1) for _ in range(n + 1)]
    trace_y = [[0] * (m + 1) for _ in range(n + 1)]

    # Store all rows for traceback
    all_m = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    all_x = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    all_y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    # Initialize
    all_m[0][0] = 0
    for j in range(1, m + 1):
        all_y[0][j] = gap_open + gap_extend * j
        all_m[0][j] = NEG_INF
        all_x[0][j] = NEG_INF
    for i in range(1, n + 1):
        all_x[i][0] = gap_open + gap_extend * i
        all_m[i][0] = NEG_INF
        all_y[i][0] = NEG_INF

    # Fill
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Score for match/mismatch
            if s1[i - 1] == s2[j - 1]:
                diag_score = match_score
            else:
                diag_score = mismatch_penalty

            # M matrix: came from a match/mismatch
            m_from_m = all_m[i - 1][j - 1] + diag_score
            m_from_x = all_x[i - 1][j - 1] + diag_score
            m_from_y = all_y[i - 1][j - 1] + diag_score
            all_m[i][j] = max(m_from_m, m_from_x, m_from_y)

            # X matrix: gap in seq2 (insertion in seq1)
            x_from_m = all_m[i - 1][j] + gap_open + gap_extend
            x_from_x = all_x[i - 1][j] + gap_extend
            all_x[i][j] = max(x_from_m, x_from_x)

            # Y matrix: gap in seq1 (insertion in seq2)
            y_from_m = all_m[i][j - 1] + gap_open + gap_extend
            y_from_y = all_y[i][j - 1] + gap_extend
            all_y[i][j] = max(y_from_m, y_from_y)

    # Best terminal score
    final_score = max(all_m[n][m], all_x[n][m], all_y[n][m])

    # Traceback to count matches
    i, j = n, m
    best_end = max(all_m[n][m], all_x[n][m], all_y[n][m])
    if best_end == all_m[n][m]:
        state = "M"
    elif best_end == all_x[n][m]:
        state = "X"
    else:
        state = "Y"

    matches = 0
    alignment_length = 0

    while i > 0 or j > 0:
        if state == "M":
            if i <= 0 or j <= 0:
                break
            if s1[i - 1] == s2[j - 1]:
                matches += 1
            alignment_length += 1
            # Determine which matrix we came from
            diag_score = match_score if s1[i - 1] == s2[j - 1] else mismatch_penalty
            prev_scores = {
                "M": all_m[i - 1][j - 1] + diag_score,
                "X": all_x[i - 1][j - 1] + diag_score,
                "Y": all_y[i - 1][j - 1] + diag_score,
            }
            state = max(prev_scores, key=prev_scores.get)  # type: ignore[arg-type]
            i -= 1
            j -= 1
        elif state == "X":
            if i <= 0:
                break
            alignment_length += 1
            x_from_m = all_m[i - 1][j] + gap_open + gap_extend
            x_from_x = all_x[i - 1][j] + gap_extend
            if x_from_m >= x_from_x:
                state = "M"
            else:
                state = "X"
            i -= 1
        elif state == "Y":
            if j <= 0:
                break
            alignment_length += 1
            y_from_m = all_m[i][j - 1] + gap_open + gap_extend
            y_from_y = all_y[i][j - 1] + gap_extend
            if y_from_m >= y_from_y:
                state = "M"
            else:
                state = "Y"
            j -= 1

    if alignment_length == 0:
        return 0.0

    return matches / alignment_length


def _quick_identity(seq1: str, seq2: str) -> float:
    """Fast approximate identity using k-mer overlap (pre-filter).

    Uses shared 8-mer fraction as a quick upper bound estimate to avoid
    expensive full alignments for clearly dissimilar sequences.

    Args:
        seq1: First sequence.
        seq2: Second sequence.

    Returns:
        Approximate identity estimate (0.0 to 1.0).
    """
    k = 8
    if len(seq1) < k or len(seq2) < k:
        # Fall back to simple comparison for short sequences
        min_len = min(len(seq1), len(seq2))
        if min_len == 0:
            return 0.0
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / max(len(seq1), len(seq2))

    kmers1 = set()
    for i in range(len(seq1) - k + 1):
        kmers1.add(seq1[i : i + k])

    kmers2 = set()
    for i in range(len(seq2) - k + 1):
        kmers2.add(seq2[i : i + k])

    if not kmers1 or not kmers2:
        return 0.0

    shared = len(kmers1 & kmers2)
    total = max(len(kmers1), len(kmers2))
    return shared / total


def cluster_otus(
    sequences: dict[str, str],
    threshold: float = 0.97,
    abundance: dict[str, int] | None = None,
    sort_by: str = "length",
    prefilter: bool = True,
) -> ClusteringResult:
    """Cluster sequences into OTUs using greedy centroid-based algorithm.

    Implements a VSEARCH/UCLUST-style greedy clustering:
    1. Sort sequences by abundance (descending) or length (descending).
    2. The first sequence becomes the centroid of OTU #1.
    3. For each subsequent sequence, compare to all existing centroids.
    4. If identity >= threshold for any centroid, assign to that OTU.
    5. Otherwise, create a new OTU with this sequence as centroid.

    Args:
        sequences: Dictionary mapping sequence IDs to nucleotide sequences.
        threshold: Minimum sequence identity to cluster (0.0 to 1.0, default 0.97).
        abundance: Optional dictionary of sequence abundances for sorting.
            If None and sort_by="abundance", sequences are sorted by length.
        sort_by: Sorting strategy - "length" (default) or "abundance".
        prefilter: If True, use k-mer pre-filter for speed (default True).

    Returns:
        ClusteringResult with OTU assignments and statistics.

    Raises:
        ValueError: If threshold is not in (0.0, 1.0] or sequences is empty.

    Examples:
        >>> seqs = {"s1": "ATCGATCGATCG", "s2": "ATCGATCGATCG", "s3": "TTTTAAAA"}
        >>> result = cluster_otus(seqs, threshold=0.97)
        >>> result.num_otus
        2
    """
    if not sequences:
        raise ValueError("Input sequences dictionary must not be empty")
    if not 0.0 < threshold <= 1.0:
        raise ValueError(f"Threshold must be in (0.0, 1.0], got {threshold}")

    # Environment override for default threshold
    env_threshold = os.environ.get(f"{_ENV_PREFIX}OTU_THRESHOLD")
    if env_threshold is not None:
        try:
            threshold = float(env_threshold)
            logger.info(f"Using threshold from environment: {threshold}")
        except ValueError:
            logger.warning(f"Invalid {_ENV_PREFIX}OTU_THRESHOLD value: {env_threshold}, using {threshold}")

    logger.info(f"Clustering {len(sequences)} sequences at {threshold:.2f} identity threshold")

    # Sort sequences
    if sort_by == "abundance" and abundance:
        sorted_ids = sorted(sequences.keys(), key=lambda x: (-abundance.get(x, 0), -len(sequences[x])))
    else:
        sorted_ids = sorted(sequences.keys(), key=lambda x: -len(sequences[x]))

    otus: list[OTU] = []
    prefilter_threshold = threshold - 0.10  # k-mer filter is a loose upper bound

    for seq_id in sorted_ids:
        seq = sequences[seq_id].upper().replace("-", "").replace(".", "")
        if not seq:
            logger.warning(f"Skipping empty sequence: {seq_id}")
            continue

        assigned = False
        best_identity = 0.0
        best_otu_idx = -1

        for idx, otu in enumerate(otus):
            centroid_seq = otu.centroid_sequence

            # Pre-filter: skip expensive alignment if k-mer overlap is too low
            if prefilter and len(seq) >= 8 and len(centroid_seq) >= 8:
                quick_id = _quick_identity(seq, centroid_seq)
                if quick_id < prefilter_threshold:
                    continue

            identity = calculate_identity(seq, centroid_seq)
            if identity >= threshold and identity > best_identity:
                best_identity = identity
                best_otu_idx = idx

        if best_otu_idx >= 0:
            otus[best_otu_idx].add_member(seq_id, seq)
            assigned = True
        else:
            new_otu = OTU(
                centroid_id=seq_id,
                centroid_sequence=seq,
                member_ids=[seq_id],
                member_sequences=[seq],
                size=1,
            )
            otus.append(new_otu)

    # Build OTU table
    otu_table = {otu.centroid_id: otu.size for otu in otus}

    result = ClusteringResult(
        otus=otus,
        threshold=threshold,
        total_sequences=len(sequences),
        num_otus=len(otus),
        otu_table=otu_table,
    )

    logger.info(f"Clustering complete: {result.num_otus} OTUs from {result.total_sequences} sequences")
    return result


def filter_chimeras(
    sequences: dict[str, str],
    reference_db: dict[str, str] | None = None,
    abundance: dict[str, int] | None = None,
    min_divergence: float = 1.0,
    min_score: float = 0.28,
) -> dict[str, bool]:
    """Detect and flag chimeric sequences using UCHIME-style de novo detection.

    The UCHIME algorithm detects chimeras by:
    1. Sorting sequences by abundance (most abundant first).
    2. For each query, finding the two best parent candidates among more
       abundant sequences.
    3. Building an optimal chimeric model by combining segments of the two parents.
    4. Computing a score based on divergence between the query and the model
       versus divergence from the parents.
    5. Flagging sequences with chimera scores above the threshold.

    When reference_db is provided, reference sequences are used as potential
    parents instead of the abundance-sorted input sequences.

    Args:
        sequences: Dictionary mapping sequence IDs to nucleotide sequences.
        reference_db: Optional reference database of known non-chimeric sequences.
        abundance: Optional abundance counts for abundance-based sorting.
        min_divergence: Minimum divergence ratio between chimeric model and parents.
        min_score: Minimum chimera score threshold to flag as chimeric (0.0-1.0).

    Returns:
        Dictionary mapping sequence IDs to boolean (True = chimeric, False = non-chimeric).

    Examples:
        >>> seqs = {"s1": "AAAAAAAAAA", "s2": "TTTTTTTTTT", "s3": "AAAAATTTTT"}
        >>> chimeras = filter_chimeras(seqs)
        >>> isinstance(chimeras, dict)
        True
    """
    if not sequences:
        return {}

    logger.info(f"Checking {len(sequences)} sequences for chimeras")

    chimera_flags: dict[str, bool] = {}

    # Sort by abundance (most abundant are least likely chimeric)
    if abundance:
        sorted_ids = sorted(sequences.keys(), key=lambda x: -abundance.get(x, 0))
    else:
        sorted_ids = sorted(sequences.keys(), key=lambda x: -len(sequences[x]))

    # Use reference db or build parent pool from more abundant sequences
    parent_pool: dict[str, str] = {}
    if reference_db:
        parent_pool = {k: v.upper() for k, v in reference_db.items()}

    for idx, seq_id in enumerate(sorted_ids):
        query = sequences[seq_id].upper()

        # If no reference DB, parents are more abundant sequences
        if not reference_db:
            if idx < 2:
                # Most abundant sequences are assumed non-chimeric
                chimera_flags[seq_id] = False
                parent_pool[seq_id] = query
                continue

        if not parent_pool:
            chimera_flags[seq_id] = False
            parent_pool[seq_id] = query
            continue

        # Find two best matching parents
        parent_scores: list[tuple[str, float]] = []
        for pid, pseq in parent_pool.items():
            if pid == seq_id:
                continue
            # Use quick identity for speed
            id_score = _quick_identity(query, pseq)
            parent_scores.append((pid, id_score))

        parent_scores.sort(key=lambda x: -x[1])

        if len(parent_scores) < 2:
            chimera_flags[seq_id] = False
            if not reference_db:
                parent_pool[seq_id] = query
            continue

        parent_a_id, score_a = parent_scores[0]
        parent_b_id, score_b = parent_scores[1]
        parent_a = parent_pool[parent_a_id]
        parent_b = parent_pool[parent_b_id]

        # Build chimeric model: try all possible breakpoints
        # A chimera is formed by joining a left segment from parent A
        # with a right segment from parent B (or vice versa)
        best_chimera_score = 0.0
        query_len = len(query)

        if query_len < 10:
            chimera_flags[seq_id] = False
            if not reference_db:
                parent_pool[seq_id] = query
            continue

        # Evaluate breakpoints at 10% intervals for efficiency
        step = max(1, query_len // 10)
        for bp in range(step, query_len - step + 1, step):
            # Model: left from A, right from B
            left_a = parent_a[:bp] if bp <= len(parent_a) else parent_a
            right_b = parent_b[bp:] if bp < len(parent_b) else ""
            model_ab = left_a + right_b

            # Model: left from B, right from A
            left_b = parent_b[:bp] if bp <= len(parent_b) else parent_b
            right_a = parent_a[bp:] if bp < len(parent_a) else ""
            model_ba = left_b + right_a

            # Score: how well does the query match the chimeric model
            # vs how well it matches each parent individually
            for model in [model_ab, model_ba]:
                if not model:
                    continue
                # Calculate identity between query and chimeric model
                min_len = min(len(query), len(model))
                if min_len == 0:
                    continue
                matches = sum(1 for a, b in zip(query[:min_len], model[:min_len]) if a == b)
                model_identity = matches / min_len

                # Calculate identity to each parent directly
                min_len_a = min(len(query), len(parent_a))
                matches_a = sum(1 for a, b in zip(query[:min_len_a], parent_a[:min_len_a]) if a == b)
                parent_a_identity = matches_a / min_len_a if min_len_a > 0 else 0.0

                min_len_b = min(len(query), len(parent_b))
                matches_b = sum(1 for a, b in zip(query[:min_len_b], parent_b[:min_len_b]) if a == b)
                parent_b_identity = matches_b / min_len_b if min_len_b > 0 else 0.0

                max_parent_identity = max(parent_a_identity, parent_b_identity)

                # Chimera score: model fits better than either parent alone
                if max_parent_identity > 0:
                    divergence = model_identity - max_parent_identity
                    chimera_score = divergence / max_parent_identity if max_parent_identity > 0 else 0.0
                else:
                    chimera_score = 0.0

                best_chimera_score = max(best_chimera_score, chimera_score)

        is_chimeric = best_chimera_score >= min_score
        chimera_flags[seq_id] = is_chimeric

        if is_chimeric:
            logger.debug(f"Chimera detected: {seq_id} (score={best_chimera_score:.3f})")

        # Add non-chimeric sequences to parent pool for de novo mode
        if not reference_db and not is_chimeric:
            parent_pool[seq_id] = query

    chimeric_count = sum(1 for v in chimera_flags.values() if v)
    logger.info(f"Chimera detection complete: {chimeric_count}/{len(chimera_flags)} flagged as chimeric")
    return chimera_flags


__all__ = [
    "OTU",
    "ClusteringResult",
    "calculate_identity",
    "cluster_otus",
    "filter_chimeras",
]
