"""Consensus sequence generation and analysis utilities.

This module provides tools for generating consensus sequences from multiple
DNA alignments, handling ambiguity codes, and quality-based consensus calling.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)


def generate_consensus(sequences: List[str], threshold: float = 0.5) -> str:
    """Generate consensus sequence from multiple aligned sequences.

    Args:
        sequences: List of aligned DNA sequences
        threshold: Minimum frequency for consensus base (default: 0.5)

    Returns:
        Consensus DNA sequence

    Raises:
        ValueError: If sequences have different lengths or are empty

    Example:
        >>> seqs = ["ATCG", "ATCG", "ATCG"]
        >>> consensus = generate_consensus(seqs)
        >>> consensus == "ATCG"
        True
    """
    if not sequences:
        raise ValueError("Must provide at least one sequence")

    # Check all sequences have same length
    seq_lengths = [len(seq) for seq in sequences]
    if len(set(seq_lengths)) > 1:
        raise ValueError("All sequences must have the same length")

    if seq_lengths[0] == 0:
        return ""

    consensus = []

    for pos in range(seq_lengths[0]):
        # Count nucleotides at this position
        base_counts = {}
        total_bases = 0

        for seq in sequences:
            base = seq[pos].upper()
            if base in 'ATCG':
                base_counts[base] = base_counts.get(base, 0) + 1
                total_bases += 1

        if not base_counts:
            consensus.append('N')  # No valid bases
            continue

        # Find base with highest frequency
        max_base, max_count = max(base_counts.items(), key=lambda x: x[1])

        # Check if it meets threshold
        frequency = max_count / total_bases
        if frequency >= threshold:
            consensus.append(max_base)
        else:
            # Use ambiguity code
            consensus.append(_get_iupac_code(base_counts))

    return ''.join(consensus)


def consensus_with_ambiguity(sequences: List[str]) -> str:
    """Generate consensus sequence with IUPAC ambiguity codes.

    Args:
        sequences: List of aligned DNA sequences

    Returns:
        Consensus sequence with ambiguity codes

    Example:
        >>> seqs = ["ATCG", "ATCG", "AGCG"]
        >>> consensus = consensus_with_ambiguity(seqs)
        >>> consensus == "ATCG"  # Position 2: A and G -> R
        False
    """
    if not sequences:
        raise ValueError("Must provide at least one sequence")

    # Check all sequences have same length
    seq_lengths = [len(seq) for seq in sequences]
    if len(set(seq_lengths)) > 1:
        raise ValueError("All sequences must have the same length")

    if seq_lengths[0] == 0:
        return ""

    consensus = []

    for pos in range(seq_lengths[0]):
        # Count nucleotides at this position
        base_counts = {}

        for seq in sequences:
            base = seq[pos].upper()
            if base in 'ATCG':
                base_counts[base] = base_counts.get(base, 0) + 1

        if not base_counts:
            consensus.append('N')
        else:
            consensus.append(_get_iupac_code(base_counts))

    return ''.join(consensus)


def _get_iupac_code(base_counts: Dict[str, int]) -> str:
    """Get IUPAC ambiguity code for base composition.

    Args:
        base_counts: Dictionary of base counts

    Returns:
        IUPAC ambiguity code
    """
    bases = set(base_counts.keys())

    # Remove gaps and invalid characters
    bases.discard('-')
    bases = bases.intersection(set('ATCG'))

    if not bases:
        return 'N'
    elif len(bases) == 1:
        return list(bases)[0]
    elif bases == {'A', 'G'}:
        return 'R'
    elif bases == {'C', 'T'}:
        return 'Y'
    elif bases == {'G', 'C'}:
        return 'S'
    elif bases == {'A', 'T'}:
        return 'W'
    elif bases == {'G', 'T'}:
        return 'K'
    elif bases == {'A', 'C'}:
        return 'M'
    elif bases == {'A', 'C', 'G'}:
        return 'V'
    elif bases == {'A', 'C', 'T'}:
        return 'H'
    elif bases == {'A', 'G', 'T'}:
        return 'D'
    elif bases == {'C', 'G', 'T'}:
        return 'B'
    else:  # All four bases
        return 'N'


def quality_weighted_consensus(sequences: List[str], qualities: List[List[int]]) -> str:
    """Generate quality-weighted consensus sequence.

    Args:
        sequences: List of aligned DNA sequences
        qualities: List of quality score lists (one per sequence)

    Returns:
        Quality-weighted consensus sequence

    Raises:
        ValueError: If dimensions don't match

    Example:
        >>> seqs = ["ATCG", "ATCG"]
        >>> quals = [[30, 30, 30, 30], [30, 30, 30, 30]]
        >>> consensus = quality_weighted_consensus(seqs, quals)
        >>> consensus == "ATCG"
        True
    """
    if len(sequences) != len(qualities):
        raise ValueError("Number of sequences and quality lists must match")

    if not sequences:
        return ""

    seq_length = len(sequences[0])

    # Check dimensions
    for i, (seq, qual) in enumerate(zip(sequences, qualities)):
        if len(seq) != seq_length or len(qual) != seq_length:
            raise ValueError(f"Sequence {i} and quality dimensions don't match")

    consensus = []

    for pos in range(seq_length):
        # Calculate weighted base frequencies
        weighted_counts = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}

        for seq, qual_list in zip(sequences, qualities):
            base = seq[pos].upper()
            quality = qual_list[pos]

            if base in 'ATCG':
                # Use quality as weight (Phred score)
                weight = 10 ** (quality / -10.0)  # Convert to probability
                weighted_counts[base] += weight

        # Find base with highest weighted count
        if all(count == 0 for count in weighted_counts.values()):
            consensus.append('N')
        else:
            best_base = max(weighted_counts.items(), key=lambda x: x[1])[0]
            consensus.append(best_base)

    return ''.join(consensus)


def consensus_statistics(sequences: List[str]) -> Dict[str, any]:
    """Calculate statistics about sequence consensus.

    Args:
        sequences: List of aligned sequences

    Returns:
        Dictionary with consensus statistics

    Example:
        >>> seqs = ["ATCG", "ATCG", "ATCG"]
        >>> stats = consensus_statistics(seqs)
        >>> stats['conservation'] == 1.0
        True
    """
    if not sequences:
        return {
            'total_positions': 0,
            'conserved_positions': 0,
            'conservation': 0.0,
            'ambiguous_positions': 0,
            'gap_positions': 0
        }

    seq_length = len(sequences[0])
    conserved_positions = 0
    ambiguous_positions = 0
    gap_positions = 0

    for pos in range(seq_length):
        bases = []

        for seq in sequences:
            if pos < len(seq):
                base = seq[pos].upper()
                if base in 'ATCG-':
                    bases.append(base)

        if not bases:
            continue

        # Count unique bases (excluding gaps)
        unique_bases = set(b for b in bases if b != '-')
        gap_count = bases.count('-')

        if gap_count > 0:
            gap_positions += 1

        if len(unique_bases) == 1:
            conserved_positions += 1
        elif len(unique_bases) > 1:
            ambiguous_positions += 1

    conservation = conserved_positions / seq_length if seq_length > 0 else 0.0

    return {
        'total_positions': seq_length,
        'conserved_positions': conserved_positions,
        'conservation': conservation,
        'ambiguous_positions': ambiguous_positions,
        'gap_positions': gap_positions,
        'sequences_count': len(sequences)
    }


def majority_consensus(sequences: List[str]) -> str:
    """Generate majority rule consensus (simple majority vote).

    Args:
        sequences: List of aligned sequences

    Returns:
        Majority consensus sequence

    Example:
        >>> seqs = ["ATCG", "ATCG", "AGCG"]
        >>> consensus = majority_consensus(seqs)
        >>> consensus == "ATCG"  # Position 2: A wins majority
        True
    """
    return generate_consensus(sequences, threshold=0.5)


def strict_consensus(sequences: List[str]) -> str:
    """Generate strict consensus (100% agreement required).

    Args:
        sequences: List of aligned sequences

    Returns:
        Strict consensus sequence (uses ambiguity codes for disagreements)

    Example:
        >>> seqs = ["ATCG", "ATCG", "AGCG"]
        >>> consensus = strict_consensus(seqs)
        >>> consensus == "ATCG"  # Position 2: A and G -> R
        False
    """
    return consensus_with_ambiguity(sequences)


def consensus_from_alignment(alignment: List[str]) -> Tuple[str, Dict[str, any]]:
    """Generate consensus from sequence alignment with comprehensive analysis.

    Args:
        alignment: List of aligned sequences

    Returns:
        Tuple of (consensus_sequence, analysis_dict)

    Example:
        >>> alignment = ["ATCG-", "ATCG-", "ATCG-"]
        >>> consensus, analysis = consensus_from_alignment(alignment)
        >>> "statistics" in analysis
        True
    """
    if not alignment:
        return "", {}

    # Generate different types of consensus
    majority_cons = majority_consensus(alignment)
    strict_cons = strict_consensus(alignment)

    # Calculate statistics
    stats = consensus_statistics(alignment)

    # Choose best consensus (prefer majority if highly conserved)
    if stats['conservation'] > 0.8:
        final_consensus = majority_cons
    else:
        final_consensus = strict_cons

    analysis = {
        'majority_consensus': majority_cons,
        'strict_consensus': strict_cons,
        'final_consensus': final_consensus,
        'statistics': stats,
        'alignment_length': len(alignment[0]) if alignment else 0,
        'sequence_count': len(alignment)
    }

    return final_consensus, analysis


def find_consensus_breaks(sequences: List[str], window_size: int = 10) -> List[Tuple[int, float]]:
    """Find regions where consensus breaks down.

    Args:
        sequences: List of aligned sequences
        window_size: Size of sliding window for analysis

    Returns:
        List of tuples (position, conservation_score)

    Example:
        >>> seqs = ["ATCGATCG", "ATCGATCG", "GCTAGCTA"]
        >>> breaks = find_consensus_breaks(seqs, window_size=4)
        >>> isinstance(breaks, list)
        True
    """
    if not sequences or window_size <= 0:
        return []

    breaks = []
    seq_length = len(sequences[0])

    for start in range(0, seq_length - window_size + 1, window_size // 2):
        end = min(start + window_size, seq_length)

        # Extract window from all sequences
        window_seqs = []
        for seq in sequences:
            window_seqs.append(seq[start:end])

        # Calculate conservation for this window
        stats = consensus_statistics(window_seqs)
        conservation = stats['conservation']

        breaks.append((start, conservation))

    return breaks


def bootstrap_consensus(sequences: List[str], n_bootstraps: int = 100) -> Tuple[str, float]:
    """Generate consensus with bootstrap confidence.

    Args:
        sequences: List of aligned sequences
        n_bootstraps: Number of bootstrap replicates

    Returns:
        Tuple of (consensus_sequence, confidence_score)

    Example:
        >>> seqs = ["ATCG", "ATCG", "ATCG"]
        >>> consensus, confidence = bootstrap_consensus(seqs, n_bootstraps=10)
        >>> confidence > 0.8
        True
    """
    if not sequences or n_bootstraps <= 0:
        return "", 0.0

    import random

    consensus_sequences = []

    # Generate bootstrap replicates
    for _ in range(n_bootstraps):
        # Sample with replacement
        bootstrap_sample = [random.choice(sequences) for _ in range(len(sequences))]
        consensus = majority_consensus(bootstrap_sample)
        consensus_sequences.append(consensus)

    # Calculate consensus from all bootstrap consensuses
    final_consensus = majority_consensus(consensus_sequences)

    # Calculate confidence (fraction of bootstraps agreeing with final consensus)
    agreements = sum(1 for cons in consensus_sequences if cons == final_consensus)
    confidence = agreements / n_bootstraps

    return final_consensus, confidence


