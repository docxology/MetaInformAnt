"""Motif discovery and analysis for DNA sequences.

This module provides tools for finding, analyzing, and scoring DNA sequence motifs,
including position weight matrices (PWMs), motif conservation analysis, and
pattern matching with IUPAC ambiguity codes.
"""

from __future__ import annotations

import re
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from metainformant.core import logging

logger = logging.get_logger(__name__)


def find_motifs(seq: str, motif_patterns: List[str]) -> Dict[str, List[int]]:
    """Find occurrences of multiple motifs in a DNA sequence.

    Supports IUPAC ambiguity codes for flexible pattern matching.

    Args:
        seq: DNA sequence to search
        motif_patterns: List of motif patterns (can include IUPAC codes)

    Returns:
        Dictionary mapping motif patterns to lists of start positions

    Example:
        >>> seq = "ATCGATCGATCG"
        >>> motifs = find_motifs(seq, ["ATC", "CG"])
        >>> len(motifs["ATC"])
        3
    """
    if not seq:
        return {}

    results = {}

    for motif in motif_patterns:
        positions = find_motif_positions(seq, motif)
        results[motif] = positions

    return results


def find_motif_positions(seq: str, motif: str) -> List[int]:
    """Find all occurrences of a motif in a DNA sequence.

    Supports IUPAC ambiguity codes:
    - R: A/G, Y: C/T, S: G/C, W: A/T, K: G/T, M: A/C
    - B: C/G/T, D: A/G/T, H: A/C/T, V: A/C/G, N: A/C/G/T

    Args:
        seq: DNA sequence to search
        motif: Motif pattern (may include IUPAC codes)

    Returns:
        List of start positions where motif occurs

    Example:
        >>> find_motif_positions("ATCGATCG", "ATC")
        [0, 4]
    """
    if not seq or not motif:
        return []

    # Convert IUPAC codes to regex patterns
    iupac_to_regex = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
    }

    # Build regex pattern
    regex_pattern = ''
    for char in motif.upper():
        if char in iupac_to_regex:
            regex_pattern += iupac_to_regex[char]
        else:
            # Treat unknown characters as literals (may not match)
            regex_pattern += re.escape(char)

    # Find all matches
    positions = []
    for match in re.finditer(regex_pattern, seq.upper()):
        positions.append(match.start())

    return positions


def create_pwm(sequences: List[str]) -> pd.DataFrame:
    """Create a position weight matrix (PWM) from aligned sequences.

    Args:
        sequences: List of aligned DNA sequences (must be same length)

    Returns:
        Position weight matrix as pandas DataFrame

    Raises:
        ValueError: If sequences have different lengths or are empty

    Example:
        >>> seqs = ["ATCG", "ATCG", "ATCG"]
        >>> pwm = create_pwm(seqs)
        >>> pwm.shape
        (4, 4)
    """
    if not sequences:
        raise ValueError("Must provide at least one sequence")

    # Check all sequences have same length
    seq_lengths = [len(seq) for seq in sequences]
    if len(set(seq_lengths)) > 1:
        raise ValueError("All sequences must have the same length")

    seq_length = seq_lengths[0]
    if seq_length == 0:
        raise ValueError("Sequences cannot be empty")

    # Initialize PWM matrix
    nucleotides = ['A', 'C', 'G', 'T']
    pwm = np.zeros((seq_length, 4))

    # Count nucleotide frequencies at each position
    for seq in sequences:
        for pos, nucleotide in enumerate(seq.upper()):
            if nucleotide in nucleotides:
                nuc_idx = nucleotides.index(nucleotide)
                pwm[pos, nuc_idx] += 1

    # Convert to frequencies, then log-odds (simple version)
    pwm = pwm / len(sequences)

    # Add pseudocounts to avoid log(0)
    pseudocount = 0.1
    pwm += pseudocount

    # Normalize
    pwm = pwm / pwm.sum(axis=1, keepdims=True)

    # Create DataFrame
    df = pd.DataFrame(
        pwm,
        columns=nucleotides,
        index=[f'pos_{i}' for i in range(seq_length)]
    )

    return df


def score_sequence_pwm(sequence: str, pwm: pd.DataFrame) -> List[float]:
    """Score a sequence using a position weight matrix.

    Args:
        sequence: DNA sequence to score
        pwm: Position weight matrix from create_pwm()

    Returns:
        List of scores for each position in the sequence

    Raises:
        ValueError: If PWM format is invalid

    Example:
        >>> seqs = ["ATCG", "ATCG", "ATCG"]
        >>> pwm = create_pwm(seqs)
        >>> scores = score_sequence_pwm("ATCGATCG", pwm)
        >>> len(scores)
        5
    """
    if pwm.empty:
        raise ValueError("PWM cannot be empty")

    if not sequence:
        return []

    seq_length = len(sequence)
    pwm_length = len(pwm)
    nucleotides = ['A', 'C', 'G', 'T']

    scores = []

    # Slide PWM along sequence
    for start_pos in range(seq_length - pwm_length + 1):
        score = 0.0

        for pwm_pos in range(pwm_length):
            seq_pos = start_pos + pwm_pos
            nucleotide = sequence[seq_pos].upper()

            if nucleotide not in nucleotides:
                score = 0.0  # Invalid nucleotide
                break

            nuc_idx = nucleotides.index(nucleotide)
            pwm_score = pwm.iloc[pwm_pos, nuc_idx]
            score += pwm_score

        scores.append(score)

    return scores


def motif_conservation_score(sequences: List[str], motif: str) -> float:
    """Calculate conservation score for a motif across sequences.

    Args:
        sequences: List of sequences containing the motif
        motif: Motif pattern to evaluate

    Returns:
        Conservation score (0.0 to 1.0, higher is more conserved)

    Example:
        >>> seqs = ["ATCG", "ATCG", "ATCG"]
        >>> score = motif_conservation_score(seqs, "ATC")
        >>> score
        1.0
    """
    if not sequences or not motif:
        return 0.0

    motif_length = len(motif)
    total_positions = 0
    conserved_positions = 0

    # For each position in the motif
    for pos in range(motif_length):
        nucleotides_at_pos = []

        # Collect nucleotides at this position across all sequences
        for seq in sequences:
            if pos < len(seq):
                nucleotides_at_pos.append(seq[pos].upper())

        if nucleotides_at_pos:
            # Count most common nucleotide
            most_common = max(set(nucleotides_at_pos), key=nucleotides_at_pos.count)
            conserved_count = nucleotides_at_pos.count(most_common)

            total_positions += len(nucleotides_at_pos)
            conserved_positions += conserved_count

    return conserved_positions / total_positions if total_positions > 0 else 0.0


def discover_motifs(sequences: List[str], motif_length: int = 6, min_occurrences: int = 2) -> List[Tuple[str, int]]:
    """Discover overrepresented motifs in a set of sequences.

    This is a simple motif discovery algorithm that finds k-mers that appear
    more frequently than expected by chance.

    Args:
        sequences: List of DNA sequences
        motif_length: Length of motifs to discover (default: 6)
        min_occurrences: Minimum number of occurrences required (default: 2)

    Returns:
        List of tuples (motif, count) sorted by frequency

    Example:
        >>> seqs = ["ATCGATCGATCG", "ATCGATCGATCG", "GCTAGCTAGCTA"]
        >>> motifs = discover_motifs(seqs, motif_length=4, min_occurrences=2)
        >>> len(motifs) >= 0
        True
    """
    if not sequences or motif_length <= 0:
        return []

    from collections import Counter

    # Extract all k-mers
    kmer_counts = Counter()

    for seq in sequences:
        seq = seq.upper()
        for i in range(len(seq) - motif_length + 1):
            kmer = seq[i:i + motif_length]
            # Only consider canonical nucleotides
            if all(c in 'ACGT' for c in kmer):
                kmer_counts[kmer] += 1

    # Filter by minimum occurrences
    filtered_motifs = [
        (motif, count) for motif, count in kmer_counts.items()
        if count >= min_occurrences
    ]

    # Sort by frequency (descending)
    filtered_motifs.sort(key=lambda x: x[1], reverse=True)

    return filtered_motifs


def motif_similarity(motif1: str, motif2: str) -> float:
    """Calculate similarity between two motifs.

    Args:
        motif1: First motif sequence
        motif2: Second motif sequence

    Returns:
        Similarity score (0.0 to 1.0)

    Example:
        >>> motif_similarity("ATCG", "ATCG")
        1.0
        >>> motif_similarity("ATCG", "GCTA")
        0.0
    """
    if not motif1 or not motif2:
        return 0.0

    min_len = min(len(motif1), len(motif2))
    if min_len == 0:
        return 0.0

    matches = 0
    for i in range(min_len):
        if motif1[i].upper() == motif2[i].upper():
            matches += 1

    return matches / min_len


def find_palindromic_motifs(seq: str, min_length: int = 6) -> List[Tuple[str, int]]:
    """Find palindromic motifs in a DNA sequence.

    Args:
        seq: DNA sequence to search
        min_length: Minimum length of palindromic motifs (default: 6)

    Returns:
        List of tuples (motif, start_position)

    Example:
        >>> find_palindromic_motifs("ATCGCGAT", min_length=4)
        [('CGCG', 2)]
    """
    if not seq or min_length < 2:
        return []

    palindromes = []
    seq_upper = seq.upper()

    for i in range(len(seq_upper) - min_length + 1):
        for length in range(min_length, len(seq_upper) - i + 1):
            substring = seq_upper[i:i + length]
            if is_palindrome(substring):
                palindromes.append((substring, i))
                break  # Take the longest palindrome starting at this position

    return palindromes


def is_palindrome(seq: str) -> bool:
    """Check if a DNA sequence is palindromic.

    A palindromic sequence reads the same forwards and backwards
    when complemented.

    Args:
        seq: DNA sequence to check

    Returns:
        True if sequence is palindromic

    Example:
        >>> is_palindrome("ATCG")
        False
        >>> is_palindrome("ATCGCGAT")  # ATCGCGAT -> ATCGCGAT (complement reads same)
        False
    """
    if not seq:
        return True

    seq_upper = seq.upper()
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Check if sequence equals its reverse complement
    reversed_complement = ''.join(complement.get(c, c) for c in reversed(seq_upper))

    return seq_upper == reversed_complement


